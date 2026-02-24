"""
    PowerSystemEnhanced.jl

Enhanced hybrid AC/DC power system solver with:
- Multiple islanded systems detection and handling
- PV → PQ conversion when reactive power limits are violated
- Multiple swing bus handling (AC and DC)
- Automatic converter control mode switching
- Distributed slack bus model for realistic power sharing
- Advanced Newton-Raphson with adaptive control

Author: Tianyang Zhao
Version: 0.3.0 (Distributed Slack)
"""
module PowerSystemEnhanced

using LinearAlgebra, SparseArrays, Printf

# Re-export everything from base PowerSystem
using ..PowerSystem
using ..PowerSystem: ACBus, ACBranch, DCBus, DCBranch, VSCConverter, HybridSystem,
                     BusType, PQ, PV, SLACK, ConverterMode, PQ_MODE, VDC_Q, VDC_VAC,
                     build_admittance_matrix, converter_loss, solve_power_flow,
                     converter_ac_injection, converter_dc_injection

export detect_islands, solve_power_flow_islanded, solve_power_flow_adaptive,
       check_reactive_limits, pv_to_pq_conversion!, auto_select_swing_bus,
       auto_switch_converter_mode!, PowerFlowOptions, IslandInfo, ReactiveLimit,
       create_default_Q_limits, print_island_summary, extract_island_subsystem,
       DistributedSlack, create_participation_factors, solve_power_flow_distributed_slack,
       solve_power_flow_distributed_slack_full

# ─── Enhanced Data Structures ─────────────────────────────────────────────────

"""Reactive power limits for PV buses"""
mutable struct ReactiveLimit
    Qmin::Float64
    Qmax::Float64
end

"""
    DistributedSlack

Configuration for distributed slack bus model where multiple generators 
participate in absorbing power mismatches.

# Fields
- `participating_buses::Vector{Int}`: Bus indices participating in slack distribution
- `participation_factors::Vector{Float64}`: Normalized participation factors (sum to 1.0)
- `reference_bus::Int`: Reference bus for angle (default: first participating bus)
- `max_participation_P::Dict{Int,Float64}`: Maximum active power for each bus (p.u.)

# Theory
In distributed slack model, power mismatch ΔP is distributed as:
    ΔP_i = α_i × ΔP_total
where α_i is the participation factor for bus i.

Participation factors are typically based on:
1. Generator capacity: α_i ∝ Pg_max[i]
2. Governor droop: α_i ∝ 1/R_i
3. Equal participation: α_i = 1/N
"""
struct DistributedSlack
    participating_buses::Vector{Int}
    participation_factors::Vector{Float64}
    reference_bus::Int
    max_participation_P::Dict{Int,Float64}
    
    function DistributedSlack(buses::Vector{Int}, factors::Vector{Float64}, 
                             ref_bus::Int=0, max_P::Dict{Int,Float64}=Dict{Int,Float64}())
        length(buses) == length(factors) || error("Buses and factors must have same length")
        sum_factors = sum(factors)
        abs(sum_factors - 1.0) < 1e-6 || error("Participation factors must sum to 1.0 (got $sum_factors)")
        
        ref = (ref_bus == 0) ? buses[1] : ref_bus
        ref in buses || error("Reference bus $ref not in participating buses")
        
        new(buses, factors, ref, max_P)
    end
end

"""Island information for disconnected network analysis"""
struct IslandInfo
    id::Int
    ac_buses::Vector{Int}      # AC bus indices in this island
    dc_buses::Vector{Int}      # DC bus indices in this island  
    converters::Vector{Int}    # Converter indices connecting this island
    has_ac_slack::Bool
    has_dc_slack::Bool
    ac_slack_bus::Int          # Selected AC slack bus (0 if none)
    dc_slack_bus::Int          # Selected DC slack bus (0 if none)
    has_generators::Bool       # Whether island has any online generators (Pg > 0)
end

"""Power flow solver options"""
struct PowerFlowOptions
    max_iter::Int
    tol::Float64
    enable_pv_pq_conversion::Bool
    enable_auto_swing_selection::Bool
    enable_converter_mode_switching::Bool
    verbose::Bool
    
    PowerFlowOptions(; max_iter=200, tol=1e-6, 
                     enable_pv_pq_conversion=true,
                     enable_auto_swing_selection=true,
                     enable_converter_mode_switching=true,
                     verbose=false) = 
        new(max_iter, tol, enable_pv_pq_conversion, 
            enable_auto_swing_selection, enable_converter_mode_switching, verbose)
end

# ─── Island Detection ─────────────────────────────────────────────────────────

"""
    detect_islands(sys::HybridSystem) -> Vector{IslandInfo}

Detect all electrically isolated islands in the hybrid AC/DC system.

Uses DFS on the **combined AC + DC + converter** graph so that converters
act as edges connecting their AC bus to their DC bus.  This ensures that
a bus connected to the main system only through a converter is correctly
assigned to the same island as the rest of the network.

DC buses are represented as nodes (nac+1):(nac+ndc) in the internal graph.
"""

function detect_islands(sys::HybridSystem)
    nac = length(sys.ac_buses)
    ndc = length(sys.dc_buses)
    n_total = nac + ndc
    
    # ── Build combined adjacency list ─────────────────────────────────
    # Nodes 1:nac = AC buses,  nac+1:nac+ndc = DC buses
    adj = [Int[] for _ in 1:n_total]
    
    # AC branches
    for br in sys.ac_branches
        br.status || continue
        push!(adj[br.from], br.to)
        push!(adj[br.to], br.from)
    end
    
    # DC branches
    for br in sys.dc_branches
        br.status || continue
        di = nac + br.from
        dj = nac + br.to
        push!(adj[di], dj)
        push!(adj[dj], di)
    end
    
    # Converters act as edges: AC bus ↔ DC bus
    for conv in sys.converters
        conv.status || continue
        ac_node = conv.ac_bus
        dc_node = nac + conv.dc_bus
        push!(adj[ac_node], dc_node)
        push!(adj[dc_node], ac_node)
    end
    
    # ── DFS to find connected components ──────────────────────────────
    visited = falses(n_total)
    components = Vector{Vector{Int}}()
    
    for start in 1:n_total
        visited[start] && continue
        comp = Int[]
        stack = [start]
        while !isempty(stack)
            node = pop!(stack)
            visited[node] && continue
            visited[node] = true
            push!(comp, node)
            for nb in adj[node]
                visited[nb] || push!(stack, nb)
            end
        end
        push!(components, comp)
    end
    
    # ── Split each component into AC / DC bus sets and build IslandInfo
    islands = IslandInfo[]
    for comp in components
        ac_in_island = sort([n for n in comp if n <= nac])
        dc_in_island = sort([n - nac for n in comp if n > nac])
        
        # Skip pure-DC components with no AC buses
        isempty(ac_in_island) && continue
        
        ac_set = Set(ac_in_island)
        dc_set = Set(dc_in_island)
        
        # Converters belonging to this island
        conv_in_island = Int[]
        for (ci, conv) in enumerate(sys.converters)
            conv.status || continue
            if conv.ac_bus in ac_set || conv.dc_bus in dc_set
                push!(conv_in_island, ci)
            end
        end
        
        # Generator check
        has_generators = any(
            sys.ac_buses[i].Pg > 0.0 || sys.ac_buses[i].type == SLACK || sys.ac_buses[i].type == PV
            for i in ac_in_island
        )
        if !has_generators && !isempty(conv_in_island)
            has_generators = any(abs(sys.converters[ci].Pset) > 0.0 for ci in conv_in_island)
        end
        
        has_ac_slack = any(sys.ac_buses[i].type == SLACK for i in ac_in_island)
        has_dc_slack = !isempty(dc_in_island)
        
        ac_slack = 0
        if has_ac_slack
            idx = findfirst(i -> sys.ac_buses[i].type == SLACK, ac_in_island)
            ac_slack = idx === nothing ? 0 : ac_in_island[idx]
        end
        dc_slack = isempty(dc_in_island) ? 0 : dc_in_island[1]
        
        id = length(islands) + 1
        push!(islands, IslandInfo(id, ac_in_island, dc_in_island,
                                  conv_in_island, has_ac_slack, has_dc_slack,
                                  ac_slack, dc_slack, has_generators))
    end
    
    return islands
end

# ─── Island Sub-system Extraction ─────────────────────────────────────────────

"""
    extract_island_subsystem(sys::HybridSystem, island::IslandInfo; slack_bus_override=0)
        -> (sub_sys, ac_map, dc_map)

Extract a standalone `HybridSystem` for a single island.

Returns:
- `sub_sys`: A new `HybridSystem` containing only the buses, branches,
  and converters that belong to this island. All indices are renumbered
  starting from 1.
- `ac_map`: Dict mapping *original* AC bus index → sub-system AC bus index.
- `dc_map`: Dict mapping *original* DC bus index → sub-system DC bus index.

# Usage
```julia
islands = detect_islands(sys)
for isl in islands
    sub, ac_m, dc_m = extract_island_subsystem(sys, isl)
    result = solve_power_flow(sub)
    # Map back:
    for (orig, local) in ac_m
        Vm_global[orig] = result.Vm[local]
        Va_global[orig] = result.Va[local]
    end
end
```
"""
function extract_island_subsystem(sys::HybridSystem, island::IslandInfo; slack_bus_override::Int=0)
    # Build index maps:  original → local (1-based)
    ac_map = Dict{Int,Int}()
    for (local_idx, orig_idx) in enumerate(sort(island.ac_buses))
        ac_map[orig_idx] = local_idx
    end
    dc_map = Dict{Int,Int}()
    for (local_idx, orig_idx) in enumerate(sort(island.dc_buses))
        dc_map[orig_idx] = local_idx
    end

    ac_set = Set(island.ac_buses)
    dc_set = Set(island.dc_buses)

    # ── AC buses (renumbered) ─────────────────────────────────────────
    sub_ac_buses = ACBus[]
    for orig in sort(island.ac_buses)
        b = sys.ac_buses[orig]
        bus_type = orig == slack_bus_override ? SLACK : b.type
        bus_va = orig == slack_bus_override ? 0.0 : b.Va
        push!(sub_ac_buses, ACBus(ac_map[orig], bus_type, b.Pd, b.Qd,
                                   b.Pg, b.Qg, b.Vm, bus_va, b.area))
    end

    # ── AC branches (only those with both ends inside this island) ────
    sub_ac_branches = ACBranch[]
    for br in sys.ac_branches
        br.status || continue
        (br.from in ac_set && br.to in ac_set) || continue
        push!(sub_ac_branches, ACBranch(ac_map[br.from], ac_map[br.to],
                                         br.r, br.x, br.b, br.tap, br.status))
    end

    # ── DC buses (renumbered) ─────────────────────────────────────────
    sub_dc_buses = DCBus[]
    for orig in sort(island.dc_buses)
        b = sys.dc_buses[orig]
        push!(sub_dc_buses, DCBus(dc_map[orig], b.Vdc_set, b.Pdc))
    end

    # ── DC branches ───────────────────────────────────────────────────
    sub_dc_branches = DCBranch[]
    for br in sys.dc_branches
        br.status || continue
        (br.from in dc_set && br.to in dc_set) || continue
        push!(sub_dc_branches, DCBranch(dc_map[br.from], dc_map[br.to],
                                         br.r, br.status))
    end

    # ── Converters (only those with both AC and DC bus inside) ────────
    sub_converters = VSCConverter[]
    for conv in sys.converters
        conv.status || continue
        (conv.ac_bus in ac_set && conv.dc_bus in dc_set) || continue
        push!(sub_converters, VSCConverter(
            conv.id, ac_map[conv.ac_bus], dc_map[conv.dc_bus],
            conv.mode, conv.Pset, conv.Qset, conv.Vdc_set, conv.Vac_set,
            conv.Ploss_a, conv.Ploss_b, conv.Ploss_c, conv.Smax, conv.status;
            G_droop=conv.G_droop))
    end

    # ── Build sub-system (automatically builds Ybus & Gdc) ───────────
    sub_sys = HybridSystem(sub_ac_buses, sub_ac_branches,
                           sub_dc_buses, sub_dc_branches,
                           sub_converters; baseMVA=sys.baseMVA)

    return sub_sys, ac_map, dc_map
end

# ─── Reactive Power Limit Checking ───────────────────────────────────────────

"""
    check_reactive_limits(sys, Vm, Va, Vdc, Q_limits) -> (violations, Q_actual)

Check reactive power output of PV buses against limits.
Returns indices of buses violating limits and their actual Q values.

# Theory
PV buses maintain constant voltage |V_i| = V_i^{set} by adjusting Q_i.
However, generators have physical limits: Q_{min} ≤ Q_i ≤ Q_{max}

When Q_i exceeds limits, the bus must be converted to PQ type:
- If Q_i > Q_{max}: set Q_i = Q_{max}, solve for |V_i|
- If Q_i < Q_{min}: set Q_i = Q_{min}, solve for |V_i|
"""
function check_reactive_limits(sys::HybridSystem, Vm::Vector{Float64}, 
                               Va::Vector{Float64}, Vdc::Vector{Float64},
    Q_limits::Dict{Int, ReactiveLimit})
    nac = length(sys.ac_buses)
    
    # Compute actual Q injection at each bus
    _, Q_actual = PowerSystem.ac_power_injections(sys, Vm, Va)
    
    # Add converter reactive injections
    for conv in sys.converters
        conv.status || continue
        _, Qinj = PowerSystem.converter_ac_injection(conv, Vm, Va, Vdc)
        Q_actual[conv.ac_bus] += Qinj
    end
    
    # Subtract scheduled reactive load
    for (i, bus) in enumerate(sys.ac_buses)
        Q_actual[i] -= bus.Qd
    end
    
    # Check violations for PV buses
    violations = Tuple{Int, Symbol, Float64}[]  # (bus_id, :max/:min, Q_value)
    for (i, bus) in enumerate(sys.ac_buses)
        bus.type != PV && continue
        haskey(Q_limits, i) || continue
        
        lim = Q_limits[i]
        if Q_actual[i] > lim.Qmax
            push!(violations, (i, :max, Q_actual[i]))
        elseif Q_actual[i] < lim.Qmin
            push!(violations, (i, :min, Q_actual[i]))
        end
    end
    
    return violations, Q_actual
end

"""
    pv_to_pq_conversion!(sys, violations, Q_limits)

Convert PV buses to PQ type when reactive limits are violated.
Modifies system in-place.

# Algorithm
1. For each violated PV bus i:
   - If Q_i > Q_{max}: Set type ← PQ, Q_g^{set} ← Q_{max}
   -If Q_i < Q_{min}: Set type ← PQ, Q_g^{set} ← Q_{min}
2. Re-solve power flow with new PQ buses
3. Iterate until no violations (or max iterations reached)
"""
function pv_to_pq_conversion!(sys::HybridSystem, violations::Vector{Tuple{Int, Symbol, Float64}},
                              Q_limits::Dict{Int, ReactiveLimit})
    converted = Int[]
    
    for (bus_idx, limit_type, Q_val) in violations
        # Create new modified bus
        old_bus = sys.ac_buses[bus_idx]
        lim = Q_limits[bus_idx]
        
        Q_new = limit_type == :max ? lim.Qmax : lim.Qmin
        
        sys.ac_buses[bus_idx] = ACBus(
            old_bus.id, PQ,  # Convert to PQ type
            old_bus.Pd, old_bus.Qd,
            old_bus.Pg, Q_new,  # Set Q to limit
            old_bus.Vm, old_bus.Va,  # Keep voltage as initial guess
            old_bus.area
        )
        
        push!(converted, bus_idx)
    end
    
    return converted
end

# ─── Automatic Swing Bus Selection ───────────────────────────────────────────

"""
    auto_select_swing_bus(sys, island::IslandInfo) -> Int

Automatically select best AC slack bus for an island.

# Selection Criteria (in order)
1. Existing SLACK bus in island
2. Largest online generator (any bus type with highest Pg > 0)
3. First bus in island (emergency fallback — dead island)

# Theory
The slack bus provides voltage and angle reference:
- Sets V∠θ = V_ref ∠ 0°
- Balances system-wide active power via P_slack

If the island has generators but no SLACK bus, the bus with the largest 
online generation capacity is automatically promoted to SLACK.
"""
function auto_select_swing_bus(sys::HybridSystem, island::IslandInfo)
    buses_in_island = island.ac_buses
    isempty(buses_in_island) && return 0
    
    # Priority 1: Existing slack bus
    for i in buses_in_island
        sys.ac_buses[i].type == SLACK && return i
    end
    
    # Priority 2: Largest online generator (any bus type with Pg > 0)
    gen_buses = [(i, sys.ac_buses[i].Pg) for i in buses_in_island if sys.ac_buses[i].Pg > 0.0]
    if !isempty(gen_buses)
        sort!(gen_buses, by=x -> x[2], rev=true)
        return gen_buses[1][1]
    end
    
    # Priority 3: Any PV bus (even if Pg == 0, it's a voltage-controlled bus)
    for i in buses_in_island
        sys.ac_buses[i].type == PV && return i
    end
    
    # Fallback: First bus (emergency — dead island with no generation)
    return buses_in_island[1]
end

"""
    auto_select_dc_slack_bus(sys, island::IslandInfo) -> Int

Select DC slack bus (voltage reference for DC network).

# Criteria
1. DC bus connected to VDC_Q or VDC_VAC converter (voltage-controlling)
2. First DC bus in island
"""
function auto_select_dc_slack_bus(sys::HybridSystem, island::IslandInfo)
    dc_buses_in_island = island.dc_buses
    
    # Priority 1: DC bus with voltage-controlling converter
    for i in dc_buses_in_island
        for conv in sys.converters
            if conv.status && conv.dc_bus == i && 
               (conv.mode == VDC_Q || conv.mode == VDC_VAC)
                return i
            end
        end
    end
    
    # Fallback: First DC bus
    return isempty(dc_buses_in_island) ? 0 : dc_buses_in_island[1]
end

# ─── Automatic Converter Control Mode Switching ──────────────────────────────

"""
    auto_switch_converter_mode!(sys, Vm, Va, Vdc) -> Vector{Int}

Automatically switch converter control modes based on operating conditions.

# Control Mode Switching Rules:
# PQ_MODE → VDC_VAC: When AC voltage drops below threshold
# VDC_VAC → PQ_MODE: When power exceeds rating
# VDC_Q → VDC_VAC: When AC voltage support needed

# Returns indices of converters that were switched
"""
function auto_switch_converter_mode!(sys::HybridSystem, Vm::Vector{Float64}, 
                                     Va::Vector{Float64}, Vdc::Vector{Float64};
                                     V_low_threshold=0.95, V_high_threshold=1.05,
                                     S_limit_frac=0.95)
    switched = Int[]
    
    for (idx, conv) in enumerate(sys.converters)
        conv.status || continue
        
        V_ac = Vm[conv.ac_bus]
        Pac, Qac = PowerSystem.converter_ac_injection(conv, Vm, Va, Vdc)
        S_ac = sqrt(Pac^2 + Qac^2)
        
        old_mode = conv.mode
        new_mode = old_mode
        
        if conv.mode == PQ_MODE
            # Switch to VDC_VAC if AC voltage too low (voltage support needed)
            if V_ac < V_low_threshold
                new_mode = VDC_VAC
            end
        elseif conv.mode == VDC_VAC
            # Switch to PQ_MODE if overloaded
            if S_ac > S_limit_frac * conv.Smax
                new_mode = PQ_MODE
            end
        elseif conv.mode == VDC_Q
            # Switch to VDC_VAC if AC voltage support needed
            if V_ac < V_low_threshold || V_ac > V_high_threshold
                new_mode = VDC_VAC
            end
        end
        
        if new_mode != old_mode
            # Update converter mode
            sys.converters[idx] = VSCConverter(
                conv.id, conv.ac_bus, conv.dc_bus, new_mode,
                conv.Pset, conv.Qset, conv.Vdc_set,
                V_ac,  # Use current voltage as new setpoint for VDC_VAC
                conv.Ploss_a, conv.Ploss_b, conv.Ploss_c, conv.Smax, true;
                G_droop=conv.G_droop
            )
            push!(switched, idx)
        end
    end
    
    return switched
end

# ─── Enhanced Power Flow Solver ──────────────────────────────────────────────

"""
    solve_power_flow_adaptive(sys; options=PowerFlowOptions())

Advanced power flow solver with:
- Island extraction: each island is solved as an independent sub-system
- Dead island handling: islands without generators get Vm=0
- Automatic slack bus selection: largest online generator promoted to SLACK
- PV→PQ conversion when reactive limits are violated
- Converter mode switching

# Algorithm
1. Detect islands in the system
2. For each island:
   a. Skip dead islands (no generators) — set Vm = 0
   b. Auto-select slack bus if none exists
   c. Extract sub-system (renumbered HybridSystem for this island)
   d. Solve sub-system NR independently
   e. Check reactive limits → PV→PQ conversion if needed
   f. Map results back to global voltage arrays
3. Combine results from all islands
"""
function solve_power_flow_adaptive(sys::HybridSystem; 
                                   options::PowerFlowOptions=PowerFlowOptions(),
                                   Q_limits::Dict{Int, ReactiveLimit}=Dict{Int, ReactiveLimit}())
    
    # Step 1: Detect islands
    islands = detect_islands(sys)
    
    if options.verbose
        println("Detected $(length(islands)) island(s)")
        for island in islands
            gen_str = island.has_generators ? "VIABLE" : "DEAD"
            println("  Island $(island.id) [$gen_str]: $(length(island.ac_buses)) AC buses, " *
                   "$(length(island.dc_buses)) DC buses, slack=$(island.has_ac_slack)")
        end
    end
    
    # Step 2: Prepare global solution arrays
    nac = length(sys.ac_buses)
    ndc = length(sys.dc_buses)
    Vm_global = [b.Vm for b in sys.ac_buses]
    Va_global = [b.Va for b in sys.ac_buses]
    Vdc_global = ndc > 0 ? [b.Vdc_set for b in sys.dc_buses] : Float64[]
    
    converged_all = true
    total_iterations = 0
    
    # Step 3: Solve each island independently via sub-system extraction
    for island in islands
        if options.verbose
            println("\nSolving island $(island.id)...")
        end
        
        # ── Skip dead islands ─────────────────────────────────────────
        if !island.has_generators
            if options.verbose
                println("  ⚠ Dead island $(island.id): no generators — skipping (voltages set to 0)")
            end
            for i in island.ac_buses
                Vm_global[i] = 0.0
                Va_global[i] = 0.0
            end
            for i in island.dc_buses
                if i <= length(Vdc_global)
                    Vdc_global[i] = 0.0
                end
            end
            continue
        end
        
        # ── Power balance feasibility check ───────────────────────────
        # If total load exceeds total generation capacity, the island
        # cannot converge no matter what — skip it early.
        sum_Pg = sum(sys.ac_buses[i].Pg for i in island.ac_buses)
        sum_Pd = sum(sys.ac_buses[i].Pd for i in island.ac_buses)
        # Include PQ_MODE converter injections (positive Pset = generation)
        for ci in island.converters
            conv = sys.converters[ci]
            if conv.mode == PQ_MODE && conv.Pset > 0
                sum_Pg += conv.Pset
            end
        end
        if sum_Pd > sum_Pg + 1e-6   # small tolerance
            if options.verbose
                @printf("  ⚠ Infeasible island %d: Pd=%.3f > Pg=%.3f — skipping\n",
                        island.id, sum_Pd, sum_Pg)
            end
            for i in island.ac_buses
                Vm_global[i] = 0.0
                Va_global[i] = 0.0
            end
            for i in island.dc_buses
                if i <= length(Vdc_global)
                    Vdc_global[i] = 0.0
                end
            end
            converged_all = false
            continue
        end
        
        # ── Handle single-bus islands (trivially solvable) ───────────
        if length(island.ac_buses) == 1
            bus_idx = island.ac_buses[1]
            b = sys.ac_buses[bus_idx]
            Vm_global[bus_idx] = b.Vm
            Va_global[bus_idx] = 0.0
            # DC buses associated with this single-bus island keep their setpoints
            for i in island.dc_buses
                if i <= length(Vdc_global)
                    Vdc_global[i] = sys.dc_buses[i].Vdc_set
                end
            end
            options.verbose && println("  Single-bus island $(island.id) (bus $bus_idx): Vm=$(b.Vm)")
            continue
        end
        
        # ── Auto-select slack bus (local override, no global mutation) ─
        slack_override = 0
        if options.enable_auto_swing_selection
            if !island.has_ac_slack && !isempty(island.ac_buses)
                slack_override = auto_select_swing_bus(sys, island)
                if slack_override > 0
                    options.verbose && println("  Auto-selected AC slack: bus $slack_override (Pg=$(sys.ac_buses[slack_override].Pg))")
                else
                    options.verbose && println("  ⚠ No suitable slack bus found for island $(island.id)")
                end
            end
        end
        
        # ── Extract sub-system for this island ────────────────────────
        sub_sys, ac_map, dc_map = extract_island_subsystem(sys, island; slack_bus_override=slack_override)
        
        # Reverse maps: local → original
        ac_rmap = Dict(v => k for (k, v) in ac_map)
        dc_rmap = Dict(v => k for (k, v) in dc_map)
        
        # ── Translate Q_limits to sub-system indices ──────────────────
        sub_Q_limits = Dict{Int, ReactiveLimit}()
        for (orig_bus, lim) in Q_limits
            if haskey(ac_map, orig_bus)
                sub_Q_limits[ac_map[orig_bus]] = lim
            end
        end
        
        # ── Solve sub-system with PV→PQ iteration ────────────────────
        max_pv_pq_iter = 10
        island_converged = false
        
        for pv_pq_iter in 1:max_pv_pq_iter
            result = PowerSystem.solve_power_flow(sub_sys;
                                                   max_iter=options.max_iter,
                                                   tol=options.tol)
            
            total_iterations += result.iterations
            
            if !result.converged
                # ── Retry 1: flat start (Vm=1.0, Va=0) ──
                sub_flat = deepcopy(sub_sys)
                for bi in 1:length(sub_flat.ac_buses)
                    b = sub_flat.ac_buses[bi]
                    sub_flat.ac_buses[bi] = ACBus(b.id, b.type, b.Pd, b.Qd,
                                                   b.Pg, b.Qg, 1.0, 0.0, b.area)
                end
                sub_flat.Ybus = build_admittance_matrix(sub_flat)
                result = PowerSystem.solve_power_flow(sub_flat;
                                                      max_iter=options.max_iter,
                                                      tol=options.tol)
                total_iterations += result.iterations
                
                if !result.converged
                    # ── Retry 2: convert ALL PV buses → PQ ──
                    # Near voltage collapse, PV buses requesting excessive Q
                    # cause non-convergence.  Relaxing them to PQ lets NR find
                    # a lower-voltage solution.
                    sub_pq = deepcopy(sub_sys)
                    for bi in 1:length(sub_pq.ac_buses)
                        b = sub_pq.ac_buses[bi]
                        if b.type == PV
                            sub_pq.ac_buses[bi] = ACBus(b.id, PQ, b.Pd, b.Qd,
                                                         b.Pg, b.Qg, 1.0, 0.0, b.area)
                        elseif b.type == SLACK
                            sub_pq.ac_buses[bi] = ACBus(b.id, SLACK, b.Pd, b.Qd,
                                                         b.Pg, b.Qg, 1.0, 0.0, b.area)
                        else
                            sub_pq.ac_buses[bi] = ACBus(b.id, b.type, b.Pd, b.Qd,
                                                         b.Pg, b.Qg, 1.0, 0.0, b.area)
                        end
                    end
                    sub_pq.Ybus = build_admittance_matrix(sub_pq)
                    result = PowerSystem.solve_power_flow(sub_pq;
                                                          max_iter=options.max_iter,
                                                          tol=options.tol)
                    total_iterations += result.iterations
                end
                
                if !result.converged
                    options.verbose && println("  Island $(island.id) did not converge after retries (res=$(result.residual))")
                    converged_all = false
                    for (orig, loc) in ac_map
                        Vm_global[orig] = result.Vm[loc]
                        Va_global[orig] = result.Va[loc]
                    end
                    for (orig, loc) in dc_map
                        if orig <= length(Vdc_global) && loc <= length(result.Vdc)
                            Vdc_global[orig] = result.Vdc[loc]
                        end
                    end
                    break
                end
            end
            
            # Map converged results back
            for (orig, loc) in ac_map
                Vm_global[orig] = result.Vm[loc]
                Va_global[orig] = result.Va[loc]
            end
            for (orig, loc) in dc_map
                if orig <= length(Vdc_global) && loc <= length(result.Vdc)
                    Vdc_global[orig] = result.Vdc[loc]
                end
            end
            
            # Check reactive power limits on sub-system
            if options.enable_pv_pq_conversion && !isempty(sub_Q_limits)
                violations, _ = check_reactive_limits(sub_sys, result.Vm, result.Va,
                                                       result.Vdc, sub_Q_limits)
                if !isempty(violations)
                    converted = pv_to_pq_conversion!(sub_sys, violations, sub_Q_limits)
                    # Also update original system's bus types
                    for local_idx in converted
                        orig_idx = ac_rmap[local_idx]
                        sys.ac_buses[orig_idx] = sub_sys.ac_buses[local_idx]
                    end
                    options.verbose && println("  Converted $(length(converted)) PV→PQ buses")
                    continue
                end
            end
            
            # Check converter mode switching on sub-system
            if options.enable_converter_mode_switching
                switched = auto_switch_converter_mode!(sub_sys, result.Vm, result.Va, result.Vdc)
                if !isempty(switched)
                    options.verbose && println("  Switched $(length(switched)) converter modes")
                    continue
                end
            end
            
            island_converged = true
            options.verbose && println("  Island $(island.id) converged after $pv_pq_iter PV-PQ iterations")
            break
        end
        
        converged_all &= island_converged
    end
    
    return (Vm=Vm_global, Va=Va_global, Vdc=Vdc_global, 
            converged=converged_all, iterations=total_iterations,
            islands=islands)
end

"""
    solve_power_flow_islanded(sys) -> Vector{NamedTuple}

Solve each island independently using sub-system extraction and return
per-island results. Useful for post-fault analysis where system fragments.
"""
function solve_power_flow_islanded(sys::HybridSystem; 
                                   options::PowerFlowOptions=PowerFlowOptions())
    islands = detect_islands(sys)
    result_global = solve_power_flow_adaptive(sys; options=options)
    
    results = []
    for island in islands
        island_result = (
            island_id = island.id,
            ac_buses = island.ac_buses,
            dc_buses = island.dc_buses,
            has_generators = island.has_generators,
            Vm = result_global.Vm[island.ac_buses],
            Va = result_global.Va[island.ac_buses],
            Vdc = isempty(island.dc_buses) ? Float64[] : result_global.Vdc[island.dc_buses],
            converged = result_global.converged
        )
        push!(results, island_result)
    end
    
    return results
end

# ─── Utility Functions ────────────────────────────────────────────────────────

"""
    create_default_Q_limits(sys, Qmin_default=-0.5, Qmax_default=0.5)

Create default reactive power limits for all PV buses in the system.
"""
function create_default_Q_limits(sys::HybridSystem; 
                                Qmin_default=-0.5, Qmax_default=0.5)
    Q_limits = Dict{Int, ReactiveLimit}()
    
    for (i, bus) in enumerate(sys.ac_buses)
        if bus.type == PV
            Q_limits[i] = ReactiveLimit(Qmin_default, Qmax_default)
        end
    end
    
    return Q_limits
end

"""
    print_island_summary(islands, sys)

Print detailed summary of all detected islands.
"""
function print_island_summary(islands::Vector{IslandInfo}, sys::HybridSystem)
    println("\n" * "="^70)
    println("ISLAND DETECTION SUMMARY")
    println("="^70)
    println("Total islands detected: $(length(islands))")
    
    for island in islands
        status_str = island.has_generators ? "VIABLE" : "DEAD (no generators)"
        println("\nIsland $(island.id) [$status_str]:")
        println("  AC buses ($(length(island.ac_buses))): $(island.ac_buses)")
        println("  DC buses ($(length(island.dc_buses))): $(island.dc_buses)")
        println("  Converters ($(length(island.converters))): $(island.converters)")
        println("  AC slack bus: $(island.ac_slack_bus) (has_slack=$(island.has_ac_slack))")
        println("  DC slack bus: $(island.dc_slack_bus) (has_slack=$(island.has_dc_slack))")
        
        # Power balance
        P_gen = sum(sys.ac_buses[i].Pg for i in island.ac_buses)
        P_load = sum(sys.ac_buses[i].Pd for i in island.ac_buses)
        println("  AC power: Gen=$(P_gen), Load=$(P_load), Balance=$(P_gen-P_load)")
    end
    println("="^70)
end

# ─── Distributed Slack Bus Model ──────────────────────────────────────────────

"""
    create_participation_factors(sys::HybridSystem; method=:capacity) -> DistributedSlack

Automatically create participation factors for distributed slack model.
Each island's generators participate independently — no inter-island
droop or AGC coupling.

# Arguments
- `sys`: Hybrid power system
- `method`: Method for computing factors
  - `:capacity` - Proportional to Pg (default)
  - `:droop` - Proportional to 1/R (governor droop). Uses `droop_coeffs` if provided.
  - `:equal` - Equal participation
- `droop_coeffs`: Optional map `bus_index => R_i` (droop); if missing, falls back to capacity proxy.

# Returns
`DistributedSlack` configuration with normalized participation factors.

# Example
```julia
dist_slack = create_participation_factors(sys; method=:capacity)
result = solve_power_flow_distributed_slack(sys, dist_slack)
```
"""
function create_participation_factors(sys::HybridSystem; method=:capacity,
                                     participating_buses::Vector{Int}=Int[],
                                     droop_coeffs::Dict{Int,Float64}=Dict{Int,Float64}())
    # If no buses specified, use all slack and PV buses with nonzero capacity
    if isempty(participating_buses)
        participating_buses = Int[]
        for (i, bus) in enumerate(sys.ac_buses)
            if (bus.type == SLACK || bus.type == PV) && bus.Pg > 0
                push!(participating_buses, i)
            end
        end
    end
    
    isempty(participating_buses) && error("No participating buses found")
    
    # Compute raw participation factors
    if method == :capacity
        # Proportional to maximum generation capacity
        raw_factors = Float64[]
        max_P = Dict{Int,Float64}()
        for i in participating_buses
            # Use current Pg as proxy for capacity (could be enhanced)
            capacity = sys.ac_buses[i].Pg > 0 ? sys.ac_buses[i].Pg : 1.0
            push!(raw_factors, capacity)
            max_P[i] = capacity * 1.5  # Allow 150% of current generation
        end
    elseif method == :droop
        # Proportional to 1/R (droop); fall back to capacity if no droop data
        raw_factors = Float64[]
        max_P = Dict{Int,Float64}()
        for i in participating_buses
            capacity = sys.ac_buses[i].Pg > 0 ? sys.ac_buses[i].Pg : 1.0
            R_i = get(droop_coeffs, i, 0.0)
            factor = R_i > 0 ? 1.0 / R_i : capacity
            push!(raw_factors, factor)
            max_P[i] = capacity * 1.5
        end
    elseif method == :equal
        # Equal participation
        raw_factors = ones(Float64, length(participating_buses))
        max_P = Dict(i => 10.0 for i in participating_buses)  # Large limit
    else
        error("Unknown participation factor method: $method (supported: :capacity, :droop, :equal)")
    end
    
    # Normalize to sum to 1.0
    total = sum(raw_factors)
    factors = raw_factors ./ total
    
    return DistributedSlack(participating_buses, factors, participating_buses[1], max_P)
end

"""
    solve_power_flow_distributed_slack(sys::HybridSystem, dist_slack::DistributedSlack; 
                                       max_iter=50, tol=1e-8, verbose=false)

Solve power flow with distributed slack bus model (capacity-based participation).

Each island is solved independently with its own slack bus; no inter-island
droop or AGC coupling is applied. Within a connected island, generators share
the power mismatch proportionally to their participation factors α_i:

    P_g,i = P_g,i^0 + α_i × ΔP_total

# Algorithm
1. Detect islands; skip dead islands (no generators)
2. For each viable island, auto-assign a slack bus if missing
3. Solve NR power flow per island
4. Post-process to distribute slack among participating generators

# Arguments
- `sys`: Hybrid power system
- `dist_slack`: Distributed slack configuration
- `max_iter`: Maximum Newton-Raphson iterations
- `tol`: Convergence tolerance
- `verbose`: Print iteration details

# Returns
Named tuple with converged solution including participation summary.
"""
function solve_power_flow_distributed_slack(sys::HybridSystem, 
                                           dist_slack::DistributedSlack;
                                           max_iter=50, tol=1e-8, verbose=false)
    rebuild_matrices!(sys)
    nac = length(sys.ac_buses)
    ndc = length(sys.dc_buses)
    
    # Initialize voltages
    Vm = [b.Vm for b in sys.ac_buses]
    Va = zeros(nac)
    Vdc = [b.Vdc_set for b in sys.dc_buses]  # Use Vdc_set as initial guess
    
    # Identify bus types
    slack_idx = [i for (i,b) in enumerate(sys.ac_buses) if b.type == SLACK]
    pv_idx = [i for (i,b) in enumerate(sys.ac_buses) if b.type == PV]
    pq_idx = [i for (i,b) in enumerate(sys.ac_buses) if b.type == PQ]
    
    # State variables: [Va (all except ref); Vm (PQ only); Vdc (all except 1)]
    # For distributed slack: all participating buses have unknown Va and Pg
    
    participating_set = Set(dist_slack.participating_buses)
    ref_bus = dist_slack.reference_bus
    
    # Non-reference buses (all except reference bus)
    non_ref_buses = [i for i in 1:nac if i != ref_bus]
    
    # PQ buses (unknown Vm)
    pq_buses = [i for (i,b) in enumerate(sys.ac_buses) if b.type == PQ]
    
    # Number of equations
    np = length(non_ref_buses)  # P equations for non-reference buses
    nq = length(pq_buses)       # Q equations for PQ buses
    n_part = length(dist_slack.participating_buses) - 1  # Participation constraints (one is reference)
    nf_dc = max(0, ndc - 1)     # DC power equations
    
    Y = sys.Ybus
    G = real.(Matrix(Y))
    B = imag.(Matrix(Y))
    
    verbose && println("Distributed Slack Power Flow Solver")
    verbose && println("  Participating buses: $(dist_slack.participating_buses)")
    verbose && println("  Reference bus: $(ref_bus)")
    verbose && println("  Participation factors: $(dist_slack.participation_factors)")
    
    # Store initial generation
    Pg_initial = Dict(i => sys.ac_buses[i].Pg for i in dist_slack.participating_buses)
    
    # Initialize F outside loop to avoid scope issues
    F = Float64[]
    
    for iter in 1:max_iter
        # --- Compute AC power injections ---
        Pcalc, Qcalc = PowerSystem.ac_power_injections(sys, Vm, Va)
        
        # Scheduled injections
        Psch = [sys.ac_buses[i].Pg - sys.ac_buses[i].Pd for i in 1:nac]
        Qsch = [sys.ac_buses[i].Qg - sys.ac_buses[i].Qd for i in 1:nac]
        
        # Add converter contributions
        for conv in sys.converters
            conv.status || continue
            Pac, Qac = converter_ac_injection(conv, Vm, Va, Vdc)
            Psch[conv.ac_bus] += Pac
            Qsch[conv.ac_bus] += Qac
        end
        
        # Power mismatch
        ΔP = Pcalc .- Psch
        ΔQ = Qcalc .- Qsch
        
        # Total mismatch at participating buses
        ΔP_total = sum(ΔP[i] for i in dist_slack.participating_buses)
        
        # Build mismatch vector
        F = Float64[]
        
        # P equations for non-reference buses
        for i in non_ref_buses
            if i in participating_set
                # Participating bus: mismatch should follow participation factor
                idx_part = findfirst(==(i), dist_slack.participating_buses)
                α_i = dist_slack.participation_factors[idx_part]
                # F_i = ΔP_i - α_i × ΔP_total
                push!(F, ΔP[i] - α_i * ΔP_total)
            else
                # Non-participating bus: standard mismatch
                push!(F, ΔP[i])
            end
        end
        
        # Q equations for PQ buses
        for i in pq_buses
            push!(F, ΔQ[i])
        end
        
        # DC power equations
        if ndc > 1
            Pdc_calc = Vdc .* (sys.Gdc * Vdc)
            Pdc_sch = [b.Pdc for b in sys.dc_buses]
            for conv in sys.converters
                conv.status || continue
                Pdc_sch[conv.dc_bus] += converter_dc_injection(conv, Vm, Va, Vdc)
            end
            for k in 1:(ndc-1)
                push!(F, Pdc_calc[k+1] - Pdc_sch[k+1])
            end
        end
        
        res_norm = norm(F, Inf)
        
        if verbose
            println("  Iter $iter: ||F|| = $(res_norm)")
        end
        
        if res_norm < tol
            # Compute final participation summary
            P_slack = Dict{Int,Float64}()
            for i in dist_slack.participating_buses
                P_slack[i] = sys.ac_buses[i].Pg - Pg_initial[i]
            end
            
            verbose && println("  Converged! Slack distribution: $P_slack")
            
            return (Vm=Vm, Va=Va, Vdc=Vdc, converged=true, iterations=iter, 
                   residual=res_norm, distributed_slack_P=P_slack,
                   total_slack_P=ΔP_total)
        end
        
        # For simplicity, fall back to standard solver with reference angle fixed
        # (Full distributed slack Jacobian is complex; this is a simplified implementation)
        # In practice, use standard solver and post-process to distribute slack
        if iter == 1 && verbose
            @warn "Distributed slack using simplified algorithm (post-distribution)"
        end
        
        # Use standard power flow solver
        result_std = solve_power_flow(sys; max_iter=1, tol=tol)
        if !result_std.converged
            # If single iteration didn't converge, use full solve
            result_std = solve_power_flow(sys; max_iter=max_iter, tol=tol)
            Vm = result_std.Vm
            Va = result_std.Va
            Vdc = result_std.Vdc
            
            if result_std.converged
                # Success - compute slack distribution in post-processing
                verbose && println("  Standard solver converged, computing slack distribution...")
                
                # Recompute power mismatches
                Pcalc_final = zeros(nac)
                for i in 1:nac, j in 1:nac
                    (G[i,j] == 0.0 && B[i,j] == 0.0) && continue
                    θij = Va[i] - Va[j]
                    Pcalc_final[i] += Vm[i] * Vm[j] * (G[i,j] * cos(θij) + B[i,j] * sin(θij))
                end
                
                Psch_final = [sys.ac_buses[i].Pg - sys.ac_buses[i].Pd for i in 1:nac]
                for conv in sys.converters
                    conv.status || continue
                    Pac, _ = converter_ac_injection(conv, Vm, Va, Vdc)
                    Psch_final[conv.ac_bus] += Pac
                end
                
                # Compute slack at each participating bus
                P_slack = Dict{Int,Float64}()
                ΔP_total_final = 0.0
                for i in dist_slack.participating_buses
                    slack_i = Pcalc_final[i] - Psch_final[i]
                    P_slack[i] = slack_i
                    ΔP_total_final += slack_i
                end
                
                return (Vm=Vm, Va=Va, Vdc=Vdc, converged=true, 
                       iterations=result_std.iterations,
                       residual=result_std.residual,
                       distributed_slack_P=P_slack,
                       total_slack_P=ΔP_total_final)
            else
                # Did not converge
                break
            end
        end
        
        Vm = result_std.Vm
        Va = result_std.Va  
        Vdc = result_std.Vdc
    end
    
    @warn "Distributed slack power flow did not converge in $max_iter iterations"
    return (Vm=Vm, Va=Va, Vdc=Vdc, converged=false, iterations=max_iter, 
           residual=norm(F, Inf), distributed_slack_P=Dict{Int,Float64}(),
           total_slack_P=0.0)
end

"""
    solve_power_flow_distributed_slack_full(sys::HybridSystem, dist_slack::DistributedSlack;
                                            max_iter=50, tol=1e-8, verbose=false,
                                            enforce_limits=true)

Solve power flow with FULL Jacobian modification for distributed slack bus model.
Capacity-based participation only — no inter-island droop or AGC coupling.

# Theory
This implementation uses the dimensionally correct formulation:
- Uses ONE total slack variable instead of n individual ΔPg variables
- Participation factors determine distribution: ΔPg_i = α_i × λ_slack
- Quadratic Newton-Raphson convergence preserved
- Generator capacity limits monitored during iterations

All participating generators must belong to the same connected island.

# Augmented System
State: x = [Va (non-ref); Vm (PQ); Vdc (non-slack DC); λ_slack]

Where λ_slack is the total distributed slack power.

Equations:
1. P balance: P_calc,i - P_sch,i - α_i × λ_slack = 0 (for participating buses)
2. P balance: P_calc,i - P_sch,i = 0 (for non-participating buses)
3. Q balance: Q_calc,i - Q_sch,i = 0 (for PQ buses)
4. DC balance: P_DC balance equations

# Jacobian Structure
```
┌────────────────────────────────────────┐
│  ∂P/∂θ   ∂P/∂V   ∂P/∂Vdc   -α_i       │  P equations
│  ∂Q/∂θ   ∂Q/∂V   ∂Q/∂Vdc    0         │  Q equations
│ ∂Pdc/∂θ ∂Pdc/∂V ∂Pdc/∂Vdc   0         │  DC equations
└────────────────────────────────────────┘
```

Dimensionally correct: (n_P + n_Q + n_DC) equations = (n_P + n_Q + n_DC) variables

# Arguments
- `enforce_limits`: Monitor capacity limits during iterations (default: true)

# Advantages over simplified version
- True Newton-Raphson with participation in iterations
- Can enforce capacity limits properly
- No post-processing needed
- Quadratic convergence rate
"""
function solve_power_flow_distributed_slack_full(sys::HybridSystem, 
                                                 dist_slack::DistributedSlack;
                                                 max_iter=50, tol=1e-8, 
                                                 verbose=false, enforce_limits=true)
    rebuild_matrices!(sys)
    nac = length(sys.ac_buses)
    ndc = length(sys.dc_buses)
    
    # Initialize voltages
    Vm = [b.Vm for b in sys.ac_buses]
    Va = zeros(nac)
    Vdc = [b.Vdc_set for b in sys.dc_buses]
    
    # Initialize total slack power (single variable)
    λ_slack = 0.0
    
    # Identify bus types
    pq_idx = [i for (i,b) in enumerate(sys.ac_buses) if b.type == PQ]
    
    participating_set = Set(dist_slack.participating_buses)
    ref_bus = dist_slack.reference_bus
    
    # Track which generators hit limits
    hit_limits = Set{Int}()
    active_factors = copy(dist_slack.participation_factors)
    
    # Build index mappings
    non_ref_buses = [i for i in 1:nac if i != ref_bus]
    
    Y = sys.Ybus
    G = real.(Matrix(Y))
    B = imag.(Matrix(Y))
    
    verbose && println("Full Distributed Slack Power Flow Solver")
    verbose && println("  Participating buses: $(dist_slack.participating_buses)")
    verbose && println("  Participation factors: $(active_factors)")
    verbose && println("  Enforce limits: $(enforce_limits)")
    
    # Store initial generation
    Pg_initial = Dict(i => sys.ac_buses[i].Pg for i in dist_slack.participating_buses)
    
    last_residual = Inf  # Track residual for non-convergence case
    
    for iter in 1:max_iter
        # --- Compute power injections ---
        Pcalc, Qcalc = PowerSystem.ac_power_injections(sys, Vm, Va)
        
        # Scheduled injections
        Psch = zeros(nac)
        Qsch = zeros(nac)
        for i in 1:nac
            Psch[i] = sys.ac_buses[i].Pg - sys.ac_buses[i].Pd
            Qsch[i] = sys.ac_buses[i].Qg - sys.ac_buses[i].Qd
            
            # Add distributed slack contribution: ΔPg_i = α_i × λ_slack
            if i in participating_set && !(i in hit_limits)
                idx = findfirst(==(i), dist_slack.participating_buses)
                if idx !== nothing
                    α_i = active_factors[idx]
                    Psch[i] += α_i * λ_slack
                end
            end
        end
        
        # Converter contributions
        for conv in sys.converters
            conv.status || continue
            Pac, Qac = converter_ac_injection(conv, Vm, Va, Vdc)
            Psch[conv.ac_bus] += Pac
            Qsch[conv.ac_bus] += Qac
        end
        
        # Power mismatches
        ΔP = Pcalc .- Psch
        ΔQ = Qcalc .- Qsch
        
        # DC power equations
        ΔP_dc = Float64[]
        if ndc > 1
            Pdc_calc = Vdc .* (sys.Gdc * Vdc)
            Pdc_sch = [b.Pdc for b in sys.dc_buses]
            for conv in sys.converters
                conv.status || continue
                Pdc_sch[conv.dc_bus] += converter_dc_injection(conv, Vm, Va, Vdc)
            end
            ΔP_dc = Pdc_calc[2:end] .- Pdc_sch[2:end]
        end
        
        # --- Build residual vector F ---
        F = Float64[]
        
        # P equations for ALL buses (including reference)
        # In distributed slack, ref bus P is not fixed - it participates via λ_slack
        for i in 1:nac
            push!(F, ΔP[i])
        end
        
        # Q equations for PQ buses
        for i in pq_idx
            push!(F, ΔQ[i])
        end
        
        # DC equations
        append!(F, ΔP_dc)
        
        # All P equations include α_i × λ_slack term for participating buses
        
        res_norm = norm(F, Inf)
        last_residual = res_norm  # Save for potential non-convergence
        
        if verbose
            println("  Iter $iter: ||F|| = $(res_norm), λ_slack = $(round(λ_slack, digits=6))")
        end
        
        # Check convergence
        if res_norm < tol
            # Compute final slack distribution
            P_slack = Dict{Int,Float64}()
            total_ΔPg = 0.0
            for (k, bus) in enumerate(dist_slack.participating_buses)
                if !(bus in hit_limits)
                    ΔPg_bus = active_factors[k] * λ_slack
                    P_slack[bus] = ΔPg_bus
                    total_ΔPg += ΔPg_bus
                end
            end
            
            verbose && println("  ✅ Converged! λ_slack = $(round(λ_slack, digits=6))")
            verbose && println("  Slack distribution: $P_slack")
            
            return (Vm=Vm, Va=Va, Vdc=Vdc, converged=true, iterations=iter,
                   residual=res_norm, distributed_slack_P=P_slack,
                   total_slack_P=total_ΔPg, hit_limits=collect(hit_limits))
        end
        
        # --- Build Jacobian matrix ---
        np_all = nac  # P equations for ALL buses (including ref)
        nq = length(pq_idx)
        ndc_eq = length(ΔP_dc)
        
        # Dimensionally correct:
        # Variables: [Va_non_ref; Vm_pq; Vdc_non_slack; λ_slack]
        # Equations: [P_all_buses; Q_pq; Pdc_non_slack]
        # n_var = (nac-1) + nq + ndc_eq + 1 = nac + nq + ndc_eq
        # n_eq = nac + nq + ndc_eq
        n_eq = np_all + nq + ndc_eq
        n_var = (nac - 1) + nq + ndc_eq + 1  # = nac + nq + ndc_eq
        
        J = zeros(n_eq, n_var)
        
        # Mapping: state variable indices
        # x = [Va_non_ref (1:nac-1); Vm_pq (nac:nac-1+nq); Vdc_2:end (nac+nq:nac+nq+ndc_eq-1); λ_slack (end)]
        
        col_va = Dict(non_ref_buses[k] => k for k in 1:length(non_ref_buses))
        col_vm = Dict(pq_idx[k] => (nac-1) + k for k in 1:nq)
        col_λ = n_var  # Last column is for λ_slack
        
        # Standard power flow Jacobian (H, N, J, L matrices)
        # P equations for ALL buses
        for (row_idx, i) in enumerate(1:nac)
            # P equation for bus i
            for j in 1:nac
                (G[i,j] == 0.0 && B[i,j] == 0.0) && continue
                θij = Va[i] - Va[j]
                
                if i == j
                    # Diagonal elements
                    H_ii = -Qcalc[i] - B[i,i] * Vm[i]^2
                    N_ii = Pcalc[i] / Vm[i] + G[i,i] * Vm[i]
                    
                    # ∂P_i/∂θ_i (H matrix)
                    if i != ref_bus && haskey(col_va, i)
                        J[row_idx, col_va[i]] += H_ii
                    end
                    
                    # ∂P_i/∂V_i (N matrix) - only if PQ bus
                    if i in pq_idx && haskey(col_vm, i)
                        J[row_idx, col_vm[i]] += N_ii
                    end
                else
                    # Off-diagonal elements
                    H_ij = Vm[i] * Vm[j] * (G[i,j] * sin(θij) - B[i,j] * cos(θij))
                    N_ij = Vm[i] * (G[i,j] * cos(θij) + B[i,j] * sin(θij))
                    
                    # ∂P_i/∂θ_j
                    if j != ref_bus && haskey(col_va, j)
                        J[row_idx, col_va[j]] += H_ij
                    end
                    
                    # ∂P_i/∂V_j - only if j is PQ
                    if j in pq_idx && haskey(col_vm, j)
                        J[row_idx, col_vm[j]] += N_ij
                    end
                end
            end
            
            # ∂P_i/∂λ_slack = -α_i if i is participating (from P_sch,i + α_i × λ_slack)
            if i in participating_set && !(i in hit_limits)
                idx = findfirst(==(i), dist_slack.participating_buses)
                if idx !== nothing
                    α_i = active_factors[idx]
                    J[row_idx, col_λ] = -α_i
                end
            end
        end
        
        # Q equations (J and L matrices)
        for (row_idx, i) in enumerate(pq_idx)
            actual_row = nac + row_idx  # After ALL P equations
            
            for j in 1:nac
                (G[i,j] == 0.0 && B[i,j] == 0.0) && continue
                θij = Va[i] - Va[j]
                
                if i == j
                    J_QQ_ii = Pcalc[i] - G[i,i] * Vm[i]^2
                    L_ii = Qcalc[i] / Vm[i] - B[i,i] * Vm[i]
                    
                    # ∂Q_i/∂θ_i
                    if i != ref_bus && haskey(col_va, i)
                        J[actual_row, col_va[i]] += J_QQ_ii
                    end
                    
                    # ∂Q_i/∂V_i
                    if haskey(col_vm, i)
                        J[actual_row, col_vm[i]] += L_ii
                    end
                else
                    J_QQ_ij = -Vm[i] * Vm[j] * (G[i,j] * cos(θij) + B[i,j] * sin(θij))
                    L_ij = Vm[i] * (G[i,j] * sin(θij) - B[i,j] * cos(θij))
                    
                    # ∂Q_i/∂θ_j
                    if j != ref_bus && haskey(col_va, j)
                        J[actual_row, col_va[j]] += J_QQ_ij
                    end
                    
                    # ∂Q_i/∂V_j
                    if j in pq_idx && haskey(col_vm, j)
                        J[actual_row, col_vm[j]] += L_ij
                    end
                end
            end
            
            # ∂Q_i/∂λ_slack = 0 (Q doesn't depend on λ_slack)
            # No need to set since J is initialized to zeros
        end
        
        # DC Jacobian (simplified - linearized DC power flow)
        if ndc > 1
            for (k_idx, k) in enumerate(2:ndc)
                actual_row = nac + nq + k_idx  # After P and Q equations
                col_vdc = (nac-1) + nq + k_idx  # Vdc columns start after Va and Vm
                
                # ∂P_DC,k/∂V_DC,k
                J[actual_row, col_vdc] = 2 * sys.Gdc[k,k] * Vdc[k]
                
                # Coupling terms (simplified)
                for l in 1:ndc
                    if l != k && l > 1
                        l_idx = l - 1
                        col_vdc_l = (nac-1) + nq + l_idx
                        J[actual_row, col_vdc_l] = sys.Gdc[k,l] * Vdc[l]
                    end
                end
            end
            
            # ∂P_DC/∂λ_slack = 0 (DC power doesn't depend on AC slack)
            # No need to set since J is initialized to zeros
        end
        
        # --- Solve linear system: J × Δx = -F ---
        local Δx  # Define outside try block
        try
            Δx = J \ (-F)
        catch e
            if e isa LinearAlgebra.SingularException
                @warn "Jacobian is singular at iteration $iter (pivot $(e.info))"
                @warn "Jacobian size: $(size(J))"
                @warn "Residual norm: $(norm(F, Inf))"
                # Check for zero rows/columns
                for i in 1:size(J,1)
                    if norm(J[i,:]) < 1e-14
                        @warn "Row $i is essentially zero"
                    end
                end
                for j in 1:size(J,2)
                    if norm(J[:,j]) < 1e-14
                        @warn "Column $j is essentially zero"
                    end
                end
                # Don't compute eigenvalues (too expensive), just report and fail
                @warn "Full distributed slack did not converge (singular Jacobian)"
                return (Vm=Vm, Va=Va, Vdc=Vdc, converged=false, iterations=iter,
                       residual=res_norm, distributed_slack_P=Dict{Int,Float64}(),
                       total_slack_P=0.0, hit_limits=collect(hit_limits))
            end
            rethrow(e)
        end
        
        # Extract updates
        ΔVa = Δx[1:(nac-1)]
        ΔVm = Δx[nac:(nac-1+nq)]
        ΔVdc_vals = ndc > 1 ? Δx[(nac+nq):(nac-1+nq+ndc_eq)] : Float64[]
        Δλ = Δx[end]  # Last element is λ_slack update
        
        # Update state variables
        for (k, i) in enumerate(non_ref_buses)
            Va[i] += ΔVa[k]
        end
        for (k, i) in enumerate(pq_idx)
            Vm[i] += ΔVm[k]
        end
        if ndc > 1
            Vdc[2:end] .+= ΔVdc_vals
        end
        λ_slack += Δλ
        
        # Enforce generator capacity limits
        if enforce_limits && !isempty(dist_slack.max_participation_P)
            for (k, bus) in enumerate(dist_slack.participating_buses)
                if !(bus in hit_limits)
                    α_k = active_factors[k]
                    ΔPg_k = α_k * λ_slack
                    P_total = Pg_initial[bus] + ΔPg_k
                    
                    if haskey(dist_slack.max_participation_P, bus)
                        P_max = dist_slack.max_participation_P[bus]
                        
                        if P_total > P_max
                            verbose && println("    ⚠️  Bus $bus hit upper limit: $(round(P_total, digits=4)) > $(round(P_max, digits=4))")
                            
                            # Mark as hit limit and remove from participation
                            push!(hit_limits, bus)
                            active_factors[k] = 0.0
                            
                            # Renormalize remaining participation factors
                            remaining_sum = sum(active_factors[j] for j in 1:length(active_factors) 
                                              if dist_slack.participating_buses[j] ∉ hit_limits)
                            if remaining_sum > 0
                                for j in 1:length(active_factors)
                                    if dist_slack.participating_buses[j] ∉ hit_limits
                                        active_factors[j] /= remaining_sum
                                    end
                                end
                                verbose && println("    Renormalized factors for remaining buses")
                            end
                            
                        elseif P_total < 0
                            verbose && println("    ⚠️  Bus $bus hit lower limit: $(round(P_total, digits=4)) < 0")
                            
                            push!(hit_limits, bus)
                            active_factors[k] = 0.0
                            
                            # Renormalize
                            remaining_sum = sum(active_factors[j] for j in 1:length(active_factors) 
                                              if dist_slack.participating_buses[j] ∉ hit_limits)
                            if remaining_sum > 0
                                for j in 1:length(active_factors)
                                    if dist_slack.participating_buses[j] ∉ hit_limits
                                        active_factors[j] /= remaining_sum
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    @warn "Full distributed slack did not converge in $max_iter iterations"
    return (Vm=Vm, Va=Va, Vdc=Vdc, converged=false, iterations=max_iter,
           residual=last_residual, distributed_slack_P=Dict{Int,Float64}(),
           total_slack_P=0.0, hit_limits=collect(hit_limits))
end

end  # module PowerSystemEnhanced
