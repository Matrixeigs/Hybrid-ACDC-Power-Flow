"""
Monte Carlo Stochastic Power Flow Simulation — Island-Aware Edition

Tests the improved hybrid AC/DC power flow solver under:
1. Random load levels (±20% variation)
2. Random line statuses (N-1, N-2 contingencies — may create islands)
3. Random converter control modes (PQ, VDC_Q, VDC_VAC)

Key features tested:
- Island detection with has_generators flag
- Dead island skipping (no generators → Vm = 0)
- Auto slack bus selection (largest online generator)
- Capacity-based distributed slack (no droop/AGC)
- Adaptive solver with PV→PQ conversion

Solvers compared per scenario:
A. solve_power_flow_adaptive  (island-aware, single slack per island)
B. solve_power_flow_distributed_slack  (simplified, post-processing slack)
C. solve_power_flow_distributed_slack_full  (full Jacobian, λ_slack variable)

Outputs:
- CSV files with convergence statistics and island info
- Console summary
"""

using Random
using Statistics
using Printf
using Dates
using LinearAlgebra

include("../src/HybridACDCPowerFlow.jl")
using .HybridACDCPowerFlow
using .HybridACDCPowerFlow: PQ_MODE, VDC_Q, VDC_VAC, ACBus, ACBranch, VSCConverter

# Ensure output directories exist
mkpath("results/monte_carlo")
mkpath("results/monte_carlo/csv")
mkpath("results/monte_carlo/logs")

println("="^80)
println("MONTE CARLO STOCHASTIC POWER FLOW — ISLAND-AWARE EDITION")
println("="^80)
println("Start time: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")

# ─── Monte Carlo parameters ──────────────────────────────────────────────────
const N_SCENARIOS = 100
const LOAD_VARIATION = 0.2        # ±20%
const CONTINGENCY_PROB = 0.1      # 10% probability per line
const MAX_OUTAGES = 3             # Up to N-3 contingencies (can split islands)

Random.seed!(42)

# ─── Results storage ──────────────────────────────────────────────────────────
mutable struct ScenarioResult
    scenario_id::Int
    load_factor::Float64
    n_outages::Int
    outage_lines::Vector{Int}
    converter_modes::Vector{Symbol}

    # Island info
    n_islands::Int
    n_dead_islands::Int        # Islands with no generators (skipped)
    n_auto_slack::Int          # Islands where slack was auto-assigned

    # Method A: Adaptive solver (island-aware)
    adaptive_converged::Bool
    adaptive_iterations::Int
    adaptive_time::Float64
    adaptive_max_voltage::Float64
    adaptive_min_voltage::Float64

    # Method B: Simplified distributed slack
    simplified_converged::Bool
    simplified_iterations::Int
    simplified_residual::Float64
    simplified_time::Float64
    simplified_total_slack::Float64
    simplified_max_voltage::Float64
    simplified_min_voltage::Float64

    # Method C: Full Jacobian distributed slack
    full_converged::Bool
    full_iterations::Int
    full_residual::Float64
    full_time::Float64
    full_total_slack::Float64
    full_max_voltage::Float64
    full_min_voltage::Float64

    # Comparison metrics (B vs C)
    voltage_diff_max::Float64
    voltage_diff_rms::Float64
    slack_diff::Float64

    # Error message (if any)
    error_msg::String

    ScenarioResult(id) = new(id, 0.0, 0, Int[], Symbol[],
                            0, 0, 0,
                            false, 0, 0.0, 0.0, 0.0,
                            false, 0, 0.0, 0.0, 0.0, 0.0, 0.0,
                            false, 0, 0.0, 0.0, 0.0, 0.0, 0.0,
                            0.0, 0.0, 0.0, "")
end

results = ScenarioResult[]

"""Generate random load scenario"""
function generate_random_loads(sys::HybridSystem, load_factor::Float64)
    sys_copy = deepcopy(sys)
    # Create new AC buses with modified loads (ACBus is immutable)
    new_ac_buses = ACBus[]
    for bus in sys_copy.ac_buses
        # Random variation around base load
        variation = 1.0 + load_factor * (2.0 * rand() - 1.0)  # ±load_factor
        new_bus = ACBus(
            bus.id,
            bus.type,
            bus.Pd * variation,
            bus.Qd * variation,
            bus.Pg,
            bus.Qg,
            bus.Vm,
            bus.Va,
            bus.area
        )
        push!(new_ac_buses, new_bus)
    end
    sys_copy.ac_buses = new_ac_buses
    return sys_copy
end

"""Generate random line outages"""
function generate_random_outages(sys::HybridSystem, max_outages::Int, prob::Float64)
    outages = Int[]
    n_branches = length(sys.ac_branches)
    
    for i in 1:n_branches
        if rand() < prob && length(outages) < max_outages
            push!(outages, i)
        end
    end
    
    if !isempty(outages)
        sys_copy = deepcopy(sys)
        # Create new branches with modified impedance for outages
        new_ac_branches = ACBranch[]
        for (i, branch) in enumerate(sys_copy.ac_branches)
            if i in outages
                # Create branch outage by setting status to false
                new_branch = ACBranch(
                    branch.from,
                    branch.to,
                    branch.r,
                    branch.x,
                    branch.b,
                    branch.tap,
                    false  # status = false for outage
                )
            else
                new_branch = branch
            end
            push!(new_ac_branches, new_branch)
        end
        sys_copy.ac_branches = new_ac_branches
        # Rebuild admittance matrix
        sys_copy.Ybus = build_admittance_matrix(sys_copy)
        return sys_copy, outages
    end
    
    return sys, outages
end

"""Generate random converter modes"""
function generate_random_converter_modes(sys::HybridSystem)
    sys_copy = deepcopy(sys)
    modes = Symbol[]
    
    mode_options = [PQ_MODE, VDC_Q, VDC_VAC]
    
    # Create new converters with random modes
    new_converters = VSCConverter[]
    for conv in sys_copy.converters
        mode = rand(mode_options)
        new_conv = VSCConverter(
            conv.id,
            conv.ac_bus,
            conv.dc_bus,
            mode,
            conv.Pset,
            conv.Qset,
            conv.Vdc_set,
            conv.Vac_set,
            conv.Ploss_a,
            conv.Ploss_b,
            conv.Ploss_c,
            conv.Smax,
            conv.status; G_droop=conv.G_droop
        )
        push!(new_converters, new_conv)
        push!(modes, Symbol(mode))
    end
    sys_copy.converters = new_converters
    
    return sys_copy, modes
end

"""Run power flow with timing"""
function timed_power_flow(solve_func, args...; kwargs...)
    start_time = time()
    result = solve_func(args...; kwargs...)
    elapsed = time() - start_time
    return result, elapsed
end

# ─── Main Monte Carlo loop ────────────────────────────────────────────────────
println("\n🎲 Starting Monte Carlo simulation (island-aware)...")
println("   Scenarios: $N_SCENARIOS")
println("   Load variation: ±$(Int(LOAD_VARIATION*100))%")
println("   Contingency probability: $(Int(CONTINGENCY_PROB*100))%")
println("   Max outages: N-$MAX_OUTAGES")
println()

for scenario_id in 1:N_SCENARIOS
    result = ScenarioResult(scenario_id)

    try
        # Build base system
        sys_base = build_ieee24_3area_acdc()

        # Generate random scenario
        load_factor = LOAD_VARIATION * (2.0 * rand() - 1.0)
        result.load_factor = load_factor

        sys = generate_random_loads(sys_base, load_factor)
        sys, outages = generate_random_outages(sys, MAX_OUTAGES, CONTINGENCY_PROB)
        sys, modes = generate_random_converter_modes(sys)

        result.n_outages = length(outages)
        result.outage_lines = outages
        result.converter_modes = modes

        # ── Island detection ──────────────────────────────────────────────
        islands = detect_islands(sys)
        result.n_islands = length(islands)
        result.n_dead_islands = count(isl -> !isl.has_generators, islands)
        result.n_auto_slack  = count(isl -> isl.has_generators && !isl.has_ac_slack, islands)

        # ── Method A: Adaptive solver (island-aware, auto slack) ──────────
        try
            sys_a = deepcopy(sys)
            opts = PowerFlowOptions(verbose=false, enable_auto_swing_selection=true,
                                    enable_pv_pq_conversion=true,
                                    enable_converter_mode_switching=true)
            t0 = time()
            res_a = solve_power_flow_adaptive(sys_a; options=opts)
            result.adaptive_time = time() - t0
            result.adaptive_converged = res_a.converged
            result.adaptive_iterations = res_a.iterations
            if res_a.converged
                # Only consider buses not in dead islands
                live_buses = Int[]
                for isl in res_a.islands
                    isl.has_generators && append!(live_buses, isl.ac_buses)
                end
                if !isempty(live_buses)
                    result.adaptive_max_voltage = maximum(res_a.Vm[live_buses])
                    result.adaptive_min_voltage = minimum(res_a.Vm[live_buses])
                end
            end
        catch e
            result.adaptive_converged = false
            result.error_msg *= "A:$(sprint(showerror, e)) "
        end

        # ── Method B: Simplified distributed slack (island-aware) ─────────
        try
            sys_b = deepcopy(sys)
            t0_b = time()

            nac_b = length(sys_b.ac_buses)
            ndc_b = length(sys_b.dc_buses)
            Vm_b = [b.Vm for b in sys_b.ac_buses]
            Va_b = zeros(nac_b)
            Vdc_b = ndc_b > 0 ? [b.Vdc_set for b in sys_b.dc_buses] : Float64[]
            total_slack_b = 0.0
            conv_b = true
            iters_b = 0
            res_b_residual = 0.0

            for isl in islands
                if !isl.has_generators
                    for i in isl.ac_buses; Vm_b[i] = 0.0; Va_b[i] = 0.0; end
                    for i in isl.dc_buses; i <= ndc_b && (Vdc_b[i] = 0.0); end
                    continue
                end
                # Power balance feasibility check
                sum_Pg_b = sum(sys_b.ac_buses[i].Pg for i in isl.ac_buses)
                sum_Pd_b = sum(sys_b.ac_buses[i].Pd for i in isl.ac_buses)
                if sum_Pd_b > sum_Pg_b + 1e-6
                    for i in isl.ac_buses; Vm_b[i] = 0.0; Va_b[i] = 0.0; end
                    for i in isl.dc_buses; i <= ndc_b && (Vdc_b[i] = 0.0); end
                    conv_b = false
                    continue
                end
                # Extract sub-system
                sub_b, ac_m_b, dc_m_b = extract_island_subsystem(sys_b, isl)
                # Auto-assign slack on sub-system
                if !any(b.type == SLACK for b in sub_b.ac_buses)
                    gen_buses = [(i, sub_b.ac_buses[i].Pg) for i in 1:length(sub_b.ac_buses)
                                 if sub_b.ac_buses[i].Pg > 0]
                    if !isempty(gen_buses)
                        sort!(gen_buses, by=x -> x[2], rev=true)
                        best = gen_buses[1][1]
                        old = sub_b.ac_buses[best]
                        sub_b.ac_buses[best] = ACBus(old.id, SLACK, old.Pd, old.Qd,
                                                      old.Pg, old.Qg, old.Vm, 0.0, old.area)
                    end
                end
                # Create participation factors for sub-system
                local dist_b
                try
                    dist_b = create_participation_factors(sub_b; method=:capacity)
                catch
                    # Fallback: solve with adaptive if no participating generators
                    res_sub = PowerSystem.solve_power_flow(sub_b; max_iter=100, tol=1e-8)
                    for (orig, loc) in ac_m_b; Vm_b[orig] = res_sub.Vm[loc]; Va_b[orig] = res_sub.Va[loc]; end
                    for (orig, loc) in dc_m_b; orig <= ndc_b && loc <= length(res_sub.Vdc) && (Vdc_b[orig] = res_sub.Vdc[loc]); end
                    conv_b &= res_sub.converged; iters_b += res_sub.iterations
                    continue
                end
                res_sub = solve_power_flow_distributed_slack(sub_b, dist_b;
                                                             verbose=false, max_iter=100, tol=1e-8)
                for (orig, loc) in ac_m_b; Vm_b[orig] = res_sub.Vm[loc]; Va_b[orig] = res_sub.Va[loc]; end
                for (orig, loc) in dc_m_b; orig <= ndc_b && loc <= length(res_sub.Vdc) && (Vdc_b[orig] = res_sub.Vdc[loc]); end
                conv_b &= res_sub.converged
                iters_b += res_sub.iterations
                res_b_residual = max(res_b_residual, res_sub.residual)
                total_slack_b += sum(values(res_sub.distributed_slack_P); init=0.0)
            end

            result.simplified_time = time() - t0_b
            result.simplified_converged = conv_b
            result.simplified_iterations = iters_b
            result.simplified_residual = res_b_residual
            result.simplified_total_slack = total_slack_b
            live_b = [i for isl in islands if isl.has_generators for i in isl.ac_buses]
            if !isempty(live_b)
                result.simplified_max_voltage = maximum(Vm_b[live_b])
                result.simplified_min_voltage = minimum(Vm_b[live_b])
            end
        catch e
            result.simplified_converged = false
            result.error_msg *= "B:$(sprint(showerror, e)) "
        end

        # ── Method C: Full Jacobian distributed slack (island-aware) ──────
        try
            sys_c = deepcopy(sys)
            t0_c = time()

            nac_c = length(sys_c.ac_buses)
            ndc_c = length(sys_c.dc_buses)
            Vm_c = [b.Vm for b in sys_c.ac_buses]
            Va_c = zeros(nac_c)
            Vdc_c = ndc_c > 0 ? [b.Vdc_set for b in sys_c.dc_buses] : Float64[]
            total_slack_c = 0.0
            conv_c = true
            iters_c = 0
            res_c_residual = 0.0

            for isl in islands
                if !isl.has_generators
                    for i in isl.ac_buses; Vm_c[i] = 0.0; Va_c[i] = 0.0; end
                    for i in isl.dc_buses; i <= ndc_c && (Vdc_c[i] = 0.0); end
                    continue
                end
                # Power balance feasibility check
                sum_Pg_c = sum(sys_c.ac_buses[i].Pg for i in isl.ac_buses)
                sum_Pd_c = sum(sys_c.ac_buses[i].Pd for i in isl.ac_buses)
                if sum_Pd_c > sum_Pg_c + 1e-6
                    for i in isl.ac_buses; Vm_c[i] = 0.0; Va_c[i] = 0.0; end
                    for i in isl.dc_buses; i <= ndc_c && (Vdc_c[i] = 0.0); end
                    conv_c = false
                    continue
                end
                # Extract sub-system
                sub_c, ac_m_c, dc_m_c = extract_island_subsystem(sys_c, isl)
                # Auto-assign slack on sub-system
                if !any(b.type == SLACK for b in sub_c.ac_buses)
                    gen_buses = [(i, sub_c.ac_buses[i].Pg) for i in 1:length(sub_c.ac_buses)
                                 if sub_c.ac_buses[i].Pg > 0]
                    if !isempty(gen_buses)
                        sort!(gen_buses, by=x -> x[2], rev=true)
                        best = gen_buses[1][1]
                        old = sub_c.ac_buses[best]
                        sub_c.ac_buses[best] = ACBus(old.id, SLACK, old.Pd, old.Qd,
                                                      old.Pg, old.Qg, old.Vm, 0.0, old.area)
                    end
                end
                # Create participation factors for sub-system
                local dist_c
                try
                    dist_c = create_participation_factors(sub_c; method=:capacity)
                catch
                    res_sub = PowerSystem.solve_power_flow(sub_c; max_iter=100, tol=1e-8)
                    for (orig, loc) in ac_m_c; Vm_c[orig] = res_sub.Vm[loc]; Va_c[orig] = res_sub.Va[loc]; end
                    for (orig, loc) in dc_m_c; orig <= ndc_c && loc <= length(res_sub.Vdc) && (Vdc_c[orig] = res_sub.Vdc[loc]); end
                    conv_c &= res_sub.converged; iters_c += res_sub.iterations
                    continue
                end
                res_sub = solve_power_flow_distributed_slack_full(sub_c, dist_c;
                                                                   verbose=false, max_iter=100, tol=1e-8,
                                                                   enforce_limits=true)
                for (orig, loc) in ac_m_c; Vm_c[orig] = res_sub.Vm[loc]; Va_c[orig] = res_sub.Va[loc]; end
                for (orig, loc) in dc_m_c; orig <= ndc_c && loc <= length(res_sub.Vdc) && (Vdc_c[orig] = res_sub.Vdc[loc]); end
                conv_c &= res_sub.converged
                iters_c += res_sub.iterations
                res_c_residual = max(res_c_residual, res_sub.residual)
                total_slack_c += sum(values(res_sub.distributed_slack_P); init=0.0)
            end

            result.full_time = time() - t0_c
            result.full_converged = conv_c
            result.full_iterations = iters_c
            result.full_residual = res_c_residual
            result.full_total_slack = total_slack_c
            live_c = [i for isl in islands if isl.has_generators for i in isl.ac_buses]
            if !isempty(live_c)
                result.full_max_voltage = maximum(Vm_c[live_c])
                result.full_min_voltage = minimum(Vm_c[live_c])
            end

            # Comparison B vs C (if both converged)
            if result.simplified_converged && result.full_converged
                result.slack_diff = abs(result.simplified_total_slack - result.full_total_slack)
            end
        catch e
            result.full_converged = false
            result.error_msg *= "C:$(sprint(showerror, e)) "
        end

        push!(results, result)

        # Progress indicator every 10 scenarios
        if scenario_id % 10 == 0
            n_a = count(r -> r.adaptive_converged, results)
            n_b = count(r -> r.simplified_converged, results)
            n_c = count(r -> r.full_converged, results)
            n_dead = sum(r -> r.n_dead_islands, results)
            @printf("  %3d/%d | A:%d B:%d C:%d | dead islands: %d\n",
                   scenario_id, N_SCENARIOS, n_a, n_b, n_c, n_dead)
        end

    catch e
        result.error_msg = "setup: $(sprint(showerror, e))"
        push!(results, result)
    end
end

println("\n✅ Simulation complete!")
println("\n" * "="^80)
println("RESULTS SUMMARY — ISLAND-AWARE MONTE CARLO")
println("="^80)

n_total = length(results)
n_a_conv = count(r -> r.adaptive_converged, results)
n_b_conv = count(r -> r.simplified_converged, results)
n_c_conv = count(r -> r.full_converged, results)
n_all_conv = count(r -> r.adaptive_converged && r.simplified_converged && r.full_converged, results)

# ── Island statistics ─────────────────────────────────────────────────────────
total_dead = sum(r -> r.n_dead_islands, results)
total_auto = sum(r -> r.n_auto_slack, results)
scenarios_with_islands = count(r -> r.n_islands > 1, results)
scenarios_with_dead = count(r -> r.n_dead_islands > 0, results)

println("\n🏝️  Island Statistics:")
println("  Scenarios with multiple islands: $scenarios_with_islands / $n_total")
println("  Scenarios with dead islands:     $scenarios_with_dead / $n_total")
println("  Total dead islands (all scenarios): $total_dead")
println("  Total auto-slack assignments:       $total_auto")

# ── Convergence statistics ────────────────────────────────────────────────────
println("\n📊 Convergence Statistics:")
println("  Total scenarios:        $n_total")
@printf("  A (Adaptive):           %d / %d  (%.1f%%)\n", n_a_conv, n_total, 100*n_a_conv/n_total)
@printf("  B (Simplified dist):    %d / %d  (%.1f%%)\n", n_b_conv, n_total, 100*n_b_conv/n_total)
@printf("  C (Full Jacobian dist): %d / %d  (%.1f%%)\n", n_c_conv, n_total, 100*n_c_conv/n_total)
@printf("  All three converged:    %d / %d  (%.1f%%)\n", n_all_conv, n_total, 100*n_all_conv/n_total)

# ── Convergence vs outage count ───────────────────────────────────────────────
println("\n📊 Convergence by Outage Level:")
for k in 0:MAX_OUTAGES
    subset = filter(r -> r.n_outages == k, results)
    isempty(subset) && continue
    na = count(r -> r.adaptive_converged, subset)
    nb = count(r -> r.simplified_converged, subset)
    nc = count(r -> r.full_converged, subset)
    @printf("  N-%d (%3d scenarios): A=%d  B=%d  C=%d\n", k, length(subset), na, nb, nc)
end

# ── Performance (converged cases) ─────────────────────────────────────────────
a_conv = filter(r -> r.adaptive_converged, results)
b_conv = filter(r -> r.simplified_converged, results)
c_conv = filter(r -> r.full_converged, results)

println("\n⏱️  Avg Solve Time (converged only):")
if !isempty(a_conv)
    @printf("  A (Adaptive):        %.2f ms  (avg %d iters)\n",
           mean(r -> r.adaptive_time, a_conv) * 1000,
           round(Int, mean(r -> r.adaptive_iterations, a_conv)))
end
if !isempty(b_conv)
    @printf("  B (Simplified):      %.2f ms  (avg %d iters)\n",
           mean(r -> r.simplified_time, b_conv) * 1000,
           round(Int, mean(r -> r.simplified_iterations, b_conv)))
end
if !isempty(c_conv)
    @printf("  C (Full Jacobian):   %.2f ms  (avg %d iters)\n",
           mean(r -> r.full_time, c_conv) * 1000,
           round(Int, mean(r -> r.full_iterations, c_conv)))
end

# ── Voltage statistics ────────────────────────────────────────────────────────
println("\n🔌 Voltage Range (converged, live buses only):")
if !isempty(a_conv)
    @printf("  A: Vmin=%.4f  Vmax=%.4f\n",
           minimum(r -> r.adaptive_min_voltage, a_conv),
           maximum(r -> r.adaptive_max_voltage, a_conv))
end
if !isempty(b_conv)
    @printf("  B: Vmin=%.4f  Vmax=%.4f\n",
           minimum(r -> r.simplified_min_voltage, b_conv),
           maximum(r -> r.simplified_max_voltage, b_conv))
end
if !isempty(c_conv)
    @printf("  C: Vmin=%.4f  Vmax=%.4f\n",
           minimum(r -> r.full_min_voltage, c_conv),
           maximum(r -> r.full_max_voltage, c_conv))
end

# ── Slack comparison B vs C ───────────────────────────────────────────────────
bc_conv = filter(r -> r.simplified_converged && r.full_converged, results)
if !isempty(bc_conv)
    println("\n🔍 B vs C Comparison ($(length(bc_conv)) both-converged scenarios):")
    @printf("  Avg slack diff: %.6f p.u.\n", mean(r -> r.slack_diff, bc_conv))
end

# ── Sample error messages ─────────────────────────────────────────────────────
errors = filter(r -> !isempty(r.error_msg), results)
if !isempty(errors)
    println("\n⚠️  Scenarios with errors: $(length(errors))")
    for r in errors[1:min(5, length(errors))]
        @printf("  Scenario %3d (N-%d, %.0f%% load): %s\n",
               r.scenario_id, r.n_outages, r.load_factor*100, r.error_msg)
    end
    if length(errors) > 5
        println("  ... and $(length(errors)-5) more")
    end
end

# ─── Save CSV ─────────────────────────────────────────────────────────────────
csv_file = "results/monte_carlo/csv/monte_carlo_results.csv"
println("\n💾 Saving results to: $csv_file")

open(csv_file, "w") do io
    println(io, "scenario_id,load_factor,n_outages,n_islands,n_dead_islands,n_auto_slack," *
               "adaptive_conv,adaptive_iter,adaptive_time_ms,adaptive_vmax,adaptive_vmin," *
               "simp_conv,simp_iter,simp_residual,simp_time_ms,simp_total_slack,simp_vmax,simp_vmin," *
               "full_conv,full_iter,full_residual,full_time_ms,full_total_slack,full_vmax,full_vmin," *
               "slack_diff,error")

    for r in results
        @printf(io, "%d,%.4f,%d,%d,%d,%d,%d,%d,%.3f,%.6f,%.6f,%d,%d,%.6e,%.3f,%.6f,%.6f,%.6f,%d,%d,%.6e,%.3f,%.6f,%.6f,%.6f,%.6f,%s\n",
               r.scenario_id, r.load_factor, r.n_outages, r.n_islands, r.n_dead_islands, r.n_auto_slack,
               r.adaptive_converged, r.adaptive_iterations, r.adaptive_time*1000,
               r.adaptive_max_voltage, r.adaptive_min_voltage,
               r.simplified_converged, r.simplified_iterations, r.simplified_residual, r.simplified_time*1000,
               r.simplified_total_slack, r.simplified_max_voltage, r.simplified_min_voltage,
               r.full_converged, r.full_iterations, r.full_residual, r.full_time*1000,
               r.full_total_slack, r.full_max_voltage, r.full_min_voltage,
               r.slack_diff, replace(r.error_msg, "," => ";"))
    end
end

# ─── Save summary ─────────────────────────────────────────────────────────────
summary_file = "results/monte_carlo/csv/summary_island_aware.csv"
println("💾 Saving summary to: $summary_file")

open(summary_file, "w") do io
    println(io, "Metric,Value")
    println(io, "Total Scenarios,$n_total")
    println(io, "Max Outages,N-$MAX_OUTAGES")
    println(io, "Scenarios with Multi-Islands,$scenarios_with_islands")
    println(io, "Scenarios with Dead Islands,$scenarios_with_dead")
    println(io, "Total Dead Islands,$total_dead")
    println(io, "Total Auto-Slack,$total_auto")
    @printf(io, "Adaptive Convergence Rate,%.2f%%\n", 100*n_a_conv/n_total)
    @printf(io, "Simplified Convergence Rate,%.2f%%\n", 100*n_b_conv/n_total)
    @printf(io, "Full Jacobian Convergence Rate,%.2f%%\n", 100*n_c_conv/n_total)
    @printf(io, "All Three Converged,%.2f%%\n", 100*n_all_conv/n_total)
    if !isempty(a_conv)
        @printf(io, "Avg Adaptive Time (ms),%.3f\n", mean(r -> r.adaptive_time, a_conv)*1000)
    end
    if !isempty(b_conv)
        @printf(io, "Avg Simplified Time (ms),%.3f\n", mean(r -> r.simplified_time, b_conv)*1000)
    end
    if !isempty(c_conv)
        @printf(io, "Avg Full Jacobian Time (ms),%.3f\n", mean(r -> r.full_time, c_conv)*1000)
    end
    if !isempty(bc_conv)
        @printf(io, "Avg Slack Diff (B vs C),%.6f\n", mean(r -> r.slack_diff, bc_conv))
    end
end

println("\n✅ Monte Carlo island-aware simulation complete!")
println("   End time: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
println("="^80)
