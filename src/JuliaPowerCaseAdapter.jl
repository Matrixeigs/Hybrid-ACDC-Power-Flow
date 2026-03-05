# ═══════════════════════════════════════════════════════════════════════════════
# JuliaPowerCaseAdapter.jl
#
# Provides seamless integration between JuliaPowerCase data structures and
# HybridACDCPowerFlow solver. Enables using JuliaPowerCase.HybridPowerCaseData
# directly with solve_power_flow().
#
# v0.6.0: Updated to use shared JuliaPowerCase types (keyword constructors)
# ═══════════════════════════════════════════════════════════════════════════════

# Import only needed symbols to avoid conflicts with PowerSystem, detect_islands, etc.
using JuliaPowerCase: HybridPowerCaseData, PowerCaseData, AC, DC,
                      Bus, Branch, Generator, VSCConverter,
                      BusType, PQ_BUS, PV_BUS, REF_BUS, ISOLATED_BUS,
                      GenModel, POLYNOMIAL_MODEL,
                      BusSchema, BranchSchema, GenSchema, DCBusSchema, DCBranchSchema, VSCSchema,
                      nrows, nbuses, nbranches

using .PowerSystem: ACBus, ACBranch, DCBus, DCBranch, HybridSystem,
                    PQ, PV, SLACK, PQ_MODE, VDC_Q, VDC_VAC

# Import solve_power_flow to extend it with new method
import .PowerSystem: solve_power_flow

export to_hybrid_system, update_results!

# ═══════════════════════════════════════════════════════════════════════════════
#  TYPE CONVERSION HELPERS
# ═══════════════════════════════════════════════════════════════════════════════

"""Convert JuliaPowerCase bus type (1=PQ, 2=PV, 3=ref) to BusType enum."""
function _convert_bus_type(type_int::Int)::BusType
    type_int == 1 && return PQ_BUS
    type_int == 2 && return PV_BUS
    type_int == 3 && return REF_BUS
    return PQ_BUS  # default
end

"""Convert JuliaPowerCase control mode (1=PQ, 2=VDC_Q, 3=VDC_VAC) to Symbol."""
function _convert_control_mode(mode_int::Int)::Symbol
    mode_int == 1 && return :pq
    mode_int == 2 && return :vdc_q
    mode_int == 3 && return :vdc_vac
    return :pq  # default
end

# ═══════════════════════════════════════════════════════════════════════════════
#  CONVERSION: HybridPowerCaseData → HybridSystem
# ═══════════════════════════════════════════════════════════════════════════════

"""
    to_hybrid_system(h::HybridPowerCaseData) -> HybridSystem

Convert JuliaPowerCase `HybridPowerCaseData` to internal `HybridSystem` for power flow.

# Notes
- Bus IDs are remapped to contiguous 1:N indices
- Generator power is extracted into separate Generator objects
- Per-unit values converted to MW/MVar (using base_mva)
"""
function to_hybrid_system(h::HybridPowerCaseData{T}) where T
    base_mva = Float64(h.base_mva)
    ac = h.ac
    
    # Build AC bus ID mapping using symbolic column access
    nbus = nrows(ac.bus)
    orig_bus_ids = [Int(ac.bus[i, :I]) for i in 1:nbus]
    busid_to_idx = Dict{Int, Int}(id => i for (i, id) in enumerate(orig_bus_ids))
    
    # Aggregate generator power per bus for bus type assignment
    gen_p = zeros(Float64, nbus)
    gen_q = zeros(Float64, nbus)
    gen_vm = ones(Float64, nbus)
    gen_status = zeros(Int, nbus)
    
    ngen = nrows(ac.gen)
    generators = Generator[]
    gen_id = 1
    
    if ngen > 0
        for i in 1:ngen
            gbus_orig = Int(ac.gen[i, :GEN_BUS])
            if haskey(busid_to_idx, gbus_orig)
                gbus_idx = busid_to_idx[gbus_orig]
                status = Int(ac.gen[i, :GEN_STATUS])
                if status == 1
                    Pg_mw = Float64(ac.gen[i, :PG])
                    Qg_mvar = Float64(ac.gen[i, :QG])
                    Vg = Float64(ac.gen[i, :VG])
                    
                    gen_p[gbus_idx] += Pg_mw / base_mva  # For bus type determination
                    gen_q[gbus_idx] += Qg_mvar / base_mva
                    gen_vm[gbus_idx] = Vg
                    gen_status[gbus_idx] = 1
                    
                    # Determine if this is the slack generator
                    bus_type_int = Int(ac.bus[gbus_idx, :TYPE])
                    is_slack = (bus_type_int == 3)
                    
                    # Create Generator object
                    push!(generators, Generator(
                        index          = gen_id,
                        name           = "Gen$gen_id",
                        bus            = gbus_idx,
                        in_service     = true,
                        gen_type       = :thermal,
                        generator_type = "synchronous",
                        pg_mw          = Pg_mw,
                        qg_mvar        = Qg_mvar,
                        vg_pu          = Vg,
                        mbase_mva      = base_mva,
                        cos_phi        = 0.85,
                        pmax_mw        = max(Pg_mw * 1.5, 10.0),
                        pmin_mw        = 0.0,
                        qmax_mvar      = max(abs(Pg_mw), 10.0),
                        qmin_mvar      = -max(abs(Pg_mw), 10.0),
                        is_slack       = is_slack,
                        controllable   = true,
                        ramp_up_mw_min = 9999.0,
                        ramp_down_mw_min = 9999.0,
                        ramp_agc       = 9999.0,
                        ramp_10        = 9999.0,
                        ramp_30        = 9999.0,
                        t_up_min_h     = 0.0,
                        t_down_min_h   = 0.0,
                        cost_model     = POLYNOMIAL_MODEL,
                        cost_startup   = 0.0,
                        cost_shutdown  = 0.0,
                        cost_coeffs    = (0.0, 20.0, 0.0),
                        co2_emission_rate = 0.0,
                        mtbf_hours     = 8760.0,
                        mttr_hours     = 24.0,
                        t_scheduled_h  = 0.0,
                        forced_outage_rate = 0.02,
                        vn_kv          = 110.0,
                        xd_sub_pu      = 0.2,
                        ra_pu          = 0.0,
                        xd_pu          = 0.0,
                        xq_pu          = 0.0,
                        x0_pu          = 0.0,
                        x_r            = 10.0,
                        r1_pu          = 0.0,
                        x1_pu          = 0.2,
                        r2_pu          = 0.0,
                        x2_pu          = 0.2,
                        r0_pu          = 0.0,
                        efficiency     = 1.0
                    ))
                    gen_id += 1
                end
            end
        end
    end
    
    # Build AC buses using JuliaPowerCase Bus{AC} type
    ac_buses = ACBus[]
    for i in 1:nbus
        bus_type = _convert_bus_type(Int(ac.bus[i, :TYPE]))
        Pd_mw = Float64(ac.bus[i, :PD])
        Qd_mvar = Float64(ac.bus[i, :QD])
        Vm = gen_status[i] == 1 ? gen_vm[i] : Float64(ac.bus[i, :VM])
        Va_deg = Float64(ac.bus[i, :VA])
        area = Int(ac.bus[i, :AREA])
        
        push!(ac_buses, Bus{AC}(
            index      = i,
            name       = "Bus$i",
            bus_id     = i,
            in_service = true,
            base_kv    = 110.0,
            bus_type   = bus_type,
            vm_pu      = Vm,
            va_deg     = Va_deg,
            vmax_pu    = 1.1,
            vmin_pu    = 0.9,
            pd_mw      = Pd_mw,
            qd_mvar    = Qd_mvar,
            gs_mw      = 0.0,
            bs_mvar    = 0.0,
            area       = area,
            zone       = 1,
            carbon_area = 1,
            carbon_zone = 1,
            nc         = 0,
            omega      = 0.0,
            is_load    = (Pd_mw > 0)
        ))
    end
    
    # Build AC branches using JuliaPowerCase Branch{AC} type
    nbr = nrows(ac.branch)
    ac_branches = ACBranch[]
    
    for i in 1:nbr
        f_orig = Int(ac.branch[i, :F_BUS])
        t_orig = Int(ac.branch[i, :T_BUS])
        
        haskey(busid_to_idx, f_orig) || continue
        haskey(busid_to_idx, t_orig) || continue
        
        f_idx = busid_to_idx[f_orig]
        t_idx = busid_to_idx[t_orig]
        r = Float64(ac.branch[i, :BR_R])
        x = Float64(ac.branch[i, :BR_X])
        b = Float64(ac.branch[i, :BR_B])
        tap = Float64(ac.branch[i, :TAP])
        tap = tap == 0.0 ? 1.0 : tap
        status = Int(ac.branch[i, :STATUS]) == 1
        
        push!(ac_branches, Branch{AC}(
            index        = i,
            name         = "Line$i",
            from_bus     = f_idx,
            to_bus       = t_idx,
            in_service   = status,
            branch_type  = :line,
            length_km    = 0.0,
            n_parallel   = 1,
            v_rated_kv   = 110.0,
            s_rated_mva  = 100.0,
            s_max_mva    = 100.0,
            r_pu         = r,
            x_pu         = x,
            b_pu         = b,
            r_ohm_km     = 0.0,
            x_ohm_km     = 0.0,
            c_nf_km      = 0.0,
            b_us_km      = 0.0,
            r0_pu        = 0.0,
            x0_pu        = 0.0,
            b0_pu        = 0.0,
            c0_nf_km     = 0.0,
            rate_a_mva   = 100.0,
            rate_b_mva   = 100.0,
            rate_c_mva   = 100.0,
            tap          = tap,
            shift_deg    = 0.0,
            angmin_deg   = -360.0,
            angmax_deg   = 360.0,
            mtbf_hours   = 0.0,
            mttr_hours   = 0.0,
            t_scheduled_h = 0.0,
            sw_hours     = 0.0,
            rp_hours     = 0.0
        ))
    end
    
    # Build DC buses using JuliaPowerCase Bus{DC} type
    dc_buses = DCBus[]
    dc_busid_to_idx = Dict{Int, Int}()
    
    ndcbus = nrows(h.dc_bus)
    if ndcbus > 0
        for i in 1:ndcbus
            dc_id = Int(h.dc_bus[i, :I])
            dc_busid_to_idx[dc_id] = i
            
            Vdc = Float64(h.dc_bus[i, :VDC])
            Pdc_mw = Float64(h.dc_bus[i, :PD])
            
            push!(dc_buses, Bus{DC}(
                index      = i,
                name       = "DCBus$i",
                bus_id     = i,
                in_service = true,
                base_kv    = 200.0,
                bus_type   = PQ_BUS,
                vm_pu      = Vdc,
                va_deg     = 0.0,
                vmax_pu    = 1.1,
                vmin_pu    = 0.9,
                pd_mw      = Pdc_mw,
                qd_mvar    = 0.0,
                gs_mw      = 0.0,
                bs_mvar    = 0.0,
                area       = 1,
                zone       = 1,
                carbon_area = 1,
                carbon_zone = 1,
                nc         = 0,
                omega      = 0.0,
                is_load    = (Pdc_mw != 0)
            ))
        end
    end
    
    # Build DC branches using JuliaPowerCase Branch{DC} type
    dc_branches = Branch{DC}[]
    
    ndcbr = nrows(h.dc_branch)
    if ndcbr > 0
        for i in 1:ndcbr
            f_orig = Int(h.dc_branch[i, :F_BUS])
            t_orig = Int(h.dc_branch[i, :T_BUS])
            
            haskey(dc_busid_to_idx, f_orig) || continue
            haskey(dc_busid_to_idx, t_orig) || continue
            
            f_idx = dc_busid_to_idx[f_orig]
            t_idx = dc_busid_to_idx[t_orig]
            r = Float64(h.dc_branch[i, :BR_R])
            status = Int(h.dc_branch[i, :BR_STATUS]) == 1
            
            push!(dc_branches, Branch{DC}(
                index        = i,
                name         = "DCLine$i",
                from_bus     = f_idx,
                to_bus       = t_idx,
                in_service   = status,
                branch_type  = :dc_line,
                length_km    = 0.0,
                n_parallel   = 1,
                v_rated_kv   = 200.0,
                s_rated_mva  = 100.0,
                s_max_mva    = 100.0,
                r_pu         = r,
                x_pu         = 0.0,
                b_pu         = 0.0,
                r_ohm_km     = 0.0,
                x_ohm_km     = 0.0,
                c_nf_km      = 0.0,
                b_us_km      = 0.0,
                r0_pu        = 0.0,
                x0_pu        = 0.0,
                b0_pu        = 0.0,
                c0_nf_km     = 0.0,
                rate_a_mva   = 100.0,
                rate_b_mva   = 100.0,
                rate_c_mva   = 100.0,
                tap          = 1.0,
                shift_deg    = 0.0,
                angmin_deg   = -360.0,
                angmax_deg   = 360.0,
                mtbf_hours   = 0.0,
                mttr_hours   = 0.0,
                t_scheduled_h = 0.0,
                sw_hours     = 0.0,
                rp_hours     = 0.0
            ))
        end
    end
    
    # Build VSC converters using JuliaPowerCase VSCConverter type
    converters = VSCConverter[]
    
    nvsc = nrows(h.vsc)
    if nvsc > 0
        for i in 1:nvsc
            ac_bus_orig = Int(h.vsc[i, :BUS_AC])
            dc_bus_orig = Int(h.vsc[i, :BUS_DC])
            
            haskey(busid_to_idx, ac_bus_orig) || continue
            haskey(dc_busid_to_idx, dc_bus_orig) || continue
            
            ac_bus_idx = busid_to_idx[ac_bus_orig]
            dc_bus_idx = dc_busid_to_idx[dc_bus_orig]
            
            # Power setpoints (already in MW/MVar in the schema)
            P_mw = Float64(h.vsc[i, :P_MW])
            Q_mvar = Float64(h.vsc[i, :Q_MVAR])
            Vac_set = Float64(h.vsc[i, :VM_AC_PU])
            Vdc_set = Float64(h.vsc[i, :VM_DC_PU])
            
            # Loss model
            loss_pct = Float64(h.vsc[i, :LOSS_PERCENT])
            eta = 1.0 - loss_pct / 100.0
            
            # Control mode (now returns Symbol)
            mode = _convert_control_mode(Int(h.vsc[i, :CONTROL_MODE]))
            
            # Capacity
            Pmax = Float64(h.vsc[i, :PMAX])
            Pmin = Float64(h.vsc[i, :PMIN])
            Smax = max(abs(Pmax), abs(Pmin))
            
            # Droop
            droop = Float64(h.vsc[i, :DROOP_KV])
            
            # Status
            status = Int(h.vsc[i, :IN_SERVICE]) == 1
            
            push!(converters, VSCConverter(
                index        = i,
                name         = "VSC$i",
                bus_ac       = ac_bus_idx,
                bus_dc       = dc_bus_idx,
                in_service   = status,
                vsc_type     = "standard",
                p_rated_mw   = Smax,
                vn_ac_kv     = 110.0,
                vn_dc_kv     = 200.0,
                p_mw         = 0.0,
                q_mvar       = 0.0,
                vm_ac_pu     = Vac_set,
                vm_dc_pu     = Vdc_set,
                pmax_mw      = Pmax,
                pmin_mw      = Pmin,
                qmax_mvar    = Smax,
                qmin_mvar    = -Smax,
                eta          = eta,
                loss_percent = loss_pct,
                loss_mw      = 0.0,
                controllable = true,
                control_mode = mode,
                p_set_mw     = P_mw,
                q_set_mvar   = Q_mvar,
                v_ac_set_pu  = Vac_set,
                v_dc_set_pu  = Vdc_set,
                k_vdc        = droop,
                k_p          = 0.0,
                k_q          = 0.0,
                v_ref_pu     = 1.0,
                f_ref_hz     = 50.0,
                mtbf_hours   = 0.0,
                mttr_hours   = 0.0,
                t_scheduled_h = 0.0
            ))
        end
    end
    
    # Build HybridSystem with generators
    return HybridSystem(ac_buses, ac_branches, dc_buses, dc_branches, converters;
                        generators=generators, baseMVA=base_mva)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  UPDATE RESULTS BACK TO JuliaPowerCase
# ═══════════════════════════════════════════════════════════════════════════════

"""
    update_results!(h::HybridPowerCaseData, result) -> HybridPowerCaseData

Update JuliaPowerCase data with power flow results in-place.

# Arguments
- `h`: HybridPowerCaseData to update
- `result`: Power flow result NamedTuple with `Vm`, `Va`, `Vdc` fields
"""
function update_results!(h::HybridPowerCaseData, result)
    ac = h.ac
    nbus = nrows(ac.bus)
    
    # Update AC bus voltages using symbolic column access
    for i in 1:min(nbus, length(result.Vm))
        ac.bus[i, :VM] = result.Vm[i]
        ac.bus[i, :VA] = rad2deg(result.Va[i])
    end
    
    # Update DC bus voltages
    ndcbus = nrows(h.dc_bus)
    if ndcbus > 0 && hasproperty(result, :Vdc) && length(result.Vdc) > 0
        for i in 1:min(ndcbus, length(result.Vdc))
            h.dc_bus[i, :VDC] = result.Vdc[i]
        end
    end
    
    return h
end

# ═══════════════════════════════════════════════════════════════════════════════
#  MAIN INTERFACE: solve_power_flow(::HybridPowerCaseData)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    solve_power_flow(h::HybridPowerCaseData; max_iter=50, tol=1e-8, update=true)

Solve hybrid AC/DC power flow for JuliaPowerCase data.

# Arguments
- `h`: HybridPowerCaseData from JuliaPowerCase
- `max_iter`: Maximum Newton-Raphson iterations (default: 50)
- `tol`: Convergence tolerance (default: 1e-8)
- `update`: If true, update `h` in-place with results (default: true)

# Returns
NamedTuple with fields:
- `Vm`: AC voltage magnitudes (p.u.)
- `Va`: AC voltage angles (rad)
- `Vdc`: DC voltages (p.u.)
- `converged`: Bool
- `iterations`: Int
- `residual`: Float64

# Example
```julia
using JuliaPowerCase, HybridACDCPowerFlow

h = case_hybrid_5ac3dc()
result = solve_power_flow(h)
result.converged  # true
```
"""
function solve_power_flow(h::HybridPowerCaseData; 
                          max_iter::Int=50, 
                          tol::Float64=1e-8,
                          update::Bool=true)
    # Convert to internal format
    sys = to_hybrid_system(h)
    
    # Solve power flow (calls the HybridSystem method)
    result = solve_power_flow(sys; max_iter=max_iter, tol=tol)
    
    # Update JuliaPowerCase data if requested
    if update && result.converged
        update_results!(h, result)
    end
    
    return result
end
