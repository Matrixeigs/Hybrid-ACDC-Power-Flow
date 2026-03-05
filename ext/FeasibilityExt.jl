"""
    FeasibilityExt

Package extension for HybridACDCPowerFlow that provides feasibility checking
and validation against reference solvers (NLsolve, JuMP+Ipopt).

v0.6.0: Updated for JuliaPowerCase types.
"""
module FeasibilityExt

using HybridACDCPowerFlow
using HybridACDCPowerFlow: HybridSystem, BusType, PQ, PV, SLACK, 
                           PQ_MODE, VDC_Q, VDC_VAC,
                           power_flow_residual, FeasibilityResult,
                           solve_power_flow

using NLsolve
using JuMP
using Ipopt
using Printf
using LinearAlgebra
using SparseArrays

using JuliaPowerCase: PQ_BUS, PV_BUS, REF_BUS

# ═══════════════════════════════════════════════════════════════════════════════
#  NLSOLVE-BASED FEASIBILITY CHECK
# ═══════════════════════════════════════════════════════════════════════════════

"""
    check_power_flow_feasibility_nlsolve(sys::HybridSystem; verbose, max_iter, tol)

Check power flow feasibility using NLsolve root-finding.
Uses trust-region method to solve the same residual equations as our Newton-Raphson.
"""
function HybridACDCPowerFlow.check_power_flow_feasibility_nlsolve(sys::HybridSystem; 
                                                                  verbose::Bool=false, 
                                                                  max_iter::Int=50, 
                                                                  tol::Float64=1e-6)
    nac = length(sys.ac_buses)
    ndc = length(sys.dc_buses)
    baseMVA = sys.baseMVA
    
    # Compute load/generation statistics (in p.u.)
    total_load = sum(bus.pd_mw for bus in sys.ac_buses) / baseMVA
    total_gen = sum(sys.Pg)
    load_margin = total_gen - total_load
    
    verbose && println("\n" * "="^70)
    verbose && println("  NLSOLVE FEASIBILITY CHECK (Root-Finding)")
    verbose && println("="^70)
    verbose && @printf("    Total load: %.4f, Gen: %.4f, Margin: %.4f p.u.\n", 
                      total_load, total_gen, load_margin)
    
    # Quick pre-check
    if load_margin < -0.1
        verbose && println("  ❌ INFEASIBLE: Insufficient generation")
        return FeasibilityResult(
            false, "INSUFFICIENT_GENERATION", 0.0, 0.0, 0,
            zeros(nac), zeros(nac), zeros(ndc), zeros(nac), zeros(nac),
            total_load, total_gen, load_margin
        )
    end
    
    # Identify bus types (using JuliaPowerCase constants)
    pq_idx = findall(b -> b.bus_type == PQ_BUS, sys.ac_buses)
    pv_idx = findall(b -> b.bus_type == PV_BUS, sys.ac_buses)
    slack_idx = findall(b -> b.bus_type == REF_BUS, sys.ac_buses)
    
    # Initialize state: Va (all non-slack), Vm (PQ only), Vdc (non-slack DC)
    Vm_init = [bus.vm_pu for bus in sys.ac_buses]
    Va_init = [deg2rad(bus.va_deg) for bus in sys.ac_buses]
    Vdc_init = ndc > 0 ? [bus.vm_pu for bus in sys.dc_buses] : Float64[]
    
    # Build initial guess vector: [Va_nonslack; Vm_pq; Vdc_non_slack]
    non_slack_ac = sort(union(pq_idx, pv_idx))
    n_Va = length(non_slack_ac)
    n_Vm = length(pq_idx)
    n_Vdc = ndc > 1 ? ndc - 1 : 0
    
    x0 = zeros(n_Va + n_Vm + n_Vdc)
    x0[1:n_Va] .= Va_init[non_slack_ac]
    x0[n_Va+1:n_Va+n_Vm] .= Vm_init[pq_idx]
    if n_Vdc > 0
        x0[n_Va+n_Vm+1:end] .= Vdc_init[2:end]
    end
    
    # Add small random perturbation
    x0 .+= 0.001 * randn(length(x0))
    
    # Define residual function wrapper
    function residual_wrapper!(F, x)
        Vm = copy(Vm_init)
        Va = copy(Va_init)
        Vdc = ndc > 0 ? copy(Vdc_init) : Float64[]
        
        Va[non_slack_ac] .= x[1:n_Va]
        Vm[pq_idx] .= x[n_Va+1:n_Va+n_Vm]
        if n_Vdc > 0
            Vdc[2:end] .= x[n_Va+n_Vm+1:end]
        end
        
        try
            F_ac, F_dc = power_flow_residual(sys, Vm, Va, Vdc)
            F_full = ndc > 0 ? vcat(F_ac, F_dc) : F_ac
            F .= F_full
        catch e
            verbose && println("  Error in residual_wrapper: $e")
            F .= 1e6 .* ones(length(x))
        end
    end
    
    start_time = time()
    
    try
        result = nlsolve(residual_wrapper!, x0; 
                        method=:trust_region,
                        ftol=tol,
                        iterations=max_iter,
                        show_trace=verbose,
                        autodiff=:finite)
        
        solve_time = (time() - start_time) * 1000
        
        if result.f_converged || result.x_converged
            # Extract solution
            Vm = copy(Vm_init)
            Va = copy(Va_init)
            Vdc = ndc > 0 ? copy(Vdc_init) : Float64[]
            
            x_sol = result.zero
            Va[non_slack_ac] .= x_sol[1:n_Va]
            Vm[pq_idx] .= x_sol[n_Va+1:n_Va+n_Vm]
            if n_Vdc > 0
                Vdc[2:end] .= x_sol[n_Va+n_Vm+1:end]
            end
            
            F_ac, F_dc = power_flow_residual(sys, Vm, Va, Vdc)
            F_full = ndc > 0 ? vcat(F_ac, F_dc) : F_ac
            residual_norm = norm(F_full, Inf)
            
            verbose && @printf("  ✅ FEASIBLE: Converged in %d iterations (residual: %.3e)\n", 
                              result.iterations, residual_norm)
            
            return FeasibilityResult(
                true, "CONVERGED", residual_norm, solve_time, result.iterations,
                Vm, Va, Vdc, copy(sys.Pg), copy(sys.Qg),
                total_load, total_gen, load_margin
            )
        else
            F_ac, F_dc = power_flow_residual(sys, Vm_init, Va_init, Vdc_init)
            F_full = ndc > 0 ? vcat(F_ac, F_dc) : F_ac
            residual_norm = norm(F_full, Inf)
            verbose && @printf("  ❌ INFEASIBLE: No convergence (residual: %.3e)\n", residual_norm)
            return FeasibilityResult(
                false, "NO_CONVERGENCE", residual_norm, solve_time, result.iterations,
                zeros(nac), zeros(nac), zeros(ndc), zeros(nac), zeros(nac),
                total_load, total_gen, load_margin
            )
        end
        
    catch e
        solve_time = (time() - start_time) * 1000
        verbose && println("  ❌ INFEASIBLE: Solver error - $e")
        return FeasibilityResult(
            false, "SOLVER_ERROR", 0.0, solve_time, 0,
            zeros(nac), zeros(nac), zeros(ndc), zeros(nac), zeros(nac),
            total_load, total_gen, load_margin
        )
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
#  JUMP/IPOPT-BASED FEASIBILITY CHECK (Full AC Power Flow Formulation)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    check_power_flow_feasibility_jump(sys::HybridSystem; verbose, max_time)

Check power flow feasibility using JuMP+Ipopt nonlinear optimization.
Solves a feasibility problem for the AC power flow equations.
"""
function HybridACDCPowerFlow.check_power_flow_feasibility_jump(sys::HybridSystem; 
                                                               verbose::Bool=false, 
                                                               max_time::Float64=30.0)
    nac = length(sys.ac_buses)
    ndc = length(sys.dc_buses)
    nc = length(sys.converters)
    baseMVA = sys.baseMVA
    
    total_load = sum(bus.pd_mw for bus in sys.ac_buses) / baseMVA
    total_gen = sum(sys.Pg)
    load_margin = total_gen - total_load
    
    verbose && println("\n" * "="^70)
    verbose && println("  JUMP/IPOPT FEASIBILITY CHECK")
    verbose && println("="^70)
    verbose && @printf("    AC buses: %d, DC buses: %d, Converters: %d\n", nac, ndc, nc)
    verbose && @printf("    Total load: %.4f p.u., Gen: %.4f p.u.\n", total_load, total_gen)
    
    if load_margin < -0.1
        verbose && println("  ❌ INFEASIBLE: Insufficient generation capacity")
        return FeasibilityResult(
            false, "INSUFFICIENT_GENERATION", 0.0, 0.0, 0,
            Float64[], Float64[], Float64[], Float64[], Float64[],
            total_load, total_gen, load_margin
        )
    end
    
    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "print_level", verbose ? 5 : 0)
    set_optimizer_attribute(model, "max_iter", 500)
    set_optimizer_attribute(model, "max_cpu_time", max_time)
    set_optimizer_attribute(model, "tol", 1e-6)
    
    # Variables
    @variable(model, 0.9 <= Vm[i=1:nac] <= 1.1, start=sys.ac_buses[i].vm_pu)
    @variable(model, -π <= Va[i=1:nac] <= π, start=deg2rad(sys.ac_buses[i].va_deg))
    
    if ndc > 0
        @variable(model, 0.9 <= Vdc[d=1:ndc] <= 1.1, start=sys.dc_buses[d].vm_pu)
        # Fix DC slack bus voltage
        fix(Vdc[1], sys.dc_buses[1].vm_pu; force=true)
    end
    
    # Find PV/Slack buses and fix their voltages
    ref_bus = nothing
    for (i, bus) in enumerate(sys.ac_buses)
        if bus.bus_type == REF_BUS
            ref_bus = i
            fix(Va[i], deg2rad(bus.va_deg); force=true)
            fix(Vm[i], bus.vm_pu; force=true)
        elseif bus.bus_type == PV_BUS
            fix(Vm[i], bus.vm_pu; force=true)
        end
    end
    
    if ref_bus === nothing
        fix(Va[1], 0.0; force=true)
    end
    
    # Get admittance matrix
    G = real.(Matrix(sys.Ybus))
    B = imag.(Matrix(sys.Ybus))
    
    # Power balance equations (P and Q for each bus)
    for i in 1:nac
        bus = sys.ac_buses[i]
        Pd = bus.pd_mw / baseMVA
        Qd = bus.qd_mvar / baseMVA
        
        # Get generator injection at this bus
        Pg_bus = sys.Pg[i]
        Qg_bus = sys.Qg[i]
        
        # Converter injection (if any)
        conv_P = 0.0
        conv_Q = 0.0
        for (c, conv) in enumerate(sys.converters)
            if conv.bus_ac == i && conv.in_service
                if conv.control_mode == PQ_MODE
                    conv_P += conv.p_set_mw / baseMVA
                    conv_Q += conv.q_set_mvar / baseMVA
                end
            end
        end
        
        # P balance: Pg - Pd + Pconv = sum(Vm_i * Vm_j * (G_ij*cos + B_ij*sin))
        if bus.bus_type != REF_BUS  # Slack bus has unknown P
            @NLconstraint(model,
                Pg_bus - Pd + conv_P == 
                sum(Vm[i] * Vm[j] * (G[i,j] * cos(Va[i] - Va[j]) + B[i,j] * sin(Va[i] - Va[j]))
                    for j in 1:nac)
            )
        end
        
        # Q balance for PQ buses only
        if bus.bus_type == PQ_BUS
            @NLconstraint(model,
                Qg_bus - Qd + conv_Q ==
                sum(Vm[i] * Vm[j] * (G[i,j] * sin(Va[i] - Va[j]) - B[i,j] * cos(Va[i] - Va[j]))
                    for j in 1:nac)
            )
        end
    end
    
    # DC network equations (if present)
    if ndc > 1
        Gdc = sys.Gdc
        for d in 2:ndc
            dc_bus = sys.dc_buses[d]
            Pdc_load = dc_bus.pd_mw / baseMVA
            
            # DC power balance
            @NLconstraint(model,
                -Pdc_load == sum(Gdc[d,j] * Vdc[d] * Vdc[j] for j in 1:ndc)
            )
        end
    end
    
    # Minimize deviation from flat start (feasibility check)
    @objective(model, Min, 
        sum((Vm[i] - 1.0)^2 for i in 1:nac) + sum((Va[i])^2 for i in 1:nac)
    )
    
    verbose && println("  Solving with Ipopt...")
    start_time = time()
    optimize!(model)
    solve_time = (time() - start_time) * 1000
    
    status = termination_status(model)
    status_str = string(status)
    
    verbose && @printf("    Status: %s\n", status)
    verbose && @printf("    Solve time: %.1f ms\n", solve_time)
    
    is_feasible = status in [MOI.LOCALLY_SOLVED, MOI.OPTIMAL, MOI.ALMOST_LOCALLY_SOLVED]
    
    if is_feasible
        obj_val = objective_value(model)
        Vm_sol = [value(Vm[i]) for i in 1:nac]
        Va_sol = [value(Va[i]) for i in 1:nac]
        Vdc_sol = ndc > 0 ? [value(Vdc[d]) for d in 1:ndc] : Float64[]
        
        verbose && println("  ✅ FEASIBLE")
        verbose && @printf("    Objective: %.6f\n", obj_val)
        verbose && @printf("    Voltage range: [%.4f, %.4f] p.u.\n", minimum(Vm_sol), maximum(Vm_sol))
        
        return FeasibilityResult(
            true, status_str, obj_val, solve_time, 0,
            Vm_sol, Va_sol, Vdc_sol, copy(sys.Pg), copy(sys.Qg),
            total_load, total_gen, load_margin
        )
    else
        verbose && println("  ❌ INFEASIBLE or solver failed")
        return FeasibilityResult(
            false, status_str, 0.0, solve_time, 0,
            Float64[], Float64[], Float64[], Float64[], Float64[],
            total_load, total_gen, load_margin
        )
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
#  VALIDATION: COMPARE NEWTON-RAPHSON VS REFERENCE SOLVERS
# ═══════════════════════════════════════════════════════════════════════════════

"""
    validate_against_nlsolve(sys; verbose=false, tol=1e-4)

Validate our Newton-Raphson solution against NLsolve.
Returns (match::Bool, max_error::Float64, details::Dict)
"""
function HybridACDCPowerFlow.validate_against_nlsolve(sys::HybridSystem; verbose::Bool=false, tol::Float64=1e-4)
    verbose && println("\n" * "="^70)
    verbose && println("  VALIDATION: Newton-Raphson vs NLsolve")
    verbose && println("="^70)
    
    # Solve with our Newton-Raphson
    res_nr = solve_power_flow(sys)
    if !res_nr.converged
        verbose && println("  ❌ Newton-Raphson did not converge")
        return false, Inf, Dict("error" => "NR_NOT_CONVERGED")
    end
    
    # Solve with NLsolve
    res_nl = HybridACDCPowerFlow.check_power_flow_feasibility_nlsolve(sys; verbose=false)
    if !res_nl.feasible
        verbose && println("  ❌ NLsolve did not converge")
        return false, Inf, Dict("error" => "NLSOLVE_NOT_CONVERGED")
    end
    
    # Compare solutions
    Vm_err = maximum(abs.(res_nr.Vm .- res_nl.Vm))
    Va_err = maximum(abs.(res_nr.Va .- res_nl.Va))
    
    ndc = length(sys.dc_buses)
    Vdc_err = ndc > 0 ? maximum(abs.(res_nr.Vdc .- res_nl.Vdc)) : 0.0
    
    max_err = max(Vm_err, Va_err, Vdc_err)
    match = max_err < tol
    
    details = Dict(
        "Vm_max_error" => Vm_err,
        "Va_max_error" => Va_err,
        "Vdc_max_error" => Vdc_err,
        "NR_iterations" => res_nr.iterations,
        "NLsolve_iterations" => res_nl.iterations,
        "NR_residual" => norm(HybridACDCPowerFlow.PowerSystem.full_residual_simple(
            sys, res_nr.Vm, res_nr.Va, res_nr.Vdc), Inf),
        "NLsolve_residual" => res_nl.objective
    )
    
    if verbose
        @printf("    Vm max error: %.2e\n", Vm_err)
        @printf("    Va max error: %.2e\n", Va_err)
        @printf("    Vdc max error: %.2e\n", Vdc_err)
        @printf("    Overall: %s (tol=%.2e)\n", match ? "✅ MATCH" : "❌ MISMATCH", tol)
    end
    
    return match, max_err, details
end

"""
    validate_against_jump(sys; verbose=false, tol=1e-3)

Validate our Newton-Raphson solution against JuMP+Ipopt.
Returns (match::Bool, max_error::Float64, details::Dict)
"""
function HybridACDCPowerFlow.validate_against_jump(sys::HybridSystem; verbose::Bool=false, tol::Float64=1e-3)
    verbose && println("\n" * "="^70)
    verbose && println("  VALIDATION: Newton-Raphson vs JuMP+Ipopt")
    verbose && println("="^70)
    
    # Solve with our Newton-Raphson (measure wall time locally)
    nr_start = time_ns()
    res_nr = solve_power_flow(sys)
    nr_time_ms = (time_ns() - nr_start) / 1e6
    if !res_nr.converged
        verbose && println("  ❌ Newton-Raphson did not converge")
        return false, Inf, Dict("error" => "NR_NOT_CONVERGED")
    end
    
    # Solve with JuMP+Ipopt
    res_jump = HybridACDCPowerFlow.check_power_flow_feasibility_jump(sys; verbose=false)
    if !res_jump.feasible
        verbose && println("  ❌ JuMP+Ipopt did not converge")
        return false, Inf, Dict("error" => "JUMP_NOT_CONVERGED")
    end
    
    # Compare solutions
    Vm_err = maximum(abs.(res_nr.Vm .- res_jump.Vm))
    Va_err = maximum(abs.(res_nr.Va .- res_jump.Va))
    
    ndc = length(sys.dc_buses)
    Vdc_err = (ndc > 0 && length(res_jump.Vdc) > 0) ? 
              maximum(abs.(res_nr.Vdc .- res_jump.Vdc)) : 0.0
    
    max_err = max(Vm_err, Va_err, Vdc_err)
    match = max_err < tol
    
    details = Dict(
        "Vm_max_error" => Vm_err,
        "Va_max_error" => Va_err,
        "Vdc_max_error" => Vdc_err,
        "NR_time_ms" => nr_time_ms,
        "JuMP_time_ms" => res_jump.solve_time
    )
    
    if verbose
        @printf("    Vm max error: %.2e\n", Vm_err)
        @printf("    Va max error: %.2e\n", Va_err)
        @printf("    Vdc max error: %.2e\n", Vdc_err)
        @printf("    Overall: %s (tol=%.2e)\n", match ? "✅ MATCH" : "❌ MISMATCH", tol)
    end
    
    return match, max_err, details
end

# Extension methods are automatically available via HybridACDCPowerFlow module

end  # module FeasibilityExt
