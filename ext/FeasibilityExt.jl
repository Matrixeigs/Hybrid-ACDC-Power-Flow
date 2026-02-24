"""
    FeasibilityExt

Package extension for HybridACDCPowerFlow that provides feasibility checking
when NLsolve, JuMP, and Ipopt are loaded.
"""
module FeasibilityExt

using HybridACDCPowerFlow
using HybridACDCPowerFlow: HybridSystem, BusType, PQ, PV, SLACK, 
                           ConverterMode, PQ_MODE, VDC_Q, VDC_VAC,
                           power_flow_residual, FeasibilityResult

using NLsolve
using JuMP
using Ipopt
using Printf
using LinearAlgebra
using SparseArrays

# ═══════════════════════════════════════════════════════════════════════════════
#  NLSOLVE-BASED FEASIBILITY CHECK
# ═══════════════════════════════════════════════════════════════════════════════

"""
    check_power_flow_feasibility_nlsolve(sys::HybridSystem; verbose::Bool=false, max_iter::Int=50)

Check power flow feasibility using NLsolve root-finding on actual residual equations.

This is the PREFERRED method because it:
- Uses the exact same `power_flow_residual` function as the actual solver
- No risk of formulation errors (directly solves f(x)=0)
- Appropriate tool for root-finding problems
- Robust trust-region methods

Returns a `FeasibilityResult` indicating whether a solution exists.
"""
function HybridACDCPowerFlow.check_power_flow_feasibility_nlsolve(sys::HybridSystem; 
                                                                  verbose::Bool=false, 
                                                                  max_iter::Int=50, 
                                                                  tol::Float64=1e-6)
    nac = length(sys.ac_buses)
    ndc = length(sys.dc_buses)
    
    # Compute load/generation statistics
    total_load = sum(bus.Pd for bus in sys.ac_buses)
    total_gen_max = sum(bus.Pg for bus in sys.ac_buses if bus.type in (PV, SLACK); init=0.0)
    if total_gen_max < 1e-6
        total_gen_max = 1.5 * total_load
    end
    load_margin = total_gen_max - total_load
    
    verbose && println("\n" * "="^70)
    verbose && println("  NLSOLVE FEASIBILITY CHECK (Root-Finding)")
    verbose && println("="^70)
    verbose && @printf("    Total load: %.4f, Gen: %.4f, Margin: %.4f p.u.\n", 
                      total_load, total_gen_max, load_margin)
    
    # Quick pre-check
    if load_margin < -0.01
        verbose && println("  ❌ INFEASIBLE: Insufficient generation")
        return FeasibilityResult(
            false, "INSUFFICIENT_GENERATION", 0.0, 0.0, 0,
            zeros(nac), zeros(nac), zeros(ndc), zeros(nac), zeros(nac),
            total_load, total_gen_max, load_margin
        )
    end
    
    # Identify bus types for state vector construction
    pq_idx = findall(b -> b.type == PQ, sys.ac_buses)
    pv_idx = findall(b -> b.type == PV, sys.ac_buses)
    slack_idx = findall(b -> b.type == SLACK, sys.ac_buses)
    
    # Initialize state: Va (all non-slack), Vm (PQ only), Vdc (non-slack DC)
    Vm_init = ones(nac)
    Va_init = zeros(nac)
    Vdc_init = ndc > 0 ? [sys.dc_buses[d].Vdc_set for d in 1:ndc] : Float64[]
    
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
    x0 .+= 0.01 * randn(length(x0))
    
    # Define residual function wrapper
    function residual_wrapper!(F, x)
        Vm = ones(nac)
        Va = zeros(nac)
        Vdc = ndc > 0 ? ones(ndc) : Float64[]
        
        Va[non_slack_ac] .= x[1:n_Va]
        Vm[pq_idx] .= x[n_Va+1:n_Va+n_Vm]
        if n_Vdc > 0
            Vdc[2:end] .= x[n_Va+n_Vm+1:end]
        end
        
        for i in pv_idx
            Vm[i] = sys.ac_buses[i].Vm
        end
        for i in slack_idx
            Vm[i] = sys.ac_buses[i].Vm
            Va[i] = sys.ac_buses[i].Va
        end
        if ndc > 0
            Vdc[1] = sys.dc_buses[1].Vdc_set
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
            Vm = ones(nac)
            Va = zeros(nac)
            Vdc = ndc > 0 ? ones(ndc) : Float64[]
            
            x_sol = result.zero
            Va[non_slack_ac] .= x_sol[1:n_Va]
            Vm[pq_idx] .= x_sol[n_Va+1:n_Va+n_Vm]
            if n_Vdc > 0
                Vdc[2:end] .= x_sol[n_Va+n_Vm+1:end]
            end
            
            for i in pv_idx
                Vm[i] = sys.ac_buses[i].Vm
            end
            for i in slack_idx
                Vm[i] = sys.ac_buses[i].Vm
                Va[i] = sys.ac_buses[i].Va
            end
            if ndc > 0
                Vdc[1] = sys.dc_buses[1].Vdc_set
            end
            
            F_ac, F_dc = power_flow_residual(sys, Vm, Va, Vdc)
            F_full = ndc > 0 ? vcat(F_ac, F_dc) : F_ac
            residual_norm = norm(F_full, Inf)
            verbose && @printf("  ✅ FEASIBLE: Converged in %d iterations (residual: %.3e)\n", 
                              result.iterations, residual_norm)
            
            Pg = [sys.ac_buses[i].Pg for i in 1:nac]
            Qg = [sys.ac_buses[i].Qg for i in 1:nac]
            
            return FeasibilityResult(
                true, "CONVERGED", residual_norm, solve_time, result.iterations,
                Vm, Va, Vdc, Pg, Qg,
                total_load, total_gen_max, load_margin
            )
        else
            F_ac, F_dc = power_flow_residual(sys, ones(nac), zeros(nac), ndc > 0 ? ones(ndc) : Float64[])
            F_full = ndc > 0 ? vcat(F_ac, F_dc) : F_ac
            residual_norm = norm(F_full, Inf)
            verbose && @printf("  ❌ INFEASIBLE: No convergence (residual: %.3e)\n", residual_norm)
            return FeasibilityResult(
                false, "NO_CONVERGENCE", residual_norm, solve_time, result.iterations,
                zeros(nac), zeros(nac), zeros(ndc), zeros(nac), zeros(nac),
                total_load, total_gen_max, load_margin
            )
        end
        
    catch e
        solve_time = (time() - start_time) * 1000
        verbose && println("  ❌ INFEASIBLE: Solver error - $e")
        return FeasibilityResult(
            false, "SOLVER_ERROR", 0.0, solve_time, 0,
            zeros(nac), zeros(nac), zeros(ndc), zeros(nac), zeros(nac),
            total_load, total_gen_max, load_margin
        )
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
#  JUMP/IPOPT-BASED FEASIBILITY CHECK
# ═══════════════════════════════════════════════════════════════════════════════

"""
    check_power_flow_feasibility_jump(sys::HybridSystem; verbose::Bool=false, max_time::Float64=10.0)

Check power flow feasibility using JuMP+Ipopt optimization.
"""
function HybridACDCPowerFlow.check_power_flow_feasibility_jump(sys::HybridSystem; 
                                                               verbose::Bool=false, 
                                                               max_time::Float64=10.0)
    nac = length(sys.ac_buses)
    ndc = length(sys.dc_buses)
    nc = length(sys.converters)
    
    total_load = sum(bus.Pd for bus in sys.ac_buses)
    total_gen_max = sum(bus.Pg for bus in sys.ac_buses if bus.type in (PV, SLACK); init=0.0)
    if total_gen_max < 1e-6
        total_gen_max = 1.5 * total_load
    end
    load_margin = total_gen_max - total_load
    
    verbose && println("\n" * "="^70)
    verbose && println("  JuMP FEASIBILITY CHECK")
    verbose && println("="^70)
    verbose && @printf("    AC buses: %d, DC buses: %d, Converters: %d\n", nac, ndc, nc)
    verbose && @printf("    Total load: %.4f p.u.\n", total_load)
    verbose && @printf("    Gen capacity: %.4f p.u.\n", total_gen_max)
    verbose && @printf("    Load margin: %.4f p.u. (%.1f%%)\n", load_margin, 100*load_margin/(total_load+1e-10))
    
    if load_margin < -0.01
        verbose && println("  ❌ INFEASIBLE: Insufficient generation capacity")
        return FeasibilityResult(
            false, "INSUFFICIENT_GENERATION", 0.0, 0.0, 0,
            Float64[], Float64[], Float64[], Float64[], Float64[],
            total_load, total_gen_max, load_margin
        )
    end
    
    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "print_level", verbose ? 5 : 0)
    set_optimizer_attribute(model, "max_iter", 1000)
    set_optimizer_attribute(model, "max_cpu_time", max_time)
    set_optimizer_attribute(model, "tol", 1e-4)
    set_optimizer_attribute(model, "constr_viol_tol", 1e-3)
    set_optimizer_attribute(model, "acceptable_tol", 1e-3)
    set_optimizer_attribute(model, "acceptable_iter", 15)
    
    @variable(model, 0.85 <= Vm[i=1:nac] <= 1.15, start=1.0)
    @variable(model, -π <= Va[i=1:nac] <= π, start=0.0)
    
    if ndc > 0
        @variable(model, 0.85 <= Vdc[d=1:ndc] <= 1.15, start=1.0)
        fix(Vdc[1], sys.dc_buses[1].Vdc_set; force=true)
    end
    
    @variable(model, Pg[i=1:nac])
    @variable(model, Qg[i=1:nac])
    
    if nc > 0
        @variable(model, Pac[c=1:nc], start=0.0)
        @variable(model, Qac[c=1:nc], start=0.0)
        
        for (c, conv) in enumerate(sys.converters)
            if conv.status
                set_lower_bound(Pac[c], -conv.Smax)
                set_upper_bound(Pac[c], conv.Smax)
                set_lower_bound(Qac[c], -conv.Smax)
                set_upper_bound(Qac[c], conv.Smax)
            else
                fix(Pac[c], 0.0; force=true)
                fix(Qac[c], 0.0; force=true)
            end
        end
    end
    
    for i in 1:nac
        bus = sys.ac_buses[i]
        if bus.type == SLACK
            set_start_value(Pg[i], total_load * 0.5)
            set_start_value(Qg[i], 0.0)
            set_lower_bound(Pg[i], 0.0)
            set_upper_bound(Pg[i], total_gen_max)
            set_lower_bound(Qg[i], -total_gen_max)
            set_upper_bound(Qg[i], total_gen_max)
        elseif bus.type == PV
            set_start_value(Pg[i], bus.Pg)
            set_start_value(Qg[i], 0.0)
            set_lower_bound(Pg[i], 0.0)
            set_upper_bound(Pg[i], max(bus.Pg * 2.0, 1.0))
            set_lower_bound(Qg[i], -1.0)
            set_upper_bound(Qg[i], 1.0)
        else
            fix(Pg[i], 0.0; force=true)
            fix(Qg[i], 0.0; force=true)
        end
    end
    
    ref_bus = findfirst(b -> b.type == SLACK, sys.ac_buses)
    if ref_bus !== nothing
        fix(Va[ref_bus], 0.0; force=true)
    else
        fix(Va[1], 0.0; force=true)
    end
    
    for i in 1:nac
        if sys.ac_buses[i].type in (PV, SLACK)
            set_lower_bound(Vm[i], sys.ac_buses[i].Vm)
            set_upper_bound(Vm[i], sys.ac_buses[i].Vm)
        end
    end
    
    G = real.(Matrix(sys.Ybus))
    B = imag.(Matrix(sys.Ybus))
    
    for i in 1:nac
        bus = sys.ac_buses[i]
        convs_at_i = [c for c in 1:nc if sys.converters[c].ac_bus == i && sys.converters[c].status]
        
        if nc > 0 && !isempty(convs_at_i)
            @NLconstraint(model,
                Pg[i] - bus.Pd + sum(Pac[c] for c in convs_at_i) == 
                sum(Vm[i] * Vm[j] * (G[i,j] * cos(Va[i] - Va[j]) + B[i,j] * sin(Va[i] - Va[j]))
                    for j in 1:nac)
            )
        else
            @NLconstraint(model,
                Pg[i] - bus.Pd == 
                sum(Vm[i] * Vm[j] * (G[i,j] * cos(Va[i] - Va[j]) + B[i,j] * sin(Va[i] - Va[j]))
                    for j in 1:nac)
            )
        end
        
        if nc > 0 && !isempty(convs_at_i)
            @NLconstraint(model,
                Qg[i] - bus.Qd + sum(Qac[c] for c in convs_at_i) ==
                sum(Vm[i] * Vm[j] * (G[i,j] * sin(Va[i] - Va[j]) - B[i,j] * cos(Va[i] - Va[j]))
                    for j in 1:nac)
            )
        else
            @NLconstraint(model,
                Qg[i] - bus.Qd ==
                sum(Vm[i] * Vm[j] * (G[i,j] * sin(Va[i] - Va[j]) - B[i,j] * cos(Va[i] - Va[j]))
                    for j in 1:nac)
            )
        end
    end
    
    if nc > 0
        for c in 1:nc
            conv = sys.converters[c]
            conv.status || continue
            
            if conv.mode == PQ_MODE
                fix(Pac[c], conv.Pset; force=true)
                fix(Qac[c], conv.Qset; force=true)
            elseif conv.mode == VDC_Q
                fix(Qac[c], conv.Qset; force=true)
            end
        end
    end
    
    if ndc > 0 && nc > 0
        Gdc = sys.Gdc
        
        for d in 2:ndc
            dc_bus = sys.dc_buses[d]
            convs_at_d = [c for c in 1:nc if sys.converters[c].dc_bus == d && sys.converters[c].status]
            
            @NLconstraint(model,
                sum(-Pac[c] * 0.98 for c in convs_at_d; init=0.0) ==
                sum(Gdc[d,j] * Vdc[d] * Vdc[j] for j in 1:ndc) - dc_bus.Pdc
            )
        end
    end
    
    @objective(model, Min, 
        sum((Pg[i] - sys.ac_buses[i].Pg)^2 + (Qg[i] - sys.ac_buses[i].Qg)^2 for i in 1:nac)
    )
    
    verbose && println("\n  Solving with Ipopt...")
    start_time = time()
    optimize!(model)
    solve_time = time() - start_time
    
    status = termination_status(model)
    status_str = string(status)
    iterations = 0
    
    verbose && println("\n  Results:")
    verbose && @printf("    Status: %s\n", status)
    verbose && @printf("    Primal: %s\n", primal_status(model))
    verbose && @printf("    Solve time: %.3f s\n", solve_time)
    
    is_feasible = status in [MOI.LOCALLY_SOLVED, MOI.OPTIMAL, MOI.ALMOST_LOCALLY_SOLVED]
    
    if is_feasible
        obj_val = objective_value(model)
        Vm_sol = [value(Vm[i]) for i in 1:nac]
        Va_sol = [value(Va[i]) for i in 1:nac]
        Vdc_sol = ndc > 1 ? [value(Vdc[d]) for d in 1:ndc] : Float64[]
        Pg_sol = [value(Pg[i]) for i in 1:nac]
        Qg_sol = [value(Qg[i]) for i in 1:nac]
        
        verbose && println("  ✓ PROBLEM IS FEASIBLE")
        verbose && @printf("    Objective: %.6f\n", obj_val)
        verbose && @printf("    Voltage range: [%.4f, %.4f] p.u.\n", minimum(Vm_sol), maximum(Vm_sol))
        verbose && @printf("    Total generation: %.4f p.u.\n", sum(Pg_sol))
        
        return FeasibilityResult(
            true, status_str, obj_val, solve_time, iterations,
            Vm_sol, Va_sol, Vdc_sol, Pg_sol, Qg_sol,
            total_load, total_gen_max, load_margin
        )
    else
        verbose && println("  ❌ PROBLEM MAY BE INFEASIBLE")
        verbose && @printf("    Status: %s\n", status)
        
        return FeasibilityResult(
            false, status_str, 0.0, solve_time, iterations,
            Float64[], Float64[], Float64[], Float64[], Float64[],
            total_load, total_gen_max, load_margin
        )
    end
end

end  # module FeasibilityExt
