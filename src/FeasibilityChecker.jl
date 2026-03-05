"""
    FeasibilityChecker

Feasibility checking interface for Hybrid AC/DC Power Flow.

The actual implementations are provided by the FeasibilityExt extension
when JuMP, Ipopt, and NLsolve are loaded. This file provides:
1. The FeasibilityResult struct definition
2. Stub functions that error when extensions aren't loaded
3. The main dispatch function check_power_flow_feasibility()

To use feasibility checking, make sure to have JuMP, Ipopt, and NLsolve
in your environment:
    using HybridACDCPowerFlow
    using JuMP, Ipopt, NLsolve  # This triggers extension loading
    
    result = check_power_flow_feasibility(sys)

Author: HybridACDCPowerFlow Development Team
Version: 0.4.0
"""

# This file is included within the HybridACDCPowerFlow module

"""
Feasibility check result structure
"""
struct FeasibilityResult
    feasible::Bool
    status::String
    objective::Float64
    solve_time::Float64
    iterations::Int
    
    # Solution (if feasible)
    Vm::Vector{Float64}
    Va::Vector{Float64}
    Vdc::Vector{Float64}
    Pg::Vector{Float64}
    Qg::Vector{Float64}
    
    # Diagnostic info
    total_load::Float64
    total_generation_capacity::Float64
    load_margin::Float64
    
    FeasibilityResult(feasible, status, obj, time, iters, Vm, Va, Vdc, Pg, Qg, load, cap, margin) =
        new(feasible, status, obj, time, iters, Vm, Va, Vdc, Pg, Qg, load, cap, margin)
end

# Stub functions that will be overridden by the extension
function check_power_flow_feasibility_nlsolve(sys; kwargs...)
    error("""
        check_power_flow_feasibility_nlsolve requires NLsolve.
        Please load NLsolve first:
            using NLsolve
        Then the extension will be automatically loaded.
        """)
end

function check_power_flow_feasibility_jump(sys; kwargs...)
    error("""
        check_power_flow_feasibility_jump requires JuMP and Ipopt.
        Please load them first:
            using JuMP, Ipopt
        Then the extension will be automatically loaded.
        """)
end

"""
    check_power_flow_feasibility(sys::HybridSystem; method=:nlsolve, verbose=false, kwargs...)

Main interface for feasibility checking. Dispatches to appropriate method:
- `:nlsolve` - NLsolve root-finding (default, recommended)
  Accepts: verbose, max_iter, tol
- `:jump` - JuMP+Ipopt optimization (fallback)
  Accepts: verbose, max_time

Requires the FeasibilityExt extension to be loaded (via using JuMP, Ipopt, NLsolve).
"""
function check_power_flow_feasibility(sys; method::Symbol=:nlsolve, verbose::Bool=false, 
                                      max_iter::Int=50, max_time::Float64=10.0, tol::Float64=1e-6)
    if method == :nlsolve
        return check_power_flow_feasibility_nlsolve(sys; verbose=verbose, max_iter=max_iter, tol=tol)
    elseif method == :jump
        return check_power_flow_feasibility_jump(sys; verbose=verbose, max_time=max_time)
    else
        error("Unknown feasibility check method: $method (use :nlsolve or :jump)")
    end
end

# Validation function stubs (implemented in FeasibilityExt)
function validate_against_nlsolve(sys; kwargs...)
    error("""
        validate_against_nlsolve requires NLsolve.
        Please load NLsolve first:
            using NLsolve
        Then the extension will be automatically loaded.
        """)
end

function validate_against_jump(sys; kwargs...)
    error("""
        validate_against_jump requires JuMP and Ipopt.
        Please load them first:
            using JuMP, Ipopt
        Then the extension will be automatically loaded.
        """)
end
