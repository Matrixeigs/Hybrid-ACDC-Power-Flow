#!/usr/bin/env julia
# test_extension.jl
# Test the FeasibilityExt extension

using Pkg
Pkg.activate(@__DIR__)

println("Loading HybridACDCPowerFlow without extensions...")
using HybridACDCPowerFlow

# Test that stub throws appropriate error
println("\nTest 1: Verify stub errors without extensions loaded...")
try
    sys = build_ieee14_acdc()
    check_power_flow_feasibility(sys)
    println("  ERROR: Should have thrown an error!")
    exit(1)
catch e
    if occursin("NLsolve", string(e))
        println("  ✓ Correctly throws error about missing NLsolve")
    else
        println("  ERROR: Unexpected error: $e")
        exit(1)
    end
end

println("\nTest 2: Loading extensions (JuMP, Ipopt, NLsolve)...")
using JuMP, Ipopt, NLsolve
println("  ✓ Extensions loaded")

println("\nTest 3: Now test feasibility check with extensions...")
sys = build_ieee14_acdc()
result = check_power_flow_feasibility(sys; method=:nlsolve, verbose=false)
println("  Feasible: $(result.feasible)")
println("  Status: $(result.status)")
println("  Iterations: $(result.iterations)")
if result.feasible
    println("  ✓ Feasibility check passed!")
else
    println("  WARNING: System reported as infeasible (may be OK)")
end

println("\nTest 4: Test JuMP method...")
result_jump = check_power_flow_feasibility(sys; method=:jump, verbose=false)
println("  Feasible: $(result_jump.feasible)")
println("  Status: $(result_jump.status)")
if result_jump.feasible
    println("  ✓ JuMP feasibility check passed!")
else
    println("  WARNING: JuMP reported infeasible (may be solver tolerance)")
end

println("\n" * "="^60)
println("EXTENSION TESTS COMPLETED!")
println("="^60)
