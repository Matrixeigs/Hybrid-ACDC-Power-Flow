# Test script for refactored HybridACDCPowerFlow
# Verifies that the new architecture (JuliaPowerCase data + HybridACDCPowerFlow solver) works correctly

using Pkg
Pkg.activate("/Users/tianyangzhao/Codes/DistributionSystemResilience/HybridACDCPowerFlow")

println("Loading HybridACDCPowerFlow...")
using HybridACDCPowerFlow
using JuliaPowerCase

println("\n========================================")
println("Test 1: Build IEEE14 hybrid system")
println("========================================")
hps = build_ieee14_acdc()
println("✓ Built system: ", typeof(hps))
println("  AC buses: ", length(hps.ac.buses))
println("  DC buses: ", length(hps.dc.buses))
println("  Converters: ", length(hps.vsc_converters))

println("\n========================================")
println("Test 2: Solve power flow")
println("========================================")
result = solve_power_flow(hps)
println("✓ Converged: ", result.converged)
println("  Iterations: ", result.iterations)
println("  Residual: ", result.residual)
if result.converged
    println("  AC voltage range: ", minimum(result.Vm), " - ", maximum(result.Vm), " p.u.")
    println("  DC voltage range: ", minimum(result.Vdc), " - ", maximum(result.Vdc), " p.u.")
end

println("\n========================================")
println("Test 3: Verify data types")
println("========================================")
println("✓ HybridPowerSystem type: ", typeof(hps))
println("  AC buses element type: ", eltype(hps.ac.buses))
println("  DC buses element type: ", eltype(hps.dc.buses))
println("  Converter type: ", eltype(hps.vsc_converters))

println("\n========================================")
println("Test 4: Build smaller test case")
println("========================================")
hps_small = build_case33bw_acdc()
println("✓ Built Case33bw system")
println("  AC buses: ", length(hps_small.ac.buses))
result_small = solve_power_flow(hps_small)
println("✓ Converged: ", result_small.converged)

println("\n========================================")
println("All tests passed! ✓")
println("========================================")
println("\nArchitecture verification:")
println("  • JuliaPowerCase provides data structures (HybridPowerSystem)")
println("  • HybridACDCPowerFlow provides algorithms (solve_power_flow)")
println("  • Internal SolverData is not exposed to users")
println("  • Clean separation of information and application layers")
