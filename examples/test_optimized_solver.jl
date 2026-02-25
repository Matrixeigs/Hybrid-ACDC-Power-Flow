#!/usr/bin/env julia
# test_optimized_solver.jl
# Test the optimized sparse Jacobian solver

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

println("Loading HybridACDCPowerFlow...")
using HybridACDCPowerFlow
println("✓ Module loaded successfully\n")

# Test case 1: IEEE14 AC/DC
println("="^60)
println("TEST 1: IEEE14 AC/DC System")
println("="^60)
sys14 = build_ieee14_acdc()
result14 = solve_power_flow(sys14)
println("  Converged: $(result14.converged)")
println("  Iterations: $(result14.iterations)")
println("  Residual: $(result14.residual)")
println("  Voltage range: [$(minimum(result14.Vm)), $(maximum(result14.Vm))]")
@assert result14.converged "IEEE14 should converge"
println("  ✓ PASSED\n")

# Test case 2: case33bw (distribution)
println("="^60)
println("TEST 2: case33bw Distribution System")
println("="^60)
sys33bw = build_case33bw_acdc()
result33bw = solve_power_flow(sys33bw)
println("  Converged: $(result33bw.converged)")
println("  Iterations: $(result33bw.iterations)")
println("  Residual: $(result33bw.residual)")
println("  Voltage range: [$(minimum(result33bw.Vm)), $(maximum(result33bw.Vm))]")
@assert result33bw.converged "case33bw should converge"
println("  ✓ PASSED\n")

# Test case 3: case69 (larger distribution)
println("="^60)
println("TEST 3: case69 Distribution System")
println("="^60)
sys69 = build_case69_acdc()
result69 = solve_power_flow(sys69)
println("  Converged: $(result69.converged)")
println("  Iterations: $(result69.iterations)")
println("  Residual: $(result69.residual)")
println("  Voltage range: [$(minimum(result69.Vm)), $(maximum(result69.Vm))]")
@assert result69.converged "case69 should converge"
println("  ✓ PASSED\n")

# Test case 4: case300 (transmission)
println("="^60)
println("TEST 4: case300 Transmission System")
println("="^60)
sys300 = build_case300_acdc()
result300 = solve_power_flow(sys300)
println("  Converged: $(result300.converged)")
println("  Iterations: $(result300.iterations)")
println("  Residual: $(result300.residual)")
println("  Voltage range: [$(minimum(result300.Vm)), $(maximum(result300.Vm))]")
@assert result300.converged "case300 should converge"
println("  ✓ PASSED\n")

# Test case 5: IEEE118 (larger transmission)
println("="^60)
println("TEST 5: IEEE118 AC/DC System")
println("="^60)
sys118 = build_ieee118_acdc()
result118 = solve_power_flow(sys118)
println("  Converged: $(result118.converged)")
println("  Iterations: $(result118.iterations)")
println("  Residual: $(result118.residual)")
println("  Voltage range: [$(minimum(result118.Vm)), $(maximum(result118.Vm))]")
@assert result118.converged "IEEE118 should converge"
println("  ✓ PASSED\n")

# Summary
println("="^60)
println("ALL TESTS PASSED!")
println("="^60)
println("\nSummary:")
println("  IEEE14:   $(result14.iterations) iters, res=$(result14.residual)")
println("  case33bw: $(result33bw.iterations) iters, res=$(result33bw.residual)")
println("  case69:   $(result69.iterations) iters, res=$(result69.residual)")
println("  case300:  $(result300.iterations) iters, res=$(result300.residual)")
println("  IEEE118:  $(result118.iterations) iters, res=$(result118.residual)")
