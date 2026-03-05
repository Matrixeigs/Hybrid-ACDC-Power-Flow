#!/usr/bin/env julia
# benchmark_optimizations.jl
# Benchmark the optimized sparse Jacobian solver

using Pkg
Pkg.activate(@__DIR__)

using HybridACDCPowerFlow
using BenchmarkTools
using Printf

println("="^70)
println("BENCHMARK: Optimized Sparse Jacobian Newton-Raphson Solver")
println("="^70)
println()

# Build test systems
systems = [
    ("IEEE14 AC/DC", build_ieee14_acdc()),
    ("case33bw", build_case33bw_acdc()),
    ("case69", build_case69_acdc()),
    ("case300", build_case300_acdc()),
    ("IEEE118 AC/DC", build_ieee118_acdc()),
]

println("Warming up...")
for (name, sys) in systems
    solve_power_flow(sys)  # Warm up
end


println("\nBenchmarking power flow solve times:\n")
results = []

for (name, sys) in systems
    # Benchmark
    b = @benchmark solve_power_flow($sys) seconds=5
    
    # Get result for accuracy check
    result = solve_power_flow(sys)
    
    # Statistics
    median_time = median(b.times) / 1e6  # Convert to ms
    min_time = minimum(b.times) / 1e6
    allocs = b.allocs
    memory = b.memory / 1024  # KB
    
    push!(results, (name=name, n_bus=length(sys.ac_buses), median_ms=median_time, 
                    min_ms=min_time, allocs=allocs, memory_kb=memory,
                    iters=result.iterations, residual=result.residual))
    
    @printf("%-15s | %3d buses | %7.3f ms (median) | %5d allocs | %8.1f KB | %d iters\n",
            name, length(sys.ac_buses), median_time, allocs, memory, result.iterations)
end

println("\n" * "="^70)
println("BENCHMARK COMPLETE")
println("="^70)

# Print detailed summary
println("\nDetailed Results:")
println("-"^70)
for r in results
    @printf("%-15s: %.3f ms median, %.3f ms min, %d allocations, %.1f KB\n",
            r.name, r.median_ms, r.min_ms, r.allocs, r.memory_kb)
    @printf("               Converged in %d iterations, residual: %.2e\n",
            r.iters, r.residual)
end

# Performance metrics
println("\nPerformance Summary:")
println("-"^70)
for r in results
    time_per_bus = r.median_ms / r.n_bus * 1000  # μs per bus
    @printf("%-15s: %.2f μs/bus, %.0f bytes/bus\n",
            r.name, time_per_bus, r.memory_kb * 1024 / r.n_bus)
end
