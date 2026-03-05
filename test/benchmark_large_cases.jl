"""
Large-Scale MATPOWER Cases Benchmark

Tests HybridACDCPowerFlow Newton-Raphson solver performance across:
- case33bw (33 buses)
- case33mg (33 buses with microgrid)
- case69 (69 buses)
- case300 (300 buses)
- case_ACTIVSg2000 (2000 buses)

Compares against NLsolve and JuMP+Ipopt when available.
"""

using NLsolve
using JuMP
using Ipopt
using HybridACDCPowerFlow
using Printf
using Statistics

# =====================================================================
# Benchmark Configuration
# =====================================================================
const NUM_WARMUP_RUNS = 2
const NUM_TIMED_RUNS = 5
const COMPARE_REFERENCE_SOLVERS = true  # Set false to skip NLsolve/JuMP comparison

# =====================================================================
# Test Case Definitions
# =====================================================================
const TEST_CASES = [
    ("case33bw", HybridACDCPowerFlow.TestSystems.build_case33bw_acdc),
    ("case33mg", HybridACDCPowerFlow.TestSystems.build_case33mg_acdc),
    ("case69", HybridACDCPowerFlow.TestSystems.build_case69_acdc),
    ("case300", HybridACDCPowerFlow.TestSystems.build_case300_acdc),
    ("case2000", HybridACDCPowerFlow.TestSystems.build_case2000_acdc),
]

# =====================================================================
# Benchmark Functions
# =====================================================================

function benchmark_newton_raphson(sys; num_warmup=NUM_WARMUP_RUNS, num_runs=NUM_TIMED_RUNS)
    # Warmup
    for _ in 1:num_warmup
        solve_power_flow(sys)
    end
    
    # Timed runs
    times = Float64[]
    converged = true
    result = nothing
    
    for _ in 1:num_runs
        t_start = time_ns()
        result = solve_power_flow(sys)
        t_end = time_ns()
        push!(times, (t_end - t_start) / 1e6)  # Convert to ms
        converged &= result.converged
    end
    
    return (
        times = times,
        mean_time = mean(times),
        std_time = std(times),
        min_time = minimum(times),
        max_time = maximum(times),
        converged = converged,
        iterations = result.iterations,
        result = result
    )
end

function benchmark_nlsolve(sys; num_warmup=NUM_WARMUP_RUNS, num_runs=NUM_TIMED_RUNS)
    # Warmup
    for _ in 1:num_warmup
        try
            check_power_flow_feasibility_nlsolve(sys)
        catch
            return nothing
        end
    end
    
    # Timed runs
    times = Float64[]
    converged = true
    result = nothing
    
    for _ in 1:num_runs
        t_start = time_ns()
        result = check_power_flow_feasibility_nlsolve(sys)
        t_end = time_ns()
        push!(times, (t_end - t_start) / 1e6)
        converged &= result.feasible
    end
    
    return (
        times = times,
        mean_time = mean(times),
        std_time = std(times),
        min_time = minimum(times),
        max_time = maximum(times),
        converged = converged,
        iterations = result.iterations,
        result = (Vm = result.Vm, Va = result.Va, converged = result.feasible)
    )
end

function benchmark_jump_ipopt(sys; num_warmup=NUM_WARMUP_RUNS, num_runs=NUM_TIMED_RUNS)
    # Warmup
    for _ in 1:num_warmup
        try
            check_power_flow_feasibility_jump(sys)
        catch e
            @warn "JuMP+Ipopt failed" exception=e
            return nothing
        end
    end
    
    # Timed runs
    times = Float64[]
    converged = true
    result = nothing
    
    for _ in 1:num_runs
        t_start = time_ns()
        result = check_power_flow_feasibility_jump(sys)
        t_end = time_ns()
        push!(times, (t_end - t_start) / 1e6)
        converged &= result.feasible
    end
    
    return (
        times = times,
        mean_time = mean(times),
        std_time = std(times),
        min_time = minimum(times),
        max_time = maximum(times),
        converged = converged,
        iterations = result.iterations,
        result = (Vm = result.Vm, Va = result.Va, converged = result.feasible)
    )
end

# =====================================================================
# System Size Analysis
# =====================================================================
function analyze_system(sys)
    n_ac_buses = length(sys.ac_buses)
    n_ac_branches = length(sys.ac_branches)
    n_dc_buses = length(sys.dc_buses)
    n_dc_branches = length(sys.dc_branches)
    n_converters = length(sys.converters)
    n_generators = length(sys.generators)
    
    # Count bus types
    n_pq = count(b -> b.bus_type == HybridACDCPowerFlow.PQ, sys.ac_buses)
    n_pv = count(b -> b.bus_type == HybridACDCPowerFlow.PV, sys.ac_buses)
    n_slack = count(b -> b.bus_type == HybridACDCPowerFlow.SLACK, sys.ac_buses)
    
    # Calculate problem size (unknowns)
    # AC: 2*n_pq + n_pv (Va for all, Vm for PQ)
    # DC: n_dc - 1 (one slack)
    n_unknowns = 2*n_pq + n_pv + max(0, n_dc_buses - 1)
    
    return (
        n_ac_buses = n_ac_buses,
        n_ac_branches = n_ac_branches,
        n_dc_buses = n_dc_buses,
        n_dc_branches = n_dc_branches,
        n_converters = n_converters,
        n_generators = n_generators,
        n_pq = n_pq,
        n_pv = n_pv,
        n_slack = n_slack,
        n_unknowns = n_unknowns
    )
end

# =====================================================================
# Result Validation
# =====================================================================
function validate_results(nr_result, ref_result, ref_name)
    if isnothing(ref_result) || !nr_result.converged || !ref_result.converged
        return nothing
    end
    
    vm_nr = nr_result.result.Vm
    va_nr = nr_result.result.Va
    vm_ref = ref_result.result.Vm
    va_ref = ref_result.result.Va
    
    vm_error = maximum(abs.(vm_nr - vm_ref))
    va_error = maximum(abs.(va_nr - va_ref))
    
    return (
        vm_max_error = vm_error,
        va_max_error = va_error,
        vm_rmse = sqrt(mean((vm_nr - vm_ref).^2)),
        va_rmse = sqrt(mean((va_nr - va_ref).^2))
    )
end

# =====================================================================
# Main Benchmark Runner
# =====================================================================
function run_benchmarks()
    println("=" ^ 80)
    println("HybridACDCPowerFlow Large-Scale Benchmark")
    println("=" ^ 80)
    println()
    
    results = Dict{String, Any}()
    
    for (name, builder) in TEST_CASES
        println("-" ^ 80)
        println("Testing: $name")
        println("-" ^ 80)
        
        # Build system
        print("Building system... ")
        t_build = @elapsed sys = builder()
        println("done ($(@sprintf("%.2f", t_build*1000)) ms)")
        
        # Analyze system
        info = analyze_system(sys)
        println()
        println("System Statistics:")
        println("  AC Buses:     $(info.n_ac_buses) (PQ=$(info.n_pq), PV=$(info.n_pv), Slack=$(info.n_slack))")
        println("  AC Branches:  $(info.n_ac_branches)")
        println("  DC Buses:     $(info.n_dc_buses)")
        println("  DC Branches:  $(info.n_dc_branches)")
        println("  Converters:   $(info.n_converters)")
        println("  Generators:   $(info.n_generators)")
        println("  Unknowns:     $(info.n_unknowns)")
        println()
        
        # Newton-Raphson benchmark
        print("Newton-Raphson (warmup=$NUM_WARMUP_RUNS, runs=$NUM_TIMED_RUNS)... ")
        nr_result = benchmark_newton_raphson(sys)
        println("done")
        println("  Time: $(@sprintf("%.3f", nr_result.mean_time)) ± $(@sprintf("%.3f", nr_result.std_time)) ms")
        println("  Range: [$(@sprintf("%.3f", nr_result.min_time)), $(@sprintf("%.3f", nr_result.max_time))] ms")
        println("  Converged: $(nr_result.converged), Iterations: $(nr_result.iterations)")
        
        # Reference solver benchmarks
        nlsolve_result = nothing
        jump_result = nothing
        
        if COMPARE_REFERENCE_SOLVERS && info.n_ac_buses <= 500  # Skip for large cases
            # NLsolve benchmark
            print("NLsolve (warmup=$NUM_WARMUP_RUNS, runs=$NUM_TIMED_RUNS)... ")
            nlsolve_result = benchmark_nlsolve(sys)
            if !isnothing(nlsolve_result)
                println("done")
                println("  Time: $(@sprintf("%.3f", nlsolve_result.mean_time)) ± $(@sprintf("%.3f", nlsolve_result.std_time)) ms")
                println("  Converged: $(nlsolve_result.converged), Iterations: $(nlsolve_result.iterations)")
            else
                println("skipped (not available)")
            end
            
            # JuMP+Ipopt benchmark (even slower, skip for >100 buses)
            if info.n_ac_buses <= 100
                print("JuMP+Ipopt (warmup=$NUM_WARMUP_RUNS, runs=$NUM_TIMED_RUNS)... ")
                jump_result = benchmark_jump_ipopt(sys)
                if !isnothing(jump_result)
                    println("done")
                    println("  Time: $(@sprintf("%.3f", jump_result.mean_time)) ± $(@sprintf("%.3f", jump_result.std_time)) ms")
                    println("  Converged: $(jump_result.converged), Iterations: $(jump_result.iterations)")
                else
                    println("skipped (not available)")
                end
            end
        elseif info.n_ac_buses > 500
            println("Reference solvers skipped (system too large)")
        end
        
        # Validation
        if !isnothing(nlsolve_result)
            val = validate_results(nr_result, nlsolve_result, "NLsolve")
            if !isnothing(val)
                println()
                println("Validation vs NLsolve:")
                println("  Vm max error: $(@sprintf("%.2e", val.vm_max_error)) p.u.")
                println("  Va max error: $(@sprintf("%.2e", val.va_max_error)) rad")
            end
        end
        
        if !isnothing(jump_result)
            val = validate_results(nr_result, jump_result, "JuMP+Ipopt")
            if !isnothing(val)
                println()
                println("Validation vs JuMP+Ipopt:")
                println("  Vm max error: $(@sprintf("%.2e", val.vm_max_error)) p.u.")
                println("  Va max error: $(@sprintf("%.2e", val.va_max_error)) rad")
            end
        end
        
        # Store results
        results[name] = (
            info = info,
            nr = nr_result,
            nlsolve = nlsolve_result,
            jump = jump_result
        )
        
        println()
    end
    
    # Summary table
    println("=" ^ 80)
    println("SUMMARY")
    println("=" ^ 80)
    println()
    println(@sprintf("%-12s %8s %12s %12s %10s %8s", 
                     "Case", "Buses", "Branches", "Unknowns", "NR Time", "Speedup"))
    println("-" ^ 70)
    
    for (name, _) in TEST_CASES
        r = results[name]
        nr_time = r.nr.mean_time
        
        # Calculate speedup vs NLsolve
        speedup = if !isnothing(r.nlsolve)
            @sprintf("%.1fx", r.nlsolve.mean_time / nr_time)
        else
            "-"
        end
        
        println(@sprintf("%-12s %8d %12d %12d %8.2f ms %8s",
                        name, r.info.n_ac_buses, r.info.n_ac_branches, 
                        r.info.n_unknowns, nr_time, speedup))
    end
    
    println()
    println("All benchmarks completed!")
    
    return results
end

# =====================================================================
# Scaling Analysis
# =====================================================================
function analyze_scaling(results)
    println()
    println("=" ^ 80)
    println("SCALING ANALYSIS")
    println("=" ^ 80)
    println()
    
    sizes = Float64[]
    times = Float64[]
    
    for (name, _) in TEST_CASES
        r = results[name]
        push!(sizes, r.info.n_ac_buses)
        push!(times, r.nr.mean_time)
    end
    
    # Fit power law: time = a * size^b
    log_sizes = log.(sizes)
    log_times = log.(times)
    
    # Linear regression on log-log scale
    n = length(sizes)
    sum_x = sum(log_sizes)
    sum_y = sum(log_times)
    sum_xy = sum(log_sizes .* log_times)
    sum_x2 = sum(log_sizes.^2)
    
    b = (n * sum_xy - sum_x * sum_y) / (n * sum_x2 - sum_x^2)
    a = exp((sum_y - b * sum_x) / n)
    
    println("Power law fit: time = $(@sprintf("%.4f", a)) × n^$(@sprintf("%.2f", b))")
    println()
    println("Scaling interpretation:")
    if b < 1.5
        println("  - Sub-quadratic scaling: excellent!")
    elseif b < 2.0
        println("  - Near-quadratic scaling: good for power flow")
    elseif b < 2.5
        println("  - Super-quadratic but manageable")
    else
        println("  - High polynomial complexity")
    end
    
    println()
    println("Predicted solve times:")
    for n in [5000, 10000, 20000]
        t_pred = a * n^b
        println("  n=$n buses: $(@sprintf("%.1f", t_pred)) ms")
    end
end

# =====================================================================
# Run Benchmark
# =====================================================================
if abspath(PROGRAM_FILE) == @__FILE__
    results = run_benchmarks()
    analyze_scaling(results)
end
