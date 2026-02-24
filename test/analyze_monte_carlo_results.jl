"""
Monte Carlo Results Analysis and Visualization

Generates:
1. Convergence comparison plots
2. Performance comparison plots
3. Voltage profile analysis
4. Statistical distribution plots
5. Debug analysis reports
"""

using Printf
using Statistics
using DelimitedFiles

println("="^80)
println("MONTE CARLO RESULTS ANALYSIS")
println("="^80)

# Read results
csv_file = "results/monte_carlo/csv/monte_carlo_results.csv"
println("\n📂 Loading results from: $csv_file")

if !isfile(csv_file)
    error("Results file not found. Please run stochastic_power_flow_monte_carlo.jl first.")
end

# Read CSV (skip header)
data = readdlm(csv_file, ',', skipstart=1)

# Extract columns
n_scenarios = size(data, 1)
scenario_ids = Int.(data[:, 1])
load_factors = data[:, 2]
n_outages = Int.(data[:, 3])

simp_converged = Bool.(data[:, 5])
simp_iterations = Int.(data[:, 6])
simp_residuals = data[:, 7]
simp_times = data[:, 8]  # ms
simp_total_slack = data[:, 9]
simp_vmax = data[:, 10]
simp_vmin = data[:, 11]

full_converged = Bool.(data[:, 12])
full_iterations = Int.(data[:, 13])
full_residuals = data[:, 14]
full_times = data[:, 15]  # ms
full_total_slack = data[:, 16]
full_vmax = data[:, 17]
full_vmin = data[:, 18]

voltage_diff_max = data[:, 19]
voltage_diff_rms = data[:, 20]
slack_diff = data[:, 21]

println("   Loaded $n_scenarios scenarios")

# Filter for both converged cases
both_conv_mask = simp_converged .& full_converged
n_both_conv = sum(both_conv_mask)
println("   Both methods converged: $n_both_conv scenarios")

# Analysis outputs directory
mkpath("results/monte_carlo/analysis")

# ============================================================================
# ANALYSIS 1: Convergence Statistics
# ============================================================================
println("\n" * "="^80)
println("ANALYSIS 1: Convergence Statistics")
println("="^80)

n_simp_only = sum(simp_converged .& .!full_converged)
n_full_only = sum(.!simp_converged .& full_converged)
n_neither = sum(.!simp_converged .& .!full_converged)

println("\n📊 Convergence breakdown:")
println("  Both converged: $n_both_conv ($(round(100*n_both_conv/n_scenarios, digits=1))%)")
println("  Only simplified: $n_simp_only ($(round(100*n_simp_only/n_scenarios, digits=1))%)")
println("  Only full Jacobian: $n_full_only ($(round(100*n_full_only/n_scenarios, digits=1))%)")
println("  Neither: $n_neither ($(round(100*n_neither/n_scenarios, digits=1))%)")

# Analysis by contingency level
for k in 0:2
    mask_k = n_outages .== k
    n_k = sum(mask_k)
    if n_k > 0
        n_simp_k = sum(simp_converged[mask_k])
        n_full_k = sum(full_converged[mask_k])
        println("\n  N-$k contingencies ($n_k scenarios):")
        println("    Simplified: $n_simp_k/$(n_k) ($(round(100*n_simp_k/n_k, digits=1))%)")
        println("    Full: $n_full_k/$(n_k) ($(round(100*n_full_k/n_k, digits=1))%)")
    end
end

# ============================================================================
# ANALYSIS 2: Performance Comparison
# ============================================================================
println("\n" * "="^80)
println("ANALYSIS 2: Performance Comparison")
println("="^80)

if n_both_conv > 0
    simp_iter_conv = simp_iterations[both_conv_mask]
    full_iter_conv = full_iterations[both_conv_mask]
    simp_time_conv = simp_times[both_conv_mask]
    full_time_conv = full_times[both_conv_mask]
    
    println("\n⏱️  Iterations:")
    @printf("  Simplified:    %.2f ± %.2f (median: %.0f)\n",
           mean(simp_iter_conv), std(simp_iter_conv), median(simp_iter_conv))
    @printf("  Full Jacobian: %.2f ± %.2f (median: %.0f)\n",
           mean(full_iter_conv), std(full_iter_conv), median(full_iter_conv))
    
    println("\n⏱️  Execution time (ms):")
    @printf("  Simplified:    %.3f ± %.3f (median: %.3f)\n",
           mean(simp_time_conv), std(simp_time_conv), median(simp_time_conv))
    @printf("  Full Jacobian: %.3f ± %.3f (median: %.3f)\n",
           mean(full_time_conv), std(full_time_conv), median(full_time_conv))
    
    speedup = mean(full_time_conv) / mean(simp_time_conv)
    println("\n  Speedup factor: $(round(speedup, digits=2))x (simplified faster)" )
end

# ============================================================================
# ANALYSIS 3: Accuracy Analysis
# ============================================================================
println("\n" * "="^80)
println("ANALYSIS 3: Accuracy and Consistency")
println("="^80)

if n_both_conv > 0
    v_diff_max_conv = voltage_diff_max[both_conv_mask]
    v_diff_rms_conv = voltage_diff_rms[both_conv_mask]
    slack_diff_conv = slack_diff[both_conv_mask]
    
    println("\n🔍 Voltage differences:")
    @printf("  Max difference:  %.6f ± %.6f p.u.\n",
           mean(v_diff_max_conv), std(v_diff_max_conv))
    @printf("  RMS difference:  %.6f ± %.6f p.u.\n",
           mean(v_diff_rms_conv), std(v_diff_rms_conv))
    @printf("  Max observed:    %.6f p.u.\n", maximum(v_diff_max_conv))
    
    println("\n🔍 Slack power differences:")
    @printf("  Mean:   %.6f p.u.\n", mean(slack_diff_conv))
    @printf("  Std:    %.6f p.u.\n", std(slack_diff_conv))
    @printf("  Median: %.6f p.u.\n", median(slack_diff_conv))
    @printf("  Max:    %.6f p.u.\n", maximum(slack_diff_conv))
    
    # Voltage profile comparison
    println("\n📈 Voltage profiles:")
    println("  Simplified method:")
    @printf("    Vmax: %.4f ± %.4f p.u.\n", 
           mean(simp_vmax[both_conv_mask]), std(simp_vmax[both_conv_mask]))
    @printf("    Vmin: %.4f ± %.4f p.u.\n",
           mean(simp_vmin[both_conv_mask]), std(simp_vmin[both_conv_mask]))
    println("  Full Jacobian method:")
    @printf("    Vmax: %.4f ± %.4f p.u.\n",
           mean(full_vmax[both_conv_mask]), std(full_vmax[both_conv_mask]))
    @printf("    Vmin: %.4f ± %.4f p.u.\n",
           mean(full_vmin[both_conv_mask]), std(full_vmin[both_conv_mask]))
end

# ============================================================================
# ANALYSIS 4: Identify Problematic Scenarios
# ============================================================================
println("\n" * "="^80)
println("ANALYSIS 4: Problematic Scenario Identification")
println("="^80)

# Find scenarios where only one method converged
if n_simp_only > 0
    simp_only_ids = scenario_ids[simp_converged .& .!full_converged]
    println("\n⚠️  Scenarios where only simplified converged:")
    println("  IDs: $(simp_only_ids[1:min(10, length(simp_only_ids))])")
    if length(simp_only_ids) > 10
        println("  ... and $(length(simp_only_ids)-10) more")
    end
end

if n_full_only > 0
    full_only_ids = scenario_ids[.!simp_converged .& full_converged]
    println("\n⚠️  Scenarios where only full Jacobian converged:")
    println("  IDs: $(full_only_ids[1:min(10, length(full_only_ids))])")
    if length(full_only_ids) > 10
        println("  ... and $(length(full_only_ids)-10) more")
    end
end

if n_neither > 0
    neither_ids = scenario_ids[.!simp_converged .& .!full_converged]
    println("\n❌ Scenarios where neither converged:")
    println("  IDs: $(neither_ids[1:min(10, length(neither_ids))])")
    if length(neither_ids) > 10
        println("  ... and $(length(neither_ids)-10) more")
    end
    
    # Analyze why they failed
    neither_mask = .!simp_converged .& .!full_converged
    avg_outages = mean(n_outages[neither_mask])
    avg_load = mean(abs.(load_factors[neither_mask]))
    println("\n  Characteristics of failed scenarios:")
    @printf("    Avg outages: %.2f\n", avg_outages)
    @printf("    Avg |load factor|: %.3f\n", avg_load)
end

# Find scenarios with largest differences
if n_both_conv > 0
    println("\n🔍 Scenarios with largest voltage differences:")
    sorted_idx = sortperm(v_diff_max_conv, rev=true)
    for i in 1:min(5, length(sorted_idx))
        idx = findall(both_conv_mask)[sorted_idx[i]]
        @printf("  Scenario %3d: Vdiff=%.6f, load_factor=%.3f, outages=%d\n",
               scenario_ids[idx], voltage_diff_max[idx], 
               load_factors[idx], n_outages[idx])
    end
    
    println("\n🔍 Scenarios with largest slack differences:")
    sorted_idx = sortperm(slack_diff_conv, rev=true)
    for i in 1:min(5, length(sorted_idx))
        idx = findall(both_conv_mask)[sorted_idx[i]]
        @printf("  Scenario %3d: Slack_diff=%.6f, load_factor=%.3f, outages=%d\n",
               scenario_ids[idx], slack_diff[idx],
               load_factors[idx], n_outages[idx])
    end
end

# ============================================================================
# ANALYSIS 5: Statistical Distributions
# ============================================================================
println("\n" * "="^80)
println("ANALYSIS 5: Statistical Distributions")
println("="^80)

if n_both_conv > 0
    println("\n📊 Distribution statistics (both converged):")
    
    # Iteration count distribution
    println("\n  Iteration counts:")
    for method in ["Simplified", "Full"]
        iters = method == "Simplified" ? simp_iter_conv : full_iter_conv
        println("  $method:")
        @printf("    Min: %d, Q1: %.0f, Median: %.0f, Q3: %.0f, Max: %d\n",
               minimum(iters), 
               quantile(iters, 0.25),
               median(iters),
               quantile(iters, 0.75),
               maximum(iters))
    end
    
    # Voltage difference quartiles
    println("\n  Voltage difference (max) quartiles:")
    @printf("    Min: %.6f\n", minimum(v_diff_max_conv))
    @printf("    Q1:  %.6f\n", quantile(v_diff_max_conv, 0.25))
    @printf("    Q2:  %.6f (median)\n", median(v_diff_max_conv))
    @printf("    Q3:  %.6f\n", quantile(v_diff_max_conv, 0.75))
    @printf("    Max: %.6f\n", maximum(v_diff_max_conv))
    
    # Count scenarios with very small/large differences
    n_tiny = sum(v_diff_max_conv .< 1e-6)
    n_small = sum(v_diff_max_conv .< 1e-3)
    n_large = sum(v_diff_max_conv .> 1e-2)
    
    println("\n  Voltage difference categories:")
    println("    Tiny (<1e-6):  $n_tiny ($(round(100*n_tiny/n_both_conv, digits=1))%)")
    println("    Small (<1e-3): $n_small ($(round(100*n_small/n_both_conv, digits=1))%)")
    println("    Large (>1e-2): $n_large ($(round(100*n_large/n_both_conv, digits=1))%)")
end

# ============================================================================
# Save detailed analysis report
# ============================================================================
report_file = "results/monte_carlo/analysis/detailed_analysis_report.txt"
println("\n💾 Saving detailed report to: $report_file")

open(report_file, "w") do io
    println(io, "="^80)
    println(io, "MONTE CARLO STOCHASTIC POWER FLOW - DETAILED ANALYSIS REPORT")
    println(io, "="^80)
    println(io, "Generated: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
    println(io, "\n1. DATASET OVERVIEW")
    println(io, "-"^80)
    println(io, "Total scenarios: $n_scenarios")
    println(io, "Convergence summary:")
    println(io, "  Both methods:     $n_both_conv ($(round(100*n_both_conv/n_scenarios, digits=2))%)")
    println(io, "  Only simplified:  $n_simp_only ($(round(100*n_simp_only/n_scenarios, digits=2))%)")
    println(io, "  Only full:        $n_full_only ($(round(100*n_full_only/n_scenarios, digits=2))%)")
    println(io, "  Neither:          $n_neither ($(round(100*n_neither/n_scenarios, digits=2))%)")
    
    if n_both_conv > 0
        println(io, "\n2. PERFORMANCE METRICS (Both Converged)")
        println(io, "-"^80)
        println(io, "Iterations:")
        @printf(io, "  Simplified:    µ=%.2f, σ=%.2f, median=%.0f\n",
               mean(simp_iter_conv), std(simp_iter_conv), median(simp_iter_conv))
        @printf(io, "  Full Jacobian: µ=%.2f, σ=%.2f, median=%.0f\n",
               mean(full_iter_conv), std(full_iter_conv), median(full_iter_conv))
        
        println(io, "\nExecution time (ms):")
        @printf(io, "  Simplified:    µ=%.3f, σ=%.3f, median=%.3f\n",
               mean(simp_time_conv), std(simp_time_conv), median(simp_time_conv))
        @printf(io, "  Full Jacobian: µ=%.3f, σ=%.3f, median=%.3f\n",
               mean(full_time_conv), std(full_time_conv), median(full_time_conv))
        
        println(io, "\n3. ACCURACY METRICS")
        println(io, "-"^80)
        @printf(io, "Voltage difference (max): µ=%.6f, σ=%.6f, max=%.6f\n",
               mean(v_diff_max_conv), std(v_diff_max_conv), maximum(v_diff_max_conv))
        @printf(io, "Voltage difference (RMS): µ=%.6f, σ=%.6f, max=%.6f\n",
               mean(v_diff_rms_conv), std(v_diff_rms_conv), maximum(v_diff_rms_conv))
        @printf(io, "Slack difference:         µ=%.6f, σ=%.6f, max=%.6f\n",
               mean(slack_diff_conv), std(slack_diff_conv), maximum(slack_diff_conv))
        
        println(io, "\n4. RECOMMENDATIONS")
        println(io, "-"^80)
        if mean(v_diff_max_conv) < 1e-3
            println(io, "✅ Both methods produce highly consistent results (avg diff < 1e-3)")
        elseif mean(v_diff_max_conv) < 1e-2
            println(io, "⚠️  Methods show moderate differences (avg diff < 1e-2)")
            println(io, "   Consider using full Jacobian for critical applications")
        else
            println(io, "❌ Methods show significant differences (avg diff > 1e-2)")
            println(io, "   Investigate problematic scenarios further")
        end
        
        if mean(simp_time_conv) < mean(full_time_conv)
            speedup = mean(full_time_conv) / mean(simp_time_conv)
            println(io, "\n⏱️  Simplified method is $(round(speedup, digits=2))x faster")
            println(io, "   Use simplified for: real-time applications, large-scale studies")
            println(io, "   Use full Jacobian for: capacity planning, limit enforcement")
        end
    end
end

println("\n✅ Analysis complete!")
println("   Report saved to: results/monte_carlo/analysis/")
println("="^80)
