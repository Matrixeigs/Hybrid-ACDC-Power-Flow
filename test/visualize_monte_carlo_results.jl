"""
Monte Carlo Results Visualization

Generates publication-quality figures for stochastic power flow analysis.

Requires: Plots, StatsPlots (install with: using Pkg; Pkg.add(["Plots", "StatsPlots"]))
"""

using Printf
using Statistics
using DelimitedFiles

# Try to load plotting libraries
try
    using Plots
    using StatsPlots
    gr()  # Use GR backend for fast rendering
catch e
    println("⚠️  Plotting libraries not found!")
    println("   Please install with: using Pkg; Pkg.add([\"Plots\", \"StatsPlots\"])")
    println("   Then re-run this script.")
    exit(1)
end

println("="^80)
println("MONTE CARLO RESULTS VISUALIZATION")
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
simp_times = data[:, 8]
simp_total_slack = data[:, 9]
simp_vmax = data[:, 10]
simp_vmin = data[:, 11]

full_converged = Bool.(data[:, 12])
full_iterations = Int.(data[:, 13])
full_residuals = data[:, 14]
full_times = data[:, 15]
full_total_slack = data[:, 16]
full_vmax = data[:, 17]
full_vmin = data[:, 18]

voltage_diff_max = data[:, 19]
voltage_diff_rms = data[:, 20]
slack_diff = data[:, 21]

println("   Loaded $n_scenarios scenarios")

# Filter for both converged
both_conv_mask = simp_converged .& full_converged
n_both_conv = sum(both_conv_mask)
println("   Both methods converged: $n_both_conv scenarios")

# Create output directory
fig_dir = "results/monte_carlo/figures"
mkpath(fig_dir)
println("   Saving figures to: $fig_dir")

# Set default plot parameters for publication quality
default(
    fontfamily="Computer Modern",
    titlefontsize=12,
    guidefontsize=11,
    tickfontsize=9,
    legendfontsize=9,
    dpi=300,
    size=(800, 600),
    margin=5Plots.mm
)

# ============================================================================
# FIGURE 1: Convergence Comparison by Contingency Level
# ============================================================================
println("\n📊 Generating Figure 1: Convergence rates...")

# Calculate convergence rates by contingency level
conv_data = zeros(3, 2)  # [N-0, N-1, N-2] × [Simplified, Full]
contingency_labels = ["N-0", "N-1", "N-2"]

for k in 0:2
    mask_k = n_outages .== k
    n_k = sum(mask_k)
    if n_k > 0
        conv_data[k+1, 1] = 100 * sum(simp_converged[mask_k]) / n_k
        conv_data[k+1, 2] = 100 * sum(full_converged[mask_k]) / n_k
    end
end

p1 = groupedbar(
    contingency_labels,
    conv_data,
    bar_position=:dodge,
    bar_width=0.7,
    label=["Simplified" "Full Jacobian"],
    xlabel="Contingency Level",
    ylabel="Convergence Rate (%)",
    title="Convergence Comparison by Contingency Level",
    ylims=(0, 105),
    legend=:bottomleft,
    color=[:steelblue :coral],
    grid=true,
    gridstyle=:dot,
    gridalpha=0.3
)

# Add value labels on bars
for i in 1:3
    for j in 1:2
        val = conv_data[i, j]
        if val > 0
            annotate!(p1, i + (j-1.5)*0.2, val + 2, text(@sprintf("%.1f", val), 8))
        end
    end
end

savefig(p1, joinpath(fig_dir, "fig1_convergence_comparison.png"))
println("   ✅ Saved: fig1_convergence_comparison.png")

# ============================================================================
# FIGURE 2: Iteration Count Comparison
# ============================================================================
println("\n📊 Generating Figure 2: Iteration counts...")

if n_both_conv > 0
    simp_iter_conv = simp_iterations[both_conv_mask]
    full_iter_conv = full_iterations[both_conv_mask]
    
    p2 = plot(layout=(1, 2), size=(1200, 500))
    
    # Histogram comparison
    histogram!(p2, subplot=1, simp_iter_conv, 
              bins=0:1:maximum([simp_iter_conv; full_iter_conv])+1,
              alpha=0.7, label="Simplified", color=:steelblue,
              xlabel="Iterations", ylabel="Frequency",
              title="Iteration Count Distribution")
    histogram!(p2, subplot=1, full_iter_conv,
              bins=0:1:maximum([simp_iter_conv; full_iter_conv])+1,
              alpha=0.7, label="Full Jacobian", color=:coral)
    
    # Box plot comparison
    boxplot!(p2, subplot=2, ["Simplified"], simp_iter_conv,
            fillalpha=0.7, color=:steelblue, label="")
    boxplot!(p2, subplot=2, ["Full Jacobian"], full_iter_conv,
            fillalpha=0.7, color=:coral, label="")
    plot!(p2, subplot=2, ylabel="Iterations", title="Iteration Count Statistics",
         xtickfontsize=10)
    
    savefig(p2, joinpath(fig_dir, "fig2_iteration_comparison.png"))
    println("   ✅ Saved: fig2_iteration_comparison.png")
end

# ============================================================================
# FIGURE 3: Execution Time Comparison
# ============================================================================
println("\n📊 Generating Figure 3: Execution times...")

if n_both_conv > 0
    simp_time_conv = simp_times[both_conv_mask]
    full_time_conv = full_times[both_conv_mask]
    
    p3 = plot(layout=(1, 2), size=(1200, 500))
    
    # Histogram
    histogram!(p3, subplot=1, simp_time_conv,
              bins=20, alpha=0.7, label="Simplified", color=:steelblue,
              xlabel="Execution Time (ms)", ylabel="Frequency",
              title="Execution Time Distribution")
    histogram!(p3, subplot=1, full_time_conv,
              bins=20, alpha=0.7, label="Full Jacobian", color=:coral)
    
    # Box plot
    boxplot!(p3, subplot=2, ["Simplified"], simp_time_conv,
            fillalpha=0.7, color=:steelblue, label="")
    boxplot!(p3, subplot=2, ["Full Jacobian"], full_time_conv,
            fillalpha=0.7, color=:coral, label="")
    plot!(p3, subplot=2, ylabel="Execution Time (ms)",
         title="Execution Time Statistics", xtickfontsize=10)
    
    savefig(p3, joinpath(fig_dir, "fig3_time_comparison.png"))
    println("   ✅ Saved: fig3_time_comparison.png")
end

# ============================================================================
# FIGURE 4: Voltage Difference Analysis
# ============================================================================
println("\n📊 Generating Figure 4: Voltage differences...")

if n_both_conv > 0
    v_diff_max_conv = voltage_diff_max[both_conv_mask]
    v_diff_rms_conv = voltage_diff_rms[both_conv_mask]
    load_factors_conv = load_factors[both_conv_mask]
    n_outages_conv = n_outages[both_conv_mask]
    
    p4 = plot(layout=(2, 2), size=(1200, 1000))
    
    # Max voltage difference histogram
    histogram!(p4, subplot=1, v_diff_max_conv,
              bins=30, alpha=0.7, color=:purple,
              xlabel="Max Voltage Difference (p.u.)", ylabel="Frequency",
              title="Max Voltage Difference Distribution", label="")
    vline!(p4, subplot=1, [mean(v_diff_max_conv)],
          color=:red, linewidth=2, linestyle=:dash,
          label=@sprintf("Mean: %.6f", mean(v_diff_max_conv)))
    
    # RMS voltage difference histogram
    histogram!(p4, subplot=2, v_diff_rms_conv,
              bins=30, alpha=0.7, color=:teal,
              xlabel="RMS Voltage Difference (p.u.)", ylabel="Frequency",
              title="RMS Voltage Difference Distribution", label="")
    vline!(p4, subplot=2, [mean(v_diff_rms_conv)],
          color=:red, linewidth=2, linestyle=:dash,
          label=@sprintf("Mean: %.6f", mean(v_diff_rms_conv)))
    
    # Voltage difference vs load factor
    scatter!(p4, subplot=3, load_factors_conv, v_diff_max_conv,
            alpha=0.6, color=:steelblue, markersize=4,
            xlabel="Load Factor", ylabel="Max Voltage Diff (p.u.)",
            title="Voltage Difference vs Load Factor", label="")
    
    # Voltage difference vs number of outages
    scatter!(p4, subplot=4, n_outages_conv, v_diff_max_conv,
            alpha=0.6, color=:coral, markersize=4,
            xlabel="Number of Outages", ylabel="Max Voltage Diff (p.u.)",
            title="Voltage Difference vs Contingency", label="")
    
    savefig(p4, joinpath(fig_dir, "fig4_voltage_difference_analysis.png"))
    println("   ✅ Saved: fig4_voltage_difference_analysis.png")
end

# ============================================================================
# FIGURE 5: Slack Power Distribution Comparison
# ============================================================================
println("\n📊 Generating Figure 5: Slack power distribution...")

if n_both_conv > 0
    simp_slack_conv = simp_total_slack[both_conv_mask]
    full_slack_conv = full_total_slack[both_conv_mask]
    slack_diff_conv = slack_diff[both_conv_mask]
    
    p5 = plot(layout=(1, 3), size=(1500, 500))
    
    # Slack power scatter
    scatter!(p5, subplot=1, simp_slack_conv, full_slack_conv,
            alpha=0.6, color=:steelblue, markersize=4,
            xlabel="Simplified Slack (p.u.)", ylabel="Full Jacobian Slack (p.u.)",
            title="Slack Power Correlation", label="")
    # Add y=x line
    slack_range = [minimum([simp_slack_conv; full_slack_conv]),
                   maximum([simp_slack_conv; full_slack_conv])]
    plot!(p5, subplot=1, slack_range, slack_range,
         color=:red, linestyle=:dash, linewidth=2, label="y=x")
    
    # Slack difference histogram
    histogram!(p5, subplot=2, slack_diff_conv,
              bins=30, alpha=0.7, color=:orange,
              xlabel="Slack Difference (p.u.)", ylabel="Frequency",
              title="Slack Power Difference Distribution", label="")
    vline!(p5, subplot=2, [mean(slack_diff_conv)],
          color=:red, linewidth=2, linestyle=:dash,
          label=@sprintf("Mean: %.6f", mean(slack_diff_conv)))
    
    # Slack difference vs load factor
    scatter!(p5, subplot=3, load_factors_conv, slack_diff_conv,
            alpha=0.6, color=:purple, markersize=4,
            xlabel="Load Factor", ylabel="Slack Difference (p.u.)",
            title="Slack Diff vs Load Factor", label="")
    
    savefig(p5, joinpath(fig_dir, "fig5_slack_comparison.png"))
    println("   ✅ Saved: fig5_slack_comparison.png")
end

# ============================================================================
# FIGURE 6: Voltage Profile Comparison
# ============================================================================
println("\n📊 Generating Figure 6: Voltage profiles...")

if n_both_conv > 0
    simp_vmax_conv = simp_vmax[both_conv_mask]
    simp_vmin_conv = simp_vmin[both_conv_mask]
    full_vmax_conv = full_vmax[both_conv_mask]
    full_vmin_conv = full_vmin[both_conv_mask]
    
    p6 = plot(layout=(2, 2), size=(1200, 1000))
    
    # Vmax comparison
    scatter!(p6, subplot=1, simp_vmax_conv, full_vmax_conv,
            alpha=0.6, color=:red, markersize=4,
            xlabel="Simplified Vmax (p.u.)", ylabel="Full Jacobian Vmax (p.u.)",
            title="Maximum Voltage Comparison", label="")
    vmax_range = [minimum([simp_vmax_conv; full_vmax_conv]),
                  maximum([simp_vmax_conv; full_vmax_conv])]
    plot!(p6, subplot=1, vmax_range, vmax_range,
         color=:black, linestyle=:dash, linewidth=2, label="y=x")
    
    # Vmin comparison
    scatter!(p6, subplot=2, simp_vmin_conv, full_vmin_conv,
            alpha=0.6, color=:blue, markersize=4,
            xlabel="Simplified Vmin (p.u.)", ylabel="Full Jacobian Vmin (p.u.)",
            title="Minimum Voltage Comparison", label="")
    vmin_range = [minimum([simp_vmin_conv; full_vmin_conv]),
                  maximum([simp_vmin_conv; full_vmin_conv])]
    plot!(p6, subplot=2, vmin_range, vmin_range,
         color=:black, linestyle=:dash, linewidth=2, label="y=x")
    
    # Voltage spread (Vmax - Vmin)
    simp_spread = simp_vmax_conv .- simp_vmin_conv
    full_spread = full_vmax_conv .- full_vmin_conv
    
    boxplot!(p6, subplot=3, ["Simplified"], simp_spread,
            fillalpha=0.7, color=:steelblue, label="")
    boxplot!(p6, subplot=3, ["Full Jacobian"], full_spread,
            fillalpha=0.7, color=:coral, label="")
    plot!(p6, subplot=3, ylabel="Voltage Spread (p.u.)",
         title="Voltage Spread Statistics", xtickfontsize=10)
    
    # Voltage spread vs contingency
    scatter!(p6, subplot=4, n_outages_conv, simp_spread,
            alpha=0.6, color=:steelblue, markersize=4,
            label="Simplified")
    scatter!(p6, subplot=4, n_outages_conv .+ 0.1, full_spread,
            alpha=0.6, color=:coral, markersize=4,
            xlabel="Number of Outages", ylabel="Voltage Spread (p.u.)",
            title="Voltage Spread vs Contingency", label="Full Jacobian")
    
    savefig(p6, joinpath(fig_dir, "fig6_voltage_profile_comparison.png"))
    println("   ✅ Saved: fig6_voltage_profile_comparison.png")
end

# ============================================================================
# FIGURE 7: Performance vs Accuracy Tradeoff
# ============================================================================
println("\n📊 Generating Figure 7: Performance vs accuracy...")

if n_both_conv > 0
    p7 = plot(layout=(1, 2), size=(1200, 500))
    
    # Iterations vs voltage accuracy
    scatter!(p7, subplot=1, simp_iter_conv, v_diff_max_conv,
            alpha=0.6, color=:steelblue, markersize=5,
            xlabel="Simplified Method Iterations",
            ylabel="Max Voltage Difference (p.u.)",
            title="Simplified Iterations vs Accuracy", label="")
    
    scatter!(p7, subplot=2, full_iter_conv, v_diff_max_conv,
            alpha=0.6, color=:coral, markersize=5,
            xlabel="Full Jacobian Method Iterations",
            ylabel="Max Voltage Difference (p.u.)",
            title="Full Jacobian Iterations vs Accuracy", label="")
    
    savefig(p7, joinpath(fig_dir, "fig7_performance_accuracy_tradeoff.png"))
    println("   ✅ Saved: fig7_performance_accuracy_tradeoff.png")
end

# ============================================================================
# FIGURE 8: Comprehensive Summary Dashboard
# ============================================================================
println("\n📊 Generating Figure 8: Summary dashboard...")

if n_both_conv > 0
    p8 = plot(layout=grid(2, 3), size=(1800, 1200))
    
    # 1. Convergence rates
    conv_rates = [100*sum(simp_converged)/n_scenarios, 
                  100*sum(full_converged)/n_scenarios]
    bar!(p8, subplot=1, ["Simplified", "Full"],
        conv_rates, color=[:steelblue :coral],
        ylabel="Convergence Rate (%)", title="Overall Convergence",
        ylims=(0, 105), legend=false)
    for i in 1:2
        annotate!(p8, subplot=1, i, conv_rates[i] + 2,
                 text(@sprintf("%.1f%%", conv_rates[i]), 10))
    end
    
    # 2. Iteration comparison
    boxplot!(p8, subplot=2, ["Simp"], simp_iter_conv,
            fillalpha=0.7, color=:steelblue, label="")
    boxplot!(p8, subplot=2, ["Full"], full_iter_conv,
            fillalpha=0.7, color=:coral, label="")
    plot!(p8, subplot=2, ylabel="Iterations",
         title="Iteration Count", xtickfontsize=9)
    
    # 3. Time comparison
    boxplot!(p8, subplot=3, ["Simp"], simp_time_conv,
            fillalpha=0.7, color=:steelblue, label="")
    boxplot!(p8, subplot=3, ["Full"], full_time_conv,
            fillalpha=0.7, color=:coral, label="")
    plot!(p8, subplot=3, ylabel="Time (ms)",
         title="Execution Time", xtickfontsize=9)
    
    # 4. Voltage difference
    histogram!(p8, subplot=4, v_diff_max_conv,
              bins=25, alpha=0.7, color=:purple,
              xlabel="Max V Diff (p.u.)", ylabel="Count",
              title="Voltage Difference", legend=false)
    
    # 5. Slack difference
    histogram!(p8, subplot=5, slack_diff_conv,
              bins=25, alpha=0.7, color=:orange,
              xlabel="Slack Diff (p.u.)", ylabel="Count",
              title="Slack Power Difference", legend=false)
    
    # 6. Summary statistics table (as plot)
    stats_text = """
    SUMMARY STATISTICS
    
    Scenarios: $n_scenarios
    Both converged: $n_both_conv
    
    Iterations (mean±std):
      Simp: $(round(mean(simp_iter_conv), digits=1))±$(round(std(simp_iter_conv), digits=1))
      Full: $(round(mean(full_iter_conv), digits=1))±$(round(std(full_iter_conv), digits=1))
    
    Time (ms):
      Simp: $(round(mean(simp_time_conv), digits=2))
      Full: $(round(mean(full_time_conv), digits=2))
    
    Voltage diff: $(round(mean(v_diff_max_conv), digits=6))
    Slack diff: $(round(mean(slack_diff_conv), digits=6))
    """
    
    plot!(p8, subplot=6, [], [], axis=false, grid=false,
         annotations=(0.1, 0.5, text(stats_text, :left, 9, :courier)))
    
    savefig(p8, joinpath(fig_dir, "fig8_summary_dashboard.png"))
    println("   ✅ Saved: fig8_summary_dashboard.png")
end

# ============================================================================
# Summary
# ============================================================================
println("\n" * "="^80)
println("✅ VISUALIZATION COMPLETE!")
println("="^80)
println("\nGenerated figures:")
println("  1. fig1_convergence_comparison.png - Convergence by contingency level")
println("  2. fig2_iteration_comparison.png - Iteration count analysis")
println("  3. fig3_time_comparison.png - Execution time analysis")
println("  4. fig4_voltage_difference_analysis.png - Voltage accuracy metrics")
println("  5. fig5_slack_comparison.png - Slack power distribution")
println("  6. fig6_voltage_profile_comparison.png - Voltage profiles")
println("  7. fig7_performance_accuracy_tradeoff.png - Performance vs accuracy")
println("  8. fig8_summary_dashboard.png - Comprehensive summary")
println("\nAll figures saved to: $fig_dir")
println("="^80)
