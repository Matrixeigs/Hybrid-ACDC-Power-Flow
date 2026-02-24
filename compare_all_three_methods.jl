#!/usr/bin/env julia
#
# Three-way comparison: NLsolve vs Simplified NR vs Full Jacobian
# Compares actual solution voltages from all three methods
#

cd(@__DIR__)
using Pkg
Pkg.activate(".")

println("Loading HybridACDCPowerFlow module...")
include("src/HybridACDCPowerFlow.jl")
using .HybridACDCPowerFlow
using .HybridACDCPowerFlow: PQ_MODE, VDC_Q, VDC_VAC, ACBus, ACBranch, VSCConverter
using Random
using Printf
using Statistics
using LinearAlgebra
using CSV
using DataFrames

println("\n" * "="^80)
println("  THREE-WAY ALGORITHM COMPARISON")
println("  NLsolve vs Simplified NR vs Full Jacobian")
println("="^80)
println()

# Load Monte Carlo results to identify converged scenarios
results_file = "results/monte_carlo/csv/monte_carlo_results.csv"
if !isfile(results_file)
    println("❌ Results file not found. Please run Monte Carlo simulation first.")
    exit(1)
end

df = CSV.read(results_file, DataFrame)

# Filter for scenarios where ALL methods converged
nlsolve_feasible = Bool.(df[!, :jump_feasible])
simp_converged = Bool.(df[!, :simp_converged])
full_converged = Bool.(df[!, :full_converged])

all_converged = nlsolve_feasible .& simp_converged .& full_converged
converged_scenarios = df[all_converged, :scenario_id]
n_scenarios = length(converged_scenarios)

println("📊 Found $n_scenarios scenarios where all three methods converged")
println("   Scenarios: ", join(converged_scenarios[1:min(10, n_scenarios)], ", "), 
        n_scenarios > 10 ? ", ..." : "")
println()

if n_scenarios == 0
    println("❌ No scenarios where all three converged!")
    exit(1)
end

# Monte Carlo parameters (must match original simulation)
Random.seed!(42)
MAX_LOAD_VARIATION = 0.20
MAX_OUTAGES = 2
CONTINGENCY_PROB = 0.30

# Function to reproduce a scenario
function reproduce_scenario(scenario_id::Int)
    Random.seed!(42)
    
    for sid in 1:scenario_id
        sys = build_ieee14_acdc()
        
        # Apply random load variation (ACBus is immutable, create new instances)
        load_factor = (2.0 * rand() - 1.0) * MAX_LOAD_VARIATION
        new_ac_buses = []
        for bus in sys.ac_buses
            variation = 1.0 + load_factor
            new_bus = HybridACDCPowerFlow.ACBus(
                bus.id,
                bus.type,
                bus.Pd * variation,
                bus.Qd * variation,
                bus.Pg,
                bus.Qg,
                bus.Vm,
                bus.Va,
                bus.area
            )
            push!(new_ac_buses, new_bus)
        end
        sys.ac_buses = new_ac_buses
        
        # Apply random line outages (ACBranch is immutable, create new instances)
        n_branches = length(sys.ac_branches)
        new_ac_branches = []
        n_outages = 0
        for (i, branch) in enumerate(sys.ac_branches)
            if rand() < CONTINGENCY_PROB && n_outages < MAX_OUTAGES
                new_branch = HybridACDCPowerFlow.ACBranch(
                    branch.from,
                    branch.to,
                    branch.r,
                    branch.x,
                    branch.b,
                    branch.tap,
                    false  # status = false for outage
                )
                n_outages += 1
            else
                new_branch = branch
            end
            push!(new_ac_branches, new_branch)
        end
        sys.ac_branches = new_ac_branches
        
        # Apply random converter modes (VSCConverter is also immutable, create new instances)
        new_converters = []
        for conv in sys.converters
            mode_val = rand(1:3)
            new_mode = [PQ_MODE, VDC_Q, VDC_VAC][mode_val]
            new_conv = HybridACDCPowerFlow.VSCConverter(
                conv.id,
                conv.ac_bus,
                conv.dc_bus,
                new_mode,
                conv.Pset,
                conv.Qset,
                conv.Vdc_set,
                conv.Vac_set,
                conv.Ploss_a,
                conv.Ploss_b,
                conv.Ploss_c,
                conv.Smax,
                conv.status; G_droop=conv.G_droop
            )
            push!(new_converters, new_conv)
        end
        sys.converters = new_converters
        
        # Rebuild admittance matrix with new configuration
        sys.Ybus = build_admittance_matrix(sys)
        
        if sid == scenario_id
            return sys
        end
    end
end

# Storage for comparison results
comparison_results = []

println("Running three-way comparison on $n_scenarios scenarios...")
println("(This may take a few minutes as we re-run NLsolve)")
println()

# Progress indicator
n_sample = min(n_scenarios, 20)  # Analyze first 20 or all if fewer
println("Analyzing first $n_sample scenarios for detailed comparison...")
println()

for (idx, scenario_id) in enumerate(converged_scenarios[1:n_sample])
    print("  Scenario $scenario_id ($idx/$n_sample)...")
    
    # Reproduce scenario
    sys = reproduce_scenario(scenario_id)
    
    # Run all three methods
    nlsolve_result = check_power_flow_feasibility(sys, method=:nlsolve, verbose=false, max_iter=100)
    
    dist_slack = create_participation_factors(sys; method=:capacity)
    simp_result = solve_power_flow_distributed_slack(sys, dist_slack; verbose=false, max_iter=100, tol=1e-8)
    full_result = solve_power_flow_distributed_slack_full(sys, dist_slack; verbose=false, max_iter=100, tol=1e-8)
    
    if !nlsolve_result.feasible || !simp_result.converged || !full_result.converged
        println(" SKIP (convergence mismatch)")
        continue
    end
    
    # Extract solutions
    nlsolve_Vm = nlsolve_result.Vm
    nlsolve_Va = nlsolve_result.Va
    simp_Vm = simp_result.Vm
    simp_Va = simp_result.Va
    full_Vm = full_result.Vm
    full_Va = full_result.Va
    
    # Compute three-way comparisons
    # NLsolve vs Simplified
    Vm_diff_nlsolve_simp = norm(nlsolve_Vm - simp_Vm, Inf)
    Vm_rms_nlsolve_simp = sqrt(mean((nlsolve_Vm - simp_Vm).^2))
    Va_diff_nlsolve_simp = norm(nlsolve_Va - simp_Va, Inf)
    
    # NLsolve vs Full
    Vm_diff_nlsolve_full = norm(nlsolve_Vm - full_Vm, Inf)
    Vm_rms_nlsolve_full = sqrt(mean((nlsolve_Vm - full_Vm).^2))
    Va_diff_nlsolve_full = norm(nlsolve_Va - full_Va, Inf)
    
    # Simplified vs Full (for reference)
    Vm_diff_simp_full = norm(simp_Vm - full_Vm, Inf)
    Vm_rms_simp_full = sqrt(mean((simp_Vm - full_Vm).^2))
    Va_diff_simp_full = norm(simp_Va - full_Va, Inf)
    
    push!(comparison_results, (
        scenario = scenario_id,
        # NLsolve vs Simplified
        nlsolve_simp_Vm_max = Vm_diff_nlsolve_simp,
        nlsolve_simp_Vm_rms = Vm_rms_nlsolve_simp,
        nlsolve_simp_Va_max = Va_diff_nlsolve_simp,
        # NLsolve vs Full
        nlsolve_full_Vm_max = Vm_diff_nlsolve_full,
        nlsolve_full_Vm_rms = Vm_rms_nlsolve_full,
        nlsolve_full_Va_max = Va_diff_nlsolve_full,
        # Simplified vs Full
        simp_full_Vm_max = Vm_diff_simp_full,
        simp_full_Vm_rms = Vm_rms_simp_full,
        simp_full_Va_max = Va_diff_simp_full,
        # Residuals
        nlsolve_residual = nlsolve_result.objective,
        simp_residual = simp_result.residual,
        full_residual = full_result.residual,
        # Iterations
        nlsolve_iters = nlsolve_result.iterations,
        simp_iters = simp_result.iterations,
        full_iters = full_result.iterations
    ))
    
    println(" ✓")
end

println()
println("="^80)
println("  THREE-WAY COMPARISON RESULTS")
println("="^80)
println()

# Convert to arrays for statistics
nlsolve_simp_Vm_max = [r.nlsolve_simp_Vm_max for r in comparison_results]
nlsolve_simp_Vm_rms = [r.nlsolve_simp_Vm_rms for r in comparison_results]
nlsolve_full_Vm_max = [r.nlsolve_full_Vm_max for r in comparison_results]
nlsolve_full_Vm_rms = [r.nlsolve_full_Vm_rms for r in comparison_results]
simp_full_Vm_max = [r.simp_full_Vm_max for r in comparison_results]
simp_full_Vm_rms = [r.simp_full_Vm_rms for r in comparison_results]

println("1️⃣  VOLTAGE MAGNITUDE DISCREPANCIES")
println("-"^80)
println()

@printf("  NLsolve vs Simplified NR:\n")
@printf("    Max difference:  %.6f ± %.6f p.u. (%.3f%%)\n", 
        mean(nlsolve_simp_Vm_max), std(nlsolve_simp_Vm_max), 100*mean(nlsolve_simp_Vm_max))
@printf("    RMS difference:  %.6f ± %.6f p.u.\n",
        mean(nlsolve_simp_Vm_rms), std(nlsolve_simp_Vm_rms))
@printf("    Range: [%.6f, %.6f] p.u.\n", minimum(nlsolve_simp_Vm_max), maximum(nlsolve_simp_Vm_max))
println()

@printf("  NLsolve vs Full Jacobian:\n")
@printf("    Max difference:  %.6f ± %.6f p.u. (%.3f%%)\n",
        mean(nlsolve_full_Vm_max), std(nlsolve_full_Vm_max), 100*mean(nlsolve_full_Vm_max))
@printf("    RMS difference:  %.6f ± %.6f p.u.\n",
        mean(nlsolve_full_Vm_rms), std(nlsolve_full_Vm_rms))
@printf("    Range: [%.6f, %.6f] p.u.\n", minimum(nlsolve_full_Vm_max), maximum(nlsolve_full_Vm_max))
println()

@printf("  Simplified vs Full Jacobian (reference):\n")
@printf("    Max difference:  %.6f ± %.6f p.u. (%.3f%%)\n",
        mean(simp_full_Vm_max), std(simp_full_Vm_max), 100*mean(simp_full_Vm_max))
@printf("    RMS difference:  %.6f ± %.6f p.u.\n",
        mean(simp_full_Vm_rms), std(simp_full_Vm_rms))
@printf("    Range: [%.6f, %.6f] p.u.\n", minimum(simp_full_Vm_max), maximum(simp_full_Vm_max))
println()

# Which NR method is closer to NLsolve?
nlsolve_closer_to_simp = sum(nlsolve_simp_Vm_max .< nlsolve_full_Vm_max)
nlsolve_closer_to_full = sum(nlsolve_full_Vm_max .< nlsolve_simp_Vm_max)

@printf("  NLsolve solution is closer to:\n")
@printf("    Simplified NR:   %d cases (%.1f%%)\n", 
        nlsolve_closer_to_simp, 100*nlsolve_closer_to_simp/length(comparison_results))
@printf("    Full Jacobian:   %d cases (%.1f%%)\n",
        nlsolve_closer_to_full, 100*nlsolve_closer_to_full/length(comparison_results))
println()

# Voltage angle comparison
nlsolve_simp_Va_max = [r.nlsolve_simp_Va_max for r in comparison_results]
nlsolve_full_Va_max = [r.nlsolve_full_Va_max for r in comparison_results]
simp_full_Va_max = [r.simp_full_Va_max for r in comparison_results]

println("2️⃣  VOLTAGE ANGLE DISCREPANCIES")
println("-"^80)
println()

@printf("  NLsolve vs Simplified NR:\n")
@printf("    Max difference:  %.6f rad (%.2f°)\n", 
        mean(nlsolve_simp_Va_max), rad2deg(mean(nlsolve_simp_Va_max)))
@printf("    Range: [%.6f, %.6f] rad\n", minimum(nlsolve_simp_Va_max), maximum(nlsolve_simp_Va_max))
println()

@printf("  NLsolve vs Full Jacobian:\n")
@printf("    Max difference:  %.6f rad (%.2f°)\n",
        mean(nlsolve_full_Va_max), rad2deg(mean(nlsolve_full_Va_max)))
@printf("    Range: [%.6f, %.6f] rad\n", minimum(nlsolve_full_Va_max), maximum(nlsolve_full_Va_max))
println()

@printf("  Simplified vs Full Jacobian:\n")
@printf("    Max difference:  %.6f rad (%.2f°)\n",
        mean(simp_full_Va_max), rad2deg(mean(simp_full_Va_max)))
@printf("    Range: [%.6f, %.6f] rad\n", minimum(simp_full_Va_max), maximum(simp_full_Va_max))
println()

# Residual comparison
nlsolve_residuals = [r.nlsolve_residual for r in comparison_results]
simp_residuals = [r.simp_residual for r in comparison_results]
full_residuals = [r.full_residual for r in comparison_results]

println("3️⃣  CONVERGENCE QUALITY (Residuals)")
println("-"^80)
println()

@printf("  NLsolve:         %.3e ± %.3e (median: %.3e)\n",
        mean(nlsolve_residuals), std(nlsolve_residuals), median(nlsolve_residuals))
@printf("  Simplified NR:   %.3e ± %.3e (median: %.3e)\n",
        mean(simp_residuals), std(simp_residuals), median(simp_residuals))
@printf("  Full Jacobian:   %.3e ± %.3e (median: %.3e)\n",
        mean(full_residuals), std(full_residuals), median(full_residuals))
println()

# Best residual
best_nlsolve = sum((nlsolve_residuals .< simp_residuals) .& (nlsolve_residuals .< full_residuals))
best_simp = sum((simp_residuals .< nlsolve_residuals) .& (simp_residuals .< full_residuals))
best_full = sum((full_residuals .< nlsolve_residuals) .& (full_residuals .< simp_residuals))

@printf("  Best residual achieved by:\n")
@printf("    NLsolve:       %d cases (%.1f%%)\n", best_nlsolve, 100*best_nlsolve/length(comparison_results))
@printf("    Simplified:    %d cases (%.1f%%)\n", best_simp, 100*best_simp/length(comparison_results))
@printf("    Full:          %d cases (%.1f%%)\n", best_full, 100*best_full/length(comparison_results))
println()

# Iteration comparison
nlsolve_iters = [r.nlsolve_iters for r in comparison_results]
simp_iters = [r.simp_iters for r in comparison_results]
full_iters = [r.full_iters for r in comparison_results]

println("4️⃣  ITERATION COUNT")
println("-"^80)
println()

@printf("  NLsolve:         %.1f ± %.1f (range: %d - %d)\n",
        mean(nlsolve_iters), std(nlsolve_iters), minimum(nlsolve_iters), maximum(nlsolve_iters))
@printf("  Simplified NR:   %.1f ± %.1f (range: %d - %d)\n",
        mean(simp_iters), std(simp_iters), minimum(simp_iters), maximum(simp_iters))
@printf("  Full Jacobian:   %.1f ± %.1f (range: %d - %d)\n",
        mean(full_iters), std(full_iters), minimum(full_iters), maximum(full_iters))
println()

# Top discrepancy cases
println("5️⃣  TOP DISCREPANCY CASES")
println("-"^80)
println()

# Sort by NLsolve-Simplified difference
sorted_nlsolve_simp = sort(comparison_results, by=r->r.nlsolve_simp_Vm_max, rev=true)
println("  NLsolve vs Simplified (largest voltage differences):")
for i in 1:min(5, length(sorted_nlsolve_simp))
    r = sorted_nlsolve_simp[i]
    @printf("    Scenario %3d: Vm_diff = %.6f p.u. (%.3f%%), Va_diff = %.4f rad\n",
            r.scenario, r.nlsolve_simp_Vm_max, 100*r.nlsolve_simp_Vm_max, r.nlsolve_simp_Va_max)
end
println()

# Sort by NLsolve-Full difference
sorted_nlsolve_full = sort(comparison_results, by=r->r.nlsolve_full_Vm_max, rev=true)
println("  NLsolve vs Full Jacobian (largest voltage differences):")
for i in 1:min(5, length(sorted_nlsolve_full))
    r = sorted_nlsolve_full[i]
    @printf("    Scenario %3d: Vm_diff = %.6f p.u. (%.3f%%), Va_diff = %.4f rad\n",
            r.scenario, r.nlsolve_full_Vm_max, 100*r.nlsolve_full_Vm_max, r.nlsolve_full_Va_max)
end
println()

# Relative differences
println("6️⃣  RELATIVE COMPARISON")
println("-"^80)
println()

# Ratio of differences
ratio_nlsolve_vs_simp_full = nlsolve_simp_Vm_max ./ simp_full_Vm_max
ratio_nlsolve_vs_simp_full2 = nlsolve_full_Vm_max ./ simp_full_Vm_max

@printf("  How does NLsolve compare to NR method differences?\n")
@printf("    NLsolve-Simp difference / Simp-Full difference: %.2fx\n", mean(ratio_nlsolve_vs_simp_full))
@printf("    NLsolve-Full difference / Simp-Full difference: %.2fx\n", mean(ratio_nlsolve_vs_simp_full2))
println()

if mean(ratio_nlsolve_vs_simp_full) > 2.0 || mean(ratio_nlsolve_vs_simp_full2) > 2.0
    println("  ⚠️  NLsolve produces DIFFERENT solutions than both NR methods")
    println("      This suggests multiple valid power flow solutions exist")
elseif mean(ratio_nlsolve_vs_simp_full) < 1.5 && mean(ratio_nlsolve_vs_simp_full2) < 1.5
    println("  ✅ NLsolve produces SIMILAR solutions to NR methods")
    println("      All three methods converge to the same equilibrium")
else
    println("  ℹ️  NLsolve produces MODERATELY different solutions")
    println("      Expected due to different numerical paths")
end
println()

# Save detailed results
println("7️⃣  SAVING RESULTS")
println("-"^80)

output_df = DataFrame(comparison_results)
output_file = "results/monte_carlo/csv/three_way_comparison.csv"
if isdir("results/monte_carlo/csv")
    CSV.write(output_file, output_df)
    println("  ✅ Detailed three-way comparison saved to:")
    println("     $output_file")
else
    CSV.write("three_way_comparison.csv", output_df)
    println("  ✅ Saved to: three_way_comparison.csv")
end
println()

println("="^80)
println("  SUMMARY")
println("="^80)
println()

avg_nlsolve_simp = mean(nlsolve_simp_Vm_max)
avg_nlsolve_full = mean(nlsolve_full_Vm_max)
avg_simp_full = mean(simp_full_Vm_max)

@printf("Average voltage discrepancies:\n")
@printf("  NLsolve ↔ Simplified:  %.4f p.u. (%.2f%%)\n", avg_nlsolve_simp, 100*avg_nlsolve_simp)
@printf("  NLsolve ↔ Full:        %.4f p.u. (%.2f%%)\n", avg_nlsolve_full, 100*avg_nlsolve_full)
@printf("  Simplified ↔ Full:     %.4f p.u. (%.2f%%)\n", avg_simp_full, 100*avg_simp_full)
println()

# Interpretation
if avg_nlsolve_simp < avg_simp_full && avg_nlsolve_full < avg_simp_full
    println("✅ NLsolve is CLOSER to both NR methods than they are to each other")
    println("   → NLsolve finds a \"middle ground\" solution")
elseif avg_nlsolve_simp < avg_nlsolve_full
    println("✅ NLsolve is CLOSER to Simplified NR")
    println("   → NLsolve and Simplified converge to similar equilibrium")
elseif avg_nlsolve_full < avg_nlsolve_simp
    println("✅ NLsolve is CLOSER to Full Jacobian")
    println("   → NLsolve and Full Jacobian converge to similar equilibrium")
else
    println("ℹ️  All three methods produce comparable solutions")
    println("   → Different numerical paths, similar physical results")
end
println()

if avg_nlsolve_simp < 0.01 && avg_nlsolve_full < 0.01
    println("🎯 EXCELLENT AGREEMENT: All three methods produce nearly identical results (<1%)")
elseif avg_nlsolve_simp < 0.02 && avg_nlsolve_full < 0.02
    println("✓  GOOD AGREEMENT: All differences within acceptable tolerance (<2%)")
else
    println("⚠️  MODERATE DIFFERENCES: May indicate multiple equilibrium points")
    println("   → All solutions are valid, just different")
end

println("\n" * "="^80)
