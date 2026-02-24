# Monte Carlo Stochastic Power Flow Analysis

## Overview

This document describes the comprehensive Monte Carlo simulation framework for validating and comparing the two distributed slack bus implementations:

1. **Simplified Method** (`solve_power_flow_distributed_slack`): Post-processing approach
2. **Full Jacobian Method** (`solve_power_flow_distributed_slack_full`): Full matrix formulation with capacity enforcement

## Methodology

### Stochastic Scenario Generation

The Monte Carlo simulation generates random power system scenarios by varying three key factors:

#### 1. Load Variations
- **Range**: ±20% of base case load
- **Distribution**: Uniform random distribution
- **Formula**: `P_d,new = P_d,base × (1 + λ)` where `λ ∈ [-0.2, 0.2]`

#### 2. N-k Contingencies
- **Contingency Types**: N-0 (no outage), N-1 (single contingency), N-2 (double contingency)
- **Probability**: 10% chance per transmission line
- **Maximum Outages**: 2 simultaneous line outages
- **Selection**: Random uniform selection from available lines

#### 3. Converter Control Modes
- **Available Modes**:
  - `PQ_MODE`: Fixed active and reactive power
  - `VDC_Q`: DC voltage control with reactive power control
  - `VDC_VAC`: DC voltage and AC voltage control
- **Assignment**: Random uniform selection for each converter

### Simulation Parameters

```julia
# Default Monte Carlo Configuration
N_SCENARIOS = 100              # Number of random scenarios
LOAD_VARIATION = 0.20          # ±20% load variation
CONTINGENCY_PROB = 0.10        # 10% probability per line
MAX_OUTAGES = 2                # Maximum N-2 contingencies
RANDOM_SEED = 42               # For reproducibility
```

### Participation Factor Configuration

All simulations use a **capacity-weighted participation** strategy:

```julia
α_i = (Pg_max,i - Pg_min,i) / Σ(Pg_max,j - Pg_min,j)
```

This ensures:
- Larger generators contribute more to slack distribution
- Respects generator capacity constraints
- Physically realistic power sharing

## Comparison Metrics

### 1. Convergence Metrics

| Metric | Description | Interpretation |
|--------|-------------|----------------|
| **Convergence Rate** | % of scenarios that converge | Higher is better |
| **Residual Norm** | Final mismatch norm | Lower is better (<1e-8) |
| **Iterations** | Number of Newton iterations | Lower indicates faster convergence |

### 2. Performance Metrics

| Metric | Description | Unit |
|--------|-------------|------|
| **Execution Time** | Wall-clock time for convergence | milliseconds (ms) |
| **Time per Iteration** | Average time per Newton step | ms/iteration |
| **Speedup Factor** | Ratio of execution times | dimensionless |

### 3. Accuracy Metrics

| Metric | Formula | Interpretation |
|--------|---------|----------------|
| **Max Voltage Diff** | `max\|V_simp - V_full\|` | Largest voltage discrepancy |
| **RMS Voltage Diff** | `√(Σ(V_simp - V_full)²/N)` | Average voltage discrepancy |
| **Slack Power Diff** | `\|P_slack,simp - P_slack,full\|` | Total slack power difference |

## Usage

### Running Monte Carlo Simulation

```bash
cd HybridACDCPowerFlow
julia --project=. test/stochastic_power_flow_monte_carlo.jl
```

**Outputs:**
- `results/monte_carlo/csv/monte_carlo_results.csv` - Detailed per-scenario results
- `results/monte_carlo/csv/summary_statistics.csv` - Aggregated statistics
- Console output with progress updates

### Analyzing Results

```bash
julia --project=. test/analyze_monte_carlo_results.jl
```

**Outputs:**
- `results/monte_carlo/analysis/detailed_analysis_report.txt` - Comprehensive text report
- Console output with statistical summaries

### Generating Visualizations

```bash
# Requires: Plots, StatsPlots
julia --project=. test/visualize_monte_carlo_results.jl
```

**Outputs:** 8 publication-quality figures (PNG format):

1. **fig1_convergence_comparison.png** - Convergence rates by contingency level
2. **fig2_iteration_comparison.png** - Iteration count distributions
3. **fig3_time_comparison.png** - Execution time distributions
4. **fig4_voltage_difference_analysis.png** - Voltage accuracy metrics
5. **fig5_slack_comparison.png** - Slack power distribution comparison
6. **fig6_voltage_profile_comparison.png** - Voltage profile analysis
7. **fig7_performance_accuracy_tradeoff.png** - Performance vs accuracy
8. **fig8_summary_dashboard.png** - Comprehensive summary dashboard

## Expected Results

### Typical Performance Characteristics

Based on the IEEE 24-bus test system with 4 AC-DC converters:

#### Simplified Method
- **Iterations**: 4-6 typical (mean ≈ 5)
- **Convergence**: >95% for N-0 and N-1, >85% for N-2
- **Execution Time**: 0.5-1.5 ms typical
- **Characteristics**: 
  - Fast convergence
  - Post-processing approach
  - No iterative limit enforcement

#### Full Jacobian Method
- **Iterations**: 9-13 typical (mean ≈ 11)
- **Convergence**: >95% for N-0 and N-1, >85% for N-2
- **Execution Time**: 1.0-3.0 ms typical
- **Characteristics**:
  - Quadratic convergence
  - Iterative limit enforcement
  - Dimensionally correct formulation

### Accuracy Comparison

| Metric | Typical Range | Interpretation |
|--------|---------------|----------------|
| Max Voltage Diff | 1e-6 to 1e-4 p.u. | Excellent agreement |
| RMS Voltage Diff | 1e-7 to 1e-5 p.u. | Very high consistency |
| Slack Diff | 1e-5 to 1e-3 p.u. | Acceptable for most applications |

**Note**: Differences typically increase with:
- Higher load variations
- More severe contingencies (N-2 > N-1 > N-0)
- Tighter generator capacity limits

## Interpretation Guidelines

### When to Use Simplified Method

✅ **Recommended for:**
- **Real-time applications** requiring fast computation
- **Large-scale studies** with many scenarios (Monte Carlo, optimization)
- **Preliminary analysis** and screening studies
- Cases where generators are **far from capacity limits**

### When to Use Full Jacobian Method

✅ **Recommended for:**
- **Capacity planning** where limit enforcement is critical
- **Critical infrastructure analysis** requiring high accuracy
- **Stressed system conditions** with generators near limits
- **Detailed operational studies** requiring rigorous modeling

### Validation Criteria

The methods are considered **equivalent** if:
- ✅ Convergence rates differ by <5%
- ✅ Mean voltage difference < 1e-3 p.u.
- ✅ Mean slack difference < 0.01 p.u.
- ✅ Both converge for >90% of scenarios

The methods show **significant differences** if:
- ⚠️ Convergence rate gap >10%
- ⚠️ Mean voltage difference >1e-2 p.u.
- ⚠️ Any scenario shows voltage difference >0.05 p.u.

## Debugging Failed Scenarios

### Identifying Problematic Cases

The analysis script identifies scenarios where:
1. Only one method converges
2. Neither method converges
3. Large differences occur despite convergence

### Common Failure Modes

| Failure Pattern | Likely Cause | Mitigation |
|----------------|--------------|------------|
| Both fail at high load | Infeasible power flow | Reduce LOAD_VARIATION |
| Both fail with N-2 | Network islanding | Check topology after outages |
| Only simplified fails | Capacity violations | Expected behavior |
| Only full fails | Numerical conditioning | Review Jacobian scaling |
| Large voltage diff | Control mode mismatch | Verify converter settings |

### Accessing Detailed Logs

Each scenario's complete results are saved in `monte_carlo_results.csv`:

```bash
# Find scenario 42 details
grep "^42," results/monte_carlo/csv/monte_carlo_results.csv
```

**CSV Column Reference:**
1. scenario_id
2. load_factor  
3. n_outages
4. outage_lines (JSON array)
5. simp_converged
6. simp_iterations
7. simp_residual
8. simp_time_ms
9. simp_total_slack
10. simp_vmax
11. simp_vmin
12. full_converged
13. full_iterations
14. full_residual
15. full_time_ms
16. full_total_slack
17. full_vmax
18. full_vmin
19. voltage_diff_max
20. voltage_diff_rms
21. slack_diff

## Customization

### Adjusting Scenario Parameters

Edit `stochastic_power_flow_monte_carlo.jl`:

```julia
# Increase number of scenarios for better statistics
const N_SCENARIOS = 500  # Default: 100

# Reduce load variation for less stressed conditions
const LOAD_VARIATION = 0.10  # Default: 0.20 (±20%)

# Enable N-3 contingencies
const MAX_OUTAGES = 3  # Default: 2

# Increase contingency probability
const CONTINGENCY_PROB = 0.20  # Default: 0.10
```

### Using Different Test Systems

Modify the base case loading:

```julia
# Load different test case
include("../test/test_cases.jl")
sys_base = create_ieee118_hybrid()  # Instead of IEEE 24-bus

# Or load from external file
using PowerSystems
sys_base = System("path/to/case.m")
```

### Custom Participation Factors

Replace the default capacity-weighted factors:

```julia
# Example: Equal participation
slack_config = Dict(
    :enabled => true,
    :mode => :equal,  # All generators share equally
    :capacity_limits => true
)

# Example: Manual specification
slack_config = Dict(
    :enabled => true,
    :mode => :custom,
    :participation_factors => Dict(1 => 0.4, 2 => 0.3, 3 => 0.3),
    :capacity_limits => true
)
```

## Statistical Analysis

### Hypothesis Testing

The Monte Carlo results can be used for statistical hypothesis testing:

**Null Hypothesis (H₀)**: The two methods produce statistically equivalent results

**Test Statistic**: Mean voltage difference

**Acceptance Criteria**: 
- Mean difference < 1e-3 p.u.
- 95% quantile difference < 1e-2 p.u.

### Sample Size Determination

For 95% confidence intervals with ±1% margin of error:

```
n ≥ (1.96 × σ / E)²
```

Where:
- σ = observed standard deviation
- E = desired margin of error

**Example**: If σ(iterations) = 2, for E = 0.2:
```
n ≥ (1.96 × 2 / 0.2)² = 384 scenarios
```

### Correlation Analysis

Examine correlations between:
- Load factor ↔ Convergence rate
- Number of outages ↔ Iteration count  
- Load factor ↔ Voltage difference

Use `results/monte_carlo/csv/monte_carlo_results.csv` with standard statistical tools.

## References

### Related Documentation
- [FULL_JACOBIAN_IMPLEMENTATION.md](FULL_JACOBIAN_IMPLEMENTATION.md) - Technical details
- [IMPLEMENTATION_SUMMARY.md](IMPLEMENTATION_SUMMARY.md) - Quick reference
- [FULL_JACOBIAN_ANALYSIS.md](FULL_JACOBIAN_ANALYSIS.md) - Problem analysis

### Theoretical Foundation
- Newton-Raphson power flow: Tinney & Hart (1967)
- Distributed slack bus: MATPOWER Manual, Appendix G
- Stochastic power flow: Chen et al. (2008), Hong (1998)

### Test System
- IEEE 24-Bus Reliability Test System (RTS)
- Modified with 4 AC-DC VSC converters
- See `test/test_cases.jl` for details

## Troubleshooting

### Installation Issues

**Problem**: `ERROR: ArgumentError: Package Plots not found`

**Solution**:
```bash
julia --project=.
using Pkg
Pkg.add(["Plots", "StatsPlots"])
```

### Convergence Issues

**Problem**: Low convergence rates (<80%)

**Possible causes**:
1. Excessive load variation → Reduce LOAD_VARIATION
2. Too many outages → Reduce MAX_OUTAGES or CONTINGENCY_PROB
3. Infeasible base case → Verify base case converges first
4. Numerical issues → Check for ill-conditioned scenarios

### Performance Issues

**Problem**: Simulation takes >10 minutes

**Solutions**:
- Reduce N_SCENARIOS
- Profile code to identify bottlenecks
- Use `@time` to measure individual operations
- Consider parallel execution (see below)

### Parallel Execution

For large-scale studies:

```julia
using Distributed
addprocs(4)  # Use 4 cores

@everywhere include("../src/PowerSystemEnhanced.jl")

results = @distributed (vcat) for i in 1:N_SCENARIOS
    # Run scenario i
    [scenario_result]
end
```

## Future Enhancements

Potential extensions to the Monte Carlo framework:

1. **Uncertainty Quantification**
   - Probabilistic voltage profiles
   - Risk metrics (Value-at-Risk, Conditional VaR)
   - Confidence intervals for all metrics

2. **Advanced Scenarios**
   - Time-series load profiles
   - Renewable generation uncertainty
   - Temperature-dependent line ratings

3. **Multi-Objective Analysis**
   - Pareto frontier: speed vs accuracy
   - Trade-off visualization
   - Optimization of solver parameters

4. **Machine Learning Integration**
   - Predict convergence from scenario features
   - Identify critical scenarios
   - Surrogate modeling for faster simulation

---

**Last Updated**: 2024
**Version**: 0.3.1
**Contact**: See main README.md
