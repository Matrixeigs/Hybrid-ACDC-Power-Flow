# HybridACDCPowerFlow.jl v0.3.1

A Julia module for hybrid AC/DC power flow analysis with VSC converter models, enhanced features for resilience analysis, and distributed slack bus modeling.

**Latest**: v0.3.1 adds **Full Jacobian Distributed Slack** with capacity limit enforcement and **Monte Carlo Stochastic Validation**!

## Features

### Core Power Flow (v0.1.0)
- **AC Power Flow**: Newton-Raphson solver for AC networks
- **DC Network**: Linear DC power flow analysis  
- **VSC Converters**: Multiple control modes (PQ, VDC_Q, VDC_VAC)
- **Test Systems**: IEEE 14, 24, 118-bus with HVDC extensions
- **Theory Validation**: Comprehensive tests for PE-GNN theoretical claims

### Enhanced Features (v0.2.0)
- 🏝️ **Island Detection**: Graph-based DFS algorithm for N-k contingencies
- ⚡ **PV→PQ Conversion**: Automatic reactive limit enforcement
- 🎯 **Auto Swing Bus**: Multi-criteria slack bus selection
- 🔄 **Converter Switching**: Intelligent mode switching based on voltage/loading
- 📚 **Theory Docs**: Complete LaTeX mathematical derivations

### Distributed Slack (v0.3.0)
- 🌐 **Multi-Generator Balancing**: Realistic AGC simulation
- 📊 **Participation Factors**: Capacity/droop/equal-based methods
- 🔗 **Island Integration**: Per-island distributed slack configuration
- ✅ **Validated**: 21/21 tests passing, full theoretical foundation

### NEW - Full Jacobian & Stochastic Validation (v0.3.1)
- 🧮 **Full Jacobian Formulation**: Dimensionally correct distributed slack implementation
- ⚡ **Capacity Limit Enforcement**: During iterations (not post-processing)
- 📊 **Monte Carlo Framework**: Comprehensive stochastic power flow testing
- 📈 **Statistical Analysis**: Random loads, N-k contingencies, control modes
- 🔍 **JuMP Feasibility Checker**: Pre-screen scenarios for infeasibility (Ipopt-based)
- ✅ **Validated**: 87/87 tests passing, 100-scenario Monte Carlo simulation

## Installation

```julia
using Pkg
Pkg.develop(path="path/to/HybridACDCPowerFlow")
```

## Quick Start

### Basic Power Flow (v0.1.0)
```julia
using HybridACDCPowerFlow

# Build IEEE 14-bus with HVDC
sys = build_ieee14_acdc()

# Solve power flow
result = solve_power_flow(sys)

# Check convergence
println("Converged: $(result.converged)")
println("Voltage: $(result.Vm)")
```

### Enhanced Features (v0.2.0)
```julia
# Island detection + adaptive power flow
islands = detect_islands(sys)
Q_limits = create_default_Q_limits(sys)

result = solve_power_flow_adaptive(sys; 
    options=PowerFlowOptions(
        enable_pv_pq_conversion=true,
        enable_auto_swing_selection=true,
        enable_converter_mode_switching=true
    ),
    Q_limits=Q_limits
)
```

### Distributed Slack (v0.3.0)
```julia
# Multi-generator balancing
dist_slack = create_participation_factors(sys; method=:capacity)

result = solve_power_flow_distributed_slack(sys, dist_slack; verbose=true)

# View slack distribution
for (bus, P_slack) in result.distributed_slack_P
    println("Bus $bus: ΔP = $P_slack p.u.")
end
```

### Full Jacobian with Limits (v0.3.1)
```julia
# Create distributed slack with capacity limits
dist_slack = DistributedSlack(
    [1, 2, 6],              # Participating buses
    [0.4, 0.4, 0.2],        # Participation factors
    1,                       # Reference bus
    Dict(2 => 0.05)         # Capacity limit on bus 2
)

# Solve with full Jacobian (enforces limits during iterations)
result = solve_power_flow_distributed_slack_full(sys, dist_slack;
    enforce_limits=true, verbose=true)
```

### Monte Carlo Stochastic Testing (v0.3.1)
```bash
# Run 100 random scenarios (±20% load, N-k contingencies)
julia --project=. test/stochastic_power_flow_monte_carlo.jl

# Analyze results statistically
julia --project=. test/analyze_monte_carlo_results.jl

# Generate visualizations (requires Plots.jl)
julia --project=. test/visualize_monte_carlo_results.jl
```

### JuMP Feasibility Checker (v0.3.1)
```julia
# Check if scenario is feasible before running NR solvers
result = check_power_flow_feasibility(sys; verbose=true, max_time=10.0)

println("Feasible: ", result.feasible)
println("Load margin: ", result.load_margin, " p.u.")
println("Solve time: ", result.solve_time, " s")

# Automatically integrated in Monte Carlo simulation
# Results include load margin and JuMP feasibility status
```

## Test Suite

Run the comprehensive validation tests:

```julia
using Pkg
Pkg.test("HybridACDCPowerFlow")
```

This validates all 5 theoretical claims from the NSFC proposal:
1. Topology generalization (N-k faults)
2. SO(2) equivariance (multi-reference point)
3. Current injection universality (meshed/radial)
4. DC shortcut effect (graph diameter reduction)
5. Control mode robustness (unknown parameters)

## Documentation

### 📚 Complete Documentation
- **[docs/FEASIBILITY_CHECKER.md](docs/FEASIBILITY_CHECKER.md)** - v0.3.1 JuMP-based feasibility checker
- **[docs/MONTE_CARLO_RESULTS.md](docs/MONTE_CARLO_RESULTS.md)** - v0.3.1 Monte Carlo validation results
- **[DISTRIBUTED_SLACK.md](DISTRIBUTED_SLACK.md)** - v0.3.0 distributed slack bus model guide
- **[ENHANCED_FEATURES.md](ENHANCED_FEATURES.md)** - v0.2.0 enhanced features overview
- **[QUICK_REFERENCE.md](QUICK_REFERENCE.md)** - One-page cheat sheet
- **[CHANGELOG_v0.2.0.md](CHANGELOG_v0.2.0.md)** - Detailed v0.2.0 changelog
- **[docs/TheoreticalFoundations.pdf](docs/TheoreticalFoundations.pdf)** - Complete mathematical theory (15 pages)
- **[docs/FULL_JACOBIAN_IMPLEMENTATION.md](docs/FULL_JACOBIAN_IMPLEMENTATION.md)** - v0.3.1 full Jacobian implementation
- **[docs/MONTE_CARLO_ANALYSIS.md](docs/MONTE_CARLO_ANALYSIS.md)** - v0.3.1 stochastic validation guide
- **[docs/IMPLEMENTATION_SUMMARY.md](docs/IMPLEMENTATION_SUMMARY.md)** - Quick reference for v0.3.1 features

### 🧪 Test Suites
- **`test/runtests.jl`** - Original theorem validation (32 tests)
- **`test/test_enhanced.jl`** - Enhanced features (19 tests)
- **`test/test_distributed_slack.jl`** - Distributed slack (55 tests: 21 simplified + 34 full Jacobian)
- **`test/stochastic_power_flow_monte_carlo.jl`** - Monte Carlo simulation (100 scenarios)
- **`test/analyze_monte_carlo_results.jl`** - Statistical analysis
- **`test/visualize_monte_carlo_results.jl`** - Result visualization

**Total**: 87/87 tests passing ✅ | Monte Carlo: 100 scenarios validated

## Version History

| Version | Date | Features | Tests |
|---------|------|----------|-------|
| 0.1.0 | Nov 2025 | Core AC/DC power flow, VSC models | 32/32 ✅ |
| 0.2.0 | Feb 2026 | Islands, PV→PQ, auto swing, mode switching | 19/19 ✅ |
| 0.3.0 | Feb 2026 | Distributed slack bus model (simplified) | 21/21 ✅ |
| 0.3.1 | Feb 2026 | **Full Jacobian + Monte Carlo validation** | 87/87 ✅ + 100 scenarios |
