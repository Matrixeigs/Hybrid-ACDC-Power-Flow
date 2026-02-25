# HybridACDCPowerFlow.jl

A Julia package for hybrid AC/DC power flow analysis with VSC converter models, enhanced resilience features, and distributed slack bus modeling.

## Features

- **AC Power Flow** — Newton-Raphson solver with sparse Jacobian and pre-allocated workspace
- **DC Network** — Linear DC power flow analysis
- **VSC Converters** — Three control modes: `PQ_MODE`, `VDC_Q`, `VDC_VAC`
- **Island Detection** — Graph-based DFS for N-k contingencies, per-island solving
- **PV→PQ Conversion** — Automatic reactive limit enforcement
- **Distributed Slack** — Multi-generator balancing with capacity/droop/equal participation
- **MATPOWER Parser** — Load standard `.m` case files (case33bw, case69, case300, etc.)
- **Feasibility Extension** — Optional JuMP/Ipopt/NLsolve-based feasibility checking

## Installation

```julia
using Pkg

# Install from local path
Pkg.develop(path="/path/to/HybridACDCPowerFlow")

# Or from a Git URL
Pkg.add(url="https://github.com/<user>/HybridACDCPowerFlow.jl.git")

# Optional: enable feasibility checking extension
Pkg.add(["JuMP", "Ipopt", "NLsolve"])
```

## Quick Start

```julia
using HybridACDCPowerFlow

# Build and solve IEEE 14-bus with HVDC
sys = build_ieee14_acdc()
result = solve_power_flow(sys)

println("Converged: $(result.converged)")
println("Voltages:  $(result.Vm)")
```

### Pre-allocated Workspace (Zero-Allocation NR Loop)

```julia
ws = create_solver_workspace(sys)
result = solve_power_flow(sys; workspace=ws)

# Reuse workspace across many solves with the same topology
for scenario in scenarios
    result = solve_power_flow(sys; workspace=ws)
end
```

### MATPOWER Test Systems

```julia
sys33  = build_case33bw_acdc()    # 33-bus radial distribution
sys69  = build_case69_acdc()      # 69-bus radial distribution
sys300 = build_case300_acdc()     # 300-bus transmission
```

### Adaptive Power Flow with Island Detection

```julia
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

### Distributed Slack

```julia
dist_slack = create_participation_factors(sys; method=:capacity)
result = solve_power_flow_distributed_slack(sys, dist_slack; verbose=true)
```

### Feasibility Checking (Extension)

```julia
using HybridACDCPowerFlow
using JuMP, Ipopt, NLsolve   # triggers extension loading

result = check_power_flow_feasibility(sys; method=:nlsolve, verbose=true)
println("Feasible: ", result.feasible)
```

## Test Suite

```julia
using Pkg
Pkg.test("HybridACDCPowerFlow")
```

Validates 5 theoretical claims: topology generalization, SO(2) equivariance, current injection universality, DC shortcut effect, and control mode robustness.

## Project Structure

```
HybridACDCPowerFlow/
├── Project.toml                  # Package metadata, deps, extensions
├── src/
│   ├── HybridACDCPowerFlow.jl    # Main module entry point
│   ├── PowerSystem.jl            # Core AC/DC solver, SolverWorkspace
│   ├── PowerSystemEnhanced.jl    # Islands, PV→PQ, distributed slack
│   ├── MatpowerParser.jl         # MATPOWER .m file parser
│   ├── TestSystems.jl            # IEEE + MATPOWER test cases
│   └── FeasibilityChecker.jl     # Stub interface for extension
├── ext/
│   └── FeasibilityExt.jl         # Optional JuMP/Ipopt/NLsolve extension
├── data/                         # MATPOWER case files (.m)
├── test/
│   ├── runtests.jl               # Main test entry point (Pkg.test)
│   ├── test_enhanced.jl          # Enhanced features tests
│   ├── test_distributed_slack.jl # Distributed slack tests
│   └── ...                       # Monte Carlo, visualization, etc.
├── examples/                     # Benchmark & comparison scripts
├── docs/                         # PDFs, Monte Carlo analysis
├── README.md
├── QUICK_REFERENCE.md
└── LICENSE
```

## Documentation

| Document | Description |
|----------|-------------|
| [README.md](README.md) | This file — overview and quick start |
| [QUICK_REFERENCE.md](QUICK_REFERENCE.md) | One-page API cheat sheet |
| [docs/TheoreticalFoundations.pdf](docs/TheoreticalFoundations.pdf) | Mathematical theory (15 pages) |
| [docs/HybridACDCPowerFlow_CompleteDocumentation.pdf](docs/HybridACDCPowerFlow_CompleteDocumentation.pdf) | Complete LaTeX documentation |
| [docs/MONTE_CARLO_ANALYSIS.md](docs/MONTE_CARLO_ANALYSIS.md) | Monte Carlo simulation methodology |
| [docs/MONTE_CARLO_RESULTS.md](docs/MONTE_CARLO_RESULTS.md) | Monte Carlo simulation results |

## Performance Benchmarks

| System | Dense (v0.3) | Sparse (v0.4) | Speedup |
|--------|-------------|---------------|---------|
| IEEE 14-bus | 0.12 ms | 0.07 ms | 1.7× |
| IEEE 24-bus | 0.18 ms | 0.12 ms | 1.5× |
| IEEE 118-bus | 2.1 ms | 0.8 ms | 2.6× |
| case300 | 8.5 ms | 2.1 ms | 4.0× |
| case2000 | 180 ms | 25 ms | 7.2× |

## Version History

| Version | Date | Highlights |
|---------|------|------------|
| 0.1.0 | Nov 2025 | Core AC/DC power flow, VSC models |
| 0.2.0 | Feb 2026 | Islands, PV→PQ, auto swing, mode switching |
| 0.3.0 | Feb 2026 | Distributed slack bus model |
| 0.3.1 | Feb 2026 | Full Jacobian + Monte Carlo validation |
| 0.4.0 | Feb 2026 | Sparse Jacobian, SolverWorkspace, extensions, MATPOWER |

## License

MIT — see [LICENSE](LICENSE).
