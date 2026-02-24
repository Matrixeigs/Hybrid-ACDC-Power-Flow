# HybridACDCPowerFlow.jl v0.4.0

A Julia module for hybrid AC/DC power flow analysis with VSC converter models, enhanced features for resilience analysis, and distributed slack bus modeling.

**Latest**: v0.4.0 adds **Optimized Sparse Jacobian**, **Pre-allocated SolverWorkspace**, **Package Extensions**, and **MATPOWER Test Systems**!

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
- ✅ **Validated**: Full theoretical foundation with Monte Carlo testing

### NEW - Performance Optimizations (v0.4.0)
- ⚡ **Sparse Jacobian**: Pre-computed sparsity pattern with COO triplet construction
- 🧠 **SolverWorkspace**: Pre-allocated buffers for zero-allocation NR iterations
- 🔧 **Sparse LU**: UMFPACK factorization with symbolic reuse
- 📦 **Package Extensions**: JuMP/Ipopt/NLsolve moved to optional `ext/FeasibilityExt.jl`
- 🔌 **MATPOWER Systems**: case33bw, case33mg, case69, case300, case2000 distribution networks
- 📁 **MatpowerParser**: Load standard MATPOWER `.m` files
- 📈 **Performance**: 2-7× speedup on large systems (case2000: 180ms → 25ms)

## Installation

```julia
using Pkg
Pkg.develop(path="path/to/HybridACDCPowerFlow")

# Optional: Enable feasibility checking
Pkg.add(["JuMP", "Ipopt", "NLsolve"])
```

## Quick Start

### Basic Power Flow
```julia
using HybridACDCPowerFlow

# Build IEEE 14-bus with HVDC
sys = build_ieee14_acdc()

# Solve power flow
result = solve_power_flow(sys)

# Check convergence
println("Converged: $(result.converged)")
println("Voltage: $(result.Vm)")

# Extract results
Vm, Va, Vdc = get_bus_voltages(result)
Pij, Qij = get_branch_flows(sys, result)
```

### Optimized Solver (v0.4.0)
```julia
# Create pre-allocated workspace (do once)
ws = create_solver_workspace(sys)

# Solve with zero allocation in NR loop
result = solve_power_flow(sys; workspace=ws)

# Reuse workspace for multiple solves (same topology)
for scenario in scenarios
    apply_load_variation!(sys, scenario)
    result = solve_power_flow(sys; workspace=ws)
end
```

### MATPOWER Test Systems (v0.4.0)
```julia
# Distribution network test cases
sys33 = build_case33bw_acdc()   # 33-bus radial
sys69 = build_case69_acdc()     # 69-bus radial
sys300 = build_case300_acdc()   # 300-bus transmission
sys2000 = build_case2000_acdc() # 2000-bus scalability test
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
# Multi-generator balancing with different methods
dist_slack = create_participation_factors(sys; method=:capacity)  # or :droop, :equal

result = solve_power_flow_distributed_slack(sys, dist_slack; verbose=true)

# View slack distribution
for (bus, P_slack) in result.distributed_slack_P
    println("Bus $bus: ΔP = $P_slack p.u.")
end
```

### Full Jacobian with Limits
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

### Feasibility Checking (Extension Required)
```julia
using HybridACDCPowerFlow
using JuMP, Ipopt, NLsolve  # Triggers extension loading

# Check if scenario is feasible
result = check_power_flow_feasibility(sys; method=:nlsolve, verbose=true)

println("Feasible: ", result.feasible)
println("Load margin: ", result.load_margin, " p.u.")
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
- **[docs/HybridACDCPowerFlow_CompleteDocumentation.pdf](docs/HybridACDCPowerFlow_CompleteDocumentation.pdf)** - Complete LaTeX documentation
- **[docs/TheoreticalFoundations.pdf](docs/TheoreticalFoundations.pdf)** - Mathematical theory (15 pages)
- **[QUICK_REFERENCE.md](QUICK_REFERENCE.md)** - One-page cheat sheet
- **[docs/FULL_JACOBIAN_IMPLEMENTATION.md](docs/FULL_JACOBIAN_IMPLEMENTATION.md)** - Full Jacobian implementation
- **[docs/FEASIBILITY_CHECKER.md](docs/FEASIBILITY_CHECKER.md)** - JuMP-based feasibility checker

### 🧪 Test Suites
- **`test/runtests.jl`** - Main test entry point
- **`test/test_enhanced.jl`** - Enhanced features tests
- **`test/test_distributed_slack.jl`** - Distributed slack tests
- **`test_optimized_solver.jl`** - Sparse Jacobian & workspace tests
- **`test_extension.jl`** - Package extension tests

## Project Structure

```
HybridACDCPowerFlow/
├── Project.toml              # v0.4.0 with [weakdeps] and [extensions]
├── src/
│   ├── HybridACDCPowerFlow.jl
│   ├── PowerSystem.jl        # Optimized core solver, SolverWorkspace
│   ├── PowerSystemEnhanced.jl
│   ├── TestSystems.jl        # IEEE + MATPOWER test cases
│   ├── MatpowerParser.jl     # MATPOWER file parser
│   └── FeasibilityChecker.jl # Stubs for extension
├── ext/
│   └── FeasibilityExt.jl     # Optional JuMP/Ipopt/NLsolve
├── data/                     # MATPOWER case files
├── docs/
└── test/
```

## Performance Benchmarks (v0.4.0)

| System | Dense (v0.3) | Sparse (v0.4) | Speedup |
|--------|--------------|---------------|---------|
| IEEE 14-bus | 0.12 ms | 0.07 ms | 1.7× |
| IEEE 24-bus | 0.18 ms | 0.12 ms | 1.5× |
| IEEE 118-bus | 2.1 ms | 0.8 ms | 2.6× |
| case300 | 8.5 ms | 2.1 ms | 4.0× |
| case2000 | 180 ms | 25 ms | 7.2× |

## Version History

| Version | Date | Features |
|---------|------|----------|
| 0.1.0 | Nov 2025 | Core AC/DC power flow, VSC models |
| 0.2.0 | Feb 2026 | Islands, PV→PQ, auto swing, mode switching |
| 0.3.0 | Feb 2026 | Distributed slack bus model |
| 0.3.1 | Feb 2026 | Full Jacobian + Monte Carlo validation |
| **0.4.0** | **Feb 2026** | **Sparse Jacobian, SolverWorkspace, Extensions, MATPOWER** |
