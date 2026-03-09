# HybridACDCPowerFlow Documentation (v0.7.0)

This folder contains the current technical documentation for the application layer package that works with JuliaPowerCase data structures.

## Architecture

**HybridACDCPowerFlow** is now a **pure application layer** that:
- Accepts data from **JuliaPowerCase** (information/data layer)
- Performs power flow computations using Newton-Raphson solvers
- Returns results without modifying input data structures

**Key Design Principle**: `JuliaPowerCase` provides `HybridPowerSystem` (data), `HybridACDCPowerFlow` provides algorithms.

## Documentation Map

- `quickstart.md`: installation and first solves
- `architecture.md`: module and data-flow design
- `api_reference.md`: public API and key signatures
- `solver_core.md`: Newton-Raphson solver details, DC bus types, pure DC solver
- `adaptive_and_islanding.md`: island detection and adaptive solving
- `distributed_slack.md`: distributed slack models
- `juliapowercase_integration.md`: HybridPowerSystem and HybridPowerCaseData workflows
- `matpower_and_cases.md`: MATPOWER parser and built-in benchmark systems
- `feasibility_extension.md`: optional JuMP/Ipopt/NLsolve extension
- `testing_and_validation.md`: test suites and validation coverage
- `troubleshooting.md`: common failure modes and fixes
- `code_quality_assessment.md`: current reassessment and improvement backlog

## Fast Start

```bash
julia --project=HybridACDCPowerFlow HybridACDCPowerFlow/test/runtests.jl
julia --project=HybridACDCPowerFlow HybridACDCPowerFlow/test/test_jpc_integration.jl
```

```julia
using JuliaPowerCase  # Data structures
using HybridACDCPowerFlow  # Algorithms

# Hybrid AC/DC power flow using JuliaPowerCase data
sys = build_ieee14_acdc()  # Returns HybridPowerSystem
result = solve_power_flow(sys)
@show result.converged result.iterations result.residual

# Pure DC power flow (DC_V = reference, DC_P = power-specified)
dc_sys = # ... HybridPowerSystem with only DC buses ...
dc_result = solve_dc_power_flow(dc_sys)
@show dc_result.Vdc
```

## Version Scope

The docs in this folder match code state `v0.7.0` (refactored to work with JuliaPowerCase v1.0 data structures).
