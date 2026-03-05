# HybridACDCPowerFlow Documentation (v0.5.0)

This folder contains the current technical documentation for the functional layer package used by JuliaPowerCase.

## Documentation Map

- `quickstart.md`: installation and first solves
- `architecture.md`: module and data-flow design
- `api_reference.md`: public API and key signatures
- `solver_core.md`: Newton-Raphson solver details
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
using HybridACDCPowerFlow
sys = build_ieee14_acdc()
result = solve_power_flow(sys)
@show result.converged result.iterations result.residual
```

## Version Scope

The docs in this folder match code state `v0.5.0` (JuliaPowerCase integration).
