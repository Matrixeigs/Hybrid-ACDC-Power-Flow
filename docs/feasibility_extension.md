# Feasibility Extension

## Extension Model

`HybridACDCPowerFlow` exposes feasibility APIs in core, but concrete implementations are loaded from `ext/FeasibilityExt.jl`.

Weak dependencies in `Project.toml`:

- `JuMP`
- `Ipopt`
- `NLsolve`

## Activation

```julia
using HybridACDCPowerFlow
using JuMP, Ipopt, NLsolve
```

Then call:

```julia
res = check_power_flow_feasibility(sys; method=:nlsolve)
```

## Methods

- `:nlsolve`
  - root-finding style feasibility check
  - returns `FeasibilityResult`
- `:jump`
  - nonlinear feasibility optimization via JuMP + Ipopt
  - returns `FeasibilityResult`

## Validation Helpers

- `validate_against_nlsolve(sys; tol=1e-4)`
- `validate_against_jump(sys; tol=1e-3)`

Both compare solver outputs and report max mismatch metrics.

## If Extension Is Not Loaded

Calling extension-dependent APIs raises an explicit error message with required packages.
