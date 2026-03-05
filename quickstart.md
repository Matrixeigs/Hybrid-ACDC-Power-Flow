# Quickstart

## 1. Requirements

- Julia `1.9+`
- Local package path: `HybridACDCPowerFlow`
- For feasibility extension only: `JuMP`, `Ipopt`, `NLsolve`

## 2. Install / Activate

```julia
using Pkg
Pkg.activate("HybridACDCPowerFlow")
Pkg.instantiate()
```

Or develop from another environment:

```julia
using Pkg
Pkg.develop(path="/absolute/path/to/HybridACDCPowerFlow")
```

## 3. Solve a Built-In System

```julia
using HybridACDCPowerFlow

sys = build_ieee14_acdc()
result = solve_power_flow(sys; max_iter=50, tol=1e-8)

println("converged=$(result.converged), iter=$(result.iterations), residual=$(result.residual)")
```

Result fields:

- `Vm`: AC voltage magnitude (p.u.)
- `Va`: AC voltage angle (rad)
- `Vdc`: DC voltage magnitude (p.u.)
- `converged`, `iterations`, `residual`

## 4. Adaptive / Island-Aware Solve

```julia
options = PowerFlowOptions(enable_pv_pq_conversion=true,
                           enable_auto_swing_selection=true,
                           enable_converter_mode_switching=true)
Q_limits = create_default_Q_limits(sys; Qmin_default=-0.5, Qmax_default=0.5)

res_adaptive = solve_power_flow_adaptive(sys; options=options, Q_limits=Q_limits)
println(res_adaptive.converged)
```

## 5. HybridPowerSystem Input (JuliaPowerCase)

```julia
using JuliaPowerCase
using HybridACDCPowerFlow

hps = HybridPowerSystem(ac=ac_sys, dc=dc_sys, vsc_converters=vscs, base_mva=100.0)
res = solve_power_flow(hps)
```

## 6. Run Tests

```bash
julia --project=HybridACDCPowerFlow HybridACDCPowerFlow/test/runtests.jl
julia --project=HybridACDCPowerFlow HybridACDCPowerFlow/test/test_jpc_integration.jl
```
