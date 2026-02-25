# Quick Reference — HybridACDCPowerFlow.jl v0.4.0

## Installation

```julia
using Pkg
Pkg.develop(path="/path/to/HybridACDCPowerFlow")    # local
# Pkg.add(url="https://github.com/<user>/HybridACDCPowerFlow.jl.git")  # remote

using HybridACDCPowerFlow
```

---

## Test Systems

```julia
sys = build_ieee14_acdc()          # 14 AC + 2 DC buses
sys = build_ieee24_3area_acdc()    # 24 AC + 4 DC buses (3 areas)
sys = build_ieee118_acdc()         # 118 AC + 8 DC buses

sys = build_case33bw_acdc()        # 33-bus radial distribution
sys = build_case33mg_acdc()        # 33-bus microgrid
sys = build_case69_acdc()          # 69-bus radial distribution
sys = build_case300_acdc()         # 300-bus transmission
sys = build_case2000_acdc()        # 2000-bus (scalability)

sys_ac = build_ac_only_version(sys)  # strip DC/converters
```

---

## Basic Power Flow

```julia
result = solve_power_flow(sys)

result.converged     # Bool
result.Vm            # AC voltage magnitudes
result.Va            # AC voltage angles (rad)
result.Vdc           # DC voltages
result.iterations    # Newton-Raphson iterations
result.residual      # Final residual norm

Vm, Va, Vdc = get_bus_voltages(result)
Pij, Qij = get_branch_flows(sys, result)
```

---

## Optimized Solver (SolverWorkspace)

```julia
ws = create_solver_workspace(sys)
result = solve_power_flow(sys; workspace=ws)

# Reuse ws across many solves with the same topology
```

---

## Island Detection

```julia
islands = detect_islands(sys)
print_island_summary(islands, sys)

for island in islands
    println("Island $(island.id): AC=$(island.ac_buses), slack=$(island.ac_slack_bus)")
end

island_results = solve_power_flow_islanded(sys)
```

---

## PV → PQ Conversion

```julia
Q_limits = create_default_Q_limits(sys; Qmin_default=-0.5, Qmax_default=1.0)

# Or manual:
Q_limits = Dict(2 => ReactiveLimit(-0.4, 0.8), 3 => ReactiveLimit(-0.3, 0.6))

violations, Q_actual = check_reactive_limits(sys, result.Vm, result.Va, result.Vdc, Q_limits)
```

---

## Adaptive Power Flow (All-in-One)

```julia
options = PowerFlowOptions(
    max_iter=20,
    tol=1e-8,
    enable_pv_pq_conversion=true,
    enable_auto_swing_selection=true,
    enable_converter_mode_switching=true,
    verbose=false
)

Q_limits = create_default_Q_limits(sys)
result = solve_power_flow_adaptive(sys; options=options, Q_limits=Q_limits)
```

---

## Converter Mode Switching

```julia
switched = auto_switch_converter_mode!(
    sys, result.Vm, result.Va, result.Vdc;
    V_low_threshold=0.95,
    V_high_threshold=1.05,
    S_limit_frac=0.90,
    hysteresis=0.02
)
```

Control modes: `PQ_MODE`, `VDC_Q`, `VDC_VAC`

---

## Distributed Slack

```julia
dist_slack = create_participation_factors(sys; method=:capacity)  # or :droop, :equal

result = solve_power_flow_distributed_slack(sys, dist_slack; verbose=true)
result = solve_power_flow_distributed_slack_full(sys, dist_slack; enforce_limits=true)
```

Manual configuration:

```julia
dist_slack = DistributedSlack(
    [1, 2, 6],           # participating buses
    [0.4, 0.4, 0.2],    # participation factors
    1,                    # reference bus
    Dict(2 => 0.05)      # capacity limits
)
```

---

## Feasibility Checking (Extension)

```julia
using JuMP, Ipopt, NLsolve  # triggers FeasibilityExt loading

result = check_power_flow_feasibility(sys; method=:nlsolve, verbose=true)
# or method=:jump for JuMP+Ipopt

result.feasible              # Bool
result.status                # String
result.load_margin           # Float64 (p.u.)
```

---

## Data Structures

| Type | Fields |
|------|--------|
| `ACBus` | `id, type, Pd, Qd, Pg, Qg, Vm, Va, area` |
| `ACBranch` | `from, to, r, x, b, tap, status` |
| `DCBus` | `id, Vdc_set, Pdc` |
| `DCBranch` | `from, to, r, status` |
| `VSCConverter` | `id, ac_bus, dc_bus, mode, Pset, Qset, Vdc_set, Vac_set, Ploss_a/b/c, Smax, G_droop, status` |
| `HybridSystem` | `ac_buses, ac_branches, baseMVA, dc_buses, dc_branches, converters, ...` |
| `BusType` | `PQ`, `PV`, `SLACK` |
| `ConverterMode` | `PQ_MODE`, `VDC_Q`, `VDC_VAC` |

---

## Running Tests

```julia
using Pkg
Pkg.test("HybridACDCPowerFlow")
```

---

## Source Layout

```
src/HybridACDCPowerFlow.jl   # module entry, exports
src/PowerSystem.jl            # core solver, sparse Jacobian
src/PowerSystemEnhanced.jl    # islands, PV→PQ, distributed slack
src/MatpowerParser.jl         # MATPOWER .m parser
src/TestSystems.jl            # test case builders
src/FeasibilityChecker.jl     # stub interface for extension
ext/FeasibilityExt.jl         # optional JuMP/Ipopt/NLsolve
```
