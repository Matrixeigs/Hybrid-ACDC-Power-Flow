# API Reference

## Core Types

- `HybridSystem`
- `SolverWorkspace`
- `ReactiveLimit`
- `PowerFlowOptions`
- `DistributedSlack`
- `FeasibilityResult`

## Core Build / Conversion

- `HybridSystem(hps::HybridPowerSystem)`
- `to_solver_system(hps::HybridPowerSystem) -> HybridSystem`
- `to_hybrid_system(h::HybridPowerCaseData) -> HybridSystem`
- `update_results!(h::HybridPowerCaseData, result)`

## Core Solvers

- `solve_power_flow(sys::HybridSystem; max_iter=50, tol=1e-8, init=nothing)`
- `solve_power_flow(hps::HybridPowerSystem; kwargs...)`
- `solve_power_flow(h::HybridPowerCaseData; max_iter=50, tol=1e-8, update=true)`

Result format:

```julia
(Vm::Vector{Float64},
 Va::Vector{Float64},
 Vdc::Vector{Float64},
 converged::Bool,
 iterations::Int,
 residual::Float64)
```

## Enhanced Solvers

- `detect_islands(sys::HybridSystem) -> Vector{IslandInfo}`
- `detect_islands(hps::HybridPowerSystem) -> Vector{IslandInfo}`
- `extract_island_subsystem(sys::HybridSystem, island::IslandInfo; slack_bus_override=0)`
- `extract_island_subsystem(hps::HybridPowerSystem, island::IslandInfo; slack_bus_override=0)`
- `solve_power_flow_adaptive(sys::HybridSystem; options=PowerFlowOptions(), Q_limits=Dict())`
- `solve_power_flow_adaptive(hps::HybridPowerSystem; options=..., Q_limits=...)`
- `solve_power_flow_islanded(sys::HybridSystem; options=PowerFlowOptions())`
- `solve_power_flow_islanded(hps::HybridPowerSystem; options=...)`

## Reactive and Converter Controls

- `create_default_Q_limits(sys; Qmin_default=-0.5, Qmax_default=0.5)`
- `check_reactive_limits(sys, Vm, Va, Vdc, Q_limits)`
- `pv_to_pq_conversion!(sys, violations, Q_limits)`
- `auto_select_swing_bus(sys, island)`
- `auto_switch_converter_mode!(sys, Vm, Va, Vdc; V_low_threshold=0.95, V_high_threshold=1.05, S_limit_frac=0.95)`

## Distributed Slack

- `create_participation_factors(sys; method=:capacity, participating_buses=Int[], droop_coeffs=Dict())`
- `solve_power_flow_distributed_slack(sys, dist_slack; max_iter=50, tol=1e-8, verbose=false)`
- `solve_power_flow_distributed_slack_full(sys, dist_slack; max_iter=50, tol=1e-8, verbose=false, enforce_limits=true)`
- `solve_power_flow_distributed_slack(hps, dist_slack; ...)`
- `solve_power_flow_distributed_slack_full(hps, dist_slack; ...)`

## Utilities

- `build_admittance_matrix(sys)`
- `rebuild_matrices!(sys)`
- `create_solver_workspace(sys)`
- `compute_power_injections!(ws, Y)`
- `build_jacobian_triplets!(ws, sys)`
- `power_flow_residual(sys, Vm, Va, Vdc)`
- `get_bus_voltages(result)`
- `get_branch_flows(sys, result)`
- `remove_ac_branch(sys, idx)`
- `extract_graph_data(sys)`

## Test-System Builders

- `build_ieee14_acdc()`
- `build_ieee24_3area_acdc()`
- `build_ieee118_acdc()`
- `build_ac_only_version(sys)`
- `build_case33bw_acdc()`
- `build_case33mg_acdc()`
- `build_case69_acdc()`
- `build_case300_acdc()`
- `build_case2000_acdc()`

## Feasibility Extension APIs

- `check_power_flow_feasibility(sys; method=:nlsolve, verbose=false, max_iter=50, max_time=10.0, tol=1e-6)`
- `check_power_flow_feasibility_nlsolve(sys; kwargs...)`
- `check_power_flow_feasibility_jump(sys; kwargs...)`
- `validate_against_nlsolve(sys; kwargs...)`
- `validate_against_jump(sys; kwargs...)`
