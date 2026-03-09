# API Reference

## Public Data Types (from JuliaPowerCase)

All data structures are imported from **JuliaPowerCase v1.0**:

- `HybridPowerSystem`: Main power system container
- `ACBus`, `ACBranch`: AC network components
- `DCBus`, `DCBranch`: DC network components  
- `Generator`: Power generation units
- `VSCConverter`: Voltage-source converters connecting AC/DC
- `IslandInfo`: Island detection results

## Solver Configuration Types

- `SolverWorkspace`: Internal workspace for Jacobian assembly (not typically used directly)
- `ReactiveLimit`: Reactive power limits for PV buses
- `PowerFlowOptions`: Solver configuration options
- `DistributedSlack`: Distributed slack configuration
- `FeasibilityResult`: Feasibility check results

**Note**: `SolverData` is an internal type and not exported. Users work with `HybridPowerSystem`.

## Core Solvers

### Basic Power Flow

```julia
solve_power_flow(sys::HybridPowerSystem; max_iter=50, tol=1e-8, init=nothing)
```

**Primary interface** - accepts `HybridPowerSystem` from JuliaPowerCase

```julia
solve_power_flow(h::HybridPowerCaseData; max_iter=50, tol=1e-8, update=true)
```

**Legacy interface** - for backward compatibility with older case data format

**Result format**:

```julia
(Vm::Vector{Float64},     # AC bus voltage magnitudes
 Va::Vector{Float64},     # AC bus voltage angles (radians)
 Vdc::Vector{Float64},    # DC bus voltages
 converged::Bool,         # Convergence status
 iterations::Int,         # Number of iterations
 residual::Float64)       # Final residual norm
```

## Enhanced Solvers

### Island Detection and Solving

```julia
detect_islands(hps::HybridPowerSystem) -> Vector{IslandInfo}
```

Detect electrical islands in the system

```julia
extract_island_subsystem(hps::HybridPowerSystem, island::IslandInfo; slack_bus_override=0)
```

Extract a subsystem for a specific island

```julia
solve_power_flow_islanded(hps::HybridPowerSystem; options=PowerFlowOptions())
```

Solve power flow treating each island independently

### Adaptive Power Flow

```julia
solve_power_flow_adaptive(hps::HybridPowerSystem; 
                          options=PowerFlowOptions(), 
                          Q_limits=Dict())
```

Adaptive solver with automatic PV→PQ conversion when reactive limits are violated

## Distributed Slack

```julia
create_participation_factors(sys; 
                             method=:capacity, 
                             participating_buses=Int[], 
                             droop_coeffs=Dict())
```

Create participation factors for distributed slack

```julia
solve_power_flow_distributed_slack(hps::HybridPowerSystem, 
                                   dist_slack::DistributedSlack; 
                                   max_iter=50, tol=1e-8, verbose=false)
```

Solve with distributed slack bus model

```julia
solve_power_flow_distributed_slack_full(hps::HybridPowerSystem, 
                                        dist_slack::DistributedSlack; 
                                        max_iter=50, tol=1e-8, 
                                        verbose=false, enforce_limits=true)
```

Solve with distributed slack and limit enforcement

## Reactive and Converter Controls

```julia
create_default_Q_limits(sys; Qmin_default=-0.5, Qmax_default=0.5)
check_reactive_limits(sys, Vm, Va, Vdc, Q_limits)
auto_select_swing_bus(sys, island)
auto_switch_converter_mode!(sys, Vm, Va, Vdc; 
                            V_low_threshold=0.95, 
                            V_high_threshold=1.05, 
                            S_limit_frac=0.95)
```

## Utilities

```julia
get_bus_voltages(result)           # Extract bus voltages from result
get_branch_flows(sys, result)       # Compute branch power flows
extract_graph_data(sys)             # Extract graph connectivity data
```

## Test-System Builders

All builders return `HybridPowerSystem` compatible with the public API:

```julia
build_ieee14_acdc()          # IEEE 14-bus with MTDC
build_ieee24_3area_acdc()    # IEEE 24-bus 3-area with MTDC
build_ieee118_acdc()         # IEEE 118-bus with MTDC
build_ac_only_version(sys)   # Remove DC components

build_case33bw_acdc()        # 33-bus distribution
build_case33mg_acdc()        # 33-bus microgrid
build_case69_acdc()          # 69-bus distribution
build_case300_acdc()         # 300-bus transmission
build_case2000_acdc()        # 2000-bus large scale
```

## Feasibility Extension APIs

```julia
check_power_flow_feasibility(sys; 
                             method=:nlsolve, 
                             verbose=false, 
                             max_iter=50, 
                             max_time=10.0, 
                             tol=1e-6)

check_power_flow_feasibility_nlsolve(sys; kwargs...)
check_power_flow_feasibility_jump(sys; kwargs...)
validate_against_nlsolve(sys; kwargs...)
validate_against_jump(sys; kwargs...)
```

## Internal Functions (Advanced Use)

These functions are exported but typically only used internally or for advanced debugging:

```julia
to_solver_data(hps::HybridPowerSystem) -> SolverData  # Convert to internal solver format
rebuild_matrices!(sys)                                 # Rebuild admittance matrices
build_admittance_matrix(sys)                          # Build Ybus matrix
create_solver_workspace(sys)                          # Create Jacobian workspace
power_flow_residual(sys, Vm, Va, Vdc)                # Compute residual vector
```
