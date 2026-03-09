# Core Solver

## Scope

The core solver in `PowerSystem.jl` solves coupled AC/DC systems using Newton-Raphson.
It supports:
- Hybrid AC/DC power flow with VSC converters
- Pure AC power flow (no DC components)
- Pure DC power flow (no AC components)

## Performance Optimizations (v0.6.0)

The solver includes several performance optimizations:

### Sin/Cos Caching
- Pre-computes `sin(Va)` and `cos(Va)` once per iteration
- Uses angle-difference identities to avoid repeated trig calls in Jacobian:
  - `sin(θi - θj) = sin(θi)·cos(θj) - cos(θi)·sin(θj)`
  - `cos(θi - θj) = cos(θi)·cos(θj) + sin(θi)·sin(θj)`

### SIMD Vectorization
- Uses `@turbo` from LoopVectorization.jl for sin/cos computation
- Enables CPU SIMD instructions for parallel floating-point operations

### StaticArrays for Converters
- Converter Jacobian blocks use `SVector{2}` for stack allocation
- Eliminates heap allocations in converter loss Jacobian computation

### Typical Performance
| System | Buses | Time | Memory |
|--------|-------|------|--------|
| IEEE14 AC/DC | 16 | 0.05 ms | 223 KB |
| case33bw | 33 | 0.10 ms | 399 KB |
| case69 | 69 | 0.25 ms | 876 KB |
| IEEE118 AC/DC | 120 | 0.36 ms | 1.2 MB |
| case300 | 300 | 1.27 ms | 2.9 MB |

## Bus Types

### AC Bus Types
- `PQ` (1): Load bus — P and Q specified
- `PV` (2): Generator bus — P and V specified  
- `SLACK` (3): Reference bus — V and θ specified

### DC Bus Types
- `DC_P` (1): Power-specified bus — P specified, V solved
- `DC_V` (3): Voltage reference bus — V fixed, P solved

DC bus types use the same underlying `BusType` enum as AC:
```julia
const DC_V = REF_BUS  # Voltage reference (like AC slack)
const DC_P = PQ_BUS   # Power specified (like AC PQ)
```

## State Variables

For `nac` AC buses and `ndc` DC buses:

- `Va` for all non-slack AC buses
- `Vm` for PQ AC buses
- `Vdc` for all DC_P buses (excludes DC_V reference)

## Residual Blocks

- AC active-power mismatch on non-slack buses
- AC reactive-power mismatch on PQ buses
- DC power mismatch on DC_P buses

Converter terms are injected into both AC and DC schedules through:

- `converter_ac_injection`
- `converter_dc_injection`

## Jacobian and Performance

- Sparsity pattern built once in `create_solver_workspace`
- Values updated in `build_jacobian_triplets!`
- COO-to-CSC index map enables in-place `nzval` updates
- Sparse LU (`UMFPACK`) is reused with `lu!` fallback rebuild

## Initialization

- Default from bus records (`vm_pu`, `va_deg`, DC `vm_pu`)
- Optional warm start:

```julia
solve_power_flow(sys; init=(Vm=Vm0, Va=Va0, Vdc=Vdc0))
```

## Convergence and Safeguards

- Infinity norm on mismatch vector
- Backtracking line search (up to 8 trial steps)
- Voltage lower bounds applied (`>= 0.05` p.u.)
- Early return on singular/non-finite solve step

## Important Unit Convention

- Inputs in data model: MW/MVar, angles in degrees, voltages in p.u.
- Internal solve: per-unit power and radians
- Returned `Va` is in radians

## Pure DC Power Flow

For systems with only DC buses and branches (no AC components), use `solve_dc_power_flow`:

```julia
result = solve_dc_power_flow(sys; max_iter=50, tol=1e-8)
```

### DC Power Flow Equations

The DC power balance at each bus:
```
Pdc_k = Vdc_k × Σⱼ G_kj × Vdc_j
```

Where `G_kj = 1/R_kj` is the DC line conductance.

### DC Bus Type Assignment

- One bus must be `DC_V` (voltage reference) — its voltage is fixed
- All other buses are `DC_P` — their power is specified, voltage solved

```julia
using HybridACDCPowerFlow

# Create DC-only system
dc_buses = [
    make_dc_bus(1, DC_V, 1.0, 0.0),    # Reference at 1.0 p.u.
    make_dc_bus(2, DC_P, 1.0, 50.0),   # 50 MW load
    make_dc_bus(3, DC_P, 1.0, 30.0)    # 30 MW load
]
# ... add DC branches ...

result = solve_dc_power_flow(sys)
# result.Vdc = [1.0, 0.992, 0.989]
```

### Backward Compatibility

If no `DC_V` bus is explicitly set in the bus_type field, the solver defaults to using DC bus 1 as the voltage reference. This maintains backward compatibility with existing systems.
