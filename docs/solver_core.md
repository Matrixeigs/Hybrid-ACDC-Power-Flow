# Core Solver

## Scope

The core solver in `PowerSystem.jl` solves a coupled AC/DC system using Newton-Raphson.

## State Variables

For `nac` AC buses and `ndc` DC buses:

- `Va` for all non-slack AC buses
- `Vm` for PQ AC buses
- `Vdc` for all non-slack DC buses (`2:ndc`)

## Residual Blocks

- AC active-power mismatch on non-slack buses
- AC reactive-power mismatch on PQ buses
- DC power mismatch on non-slack DC buses

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
