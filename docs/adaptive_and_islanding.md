# Adaptive and Islanded Solvers

## Island Detection

`detect_islands(sys::HybridSystem)` builds a combined graph:

- AC branches as AC-AC edges
- DC branches as DC-DC edges
- converters as AC-DC edges

This allows converter-connected buses to stay in the same electrical island.

## Adaptive Solver

`solve_power_flow_adaptive` workflow:

1. detect islands
2. skip dead/infeasible islands early
3. auto-select swing bus if needed
4. extract each island as a local sub-system
5. solve local NR
6. optional local PV->PQ conversion loop
7. optional local converter mode switching loop
8. map local voltages back to global arrays

Return fields include:

- `Vm`, `Va`, `Vdc`
- `converged`
- `iterations` (sum across inner solves)
- `islands`

## Islanded Solver

`solve_power_flow_islanded` returns one entry per island:

- `island_id`, `ac_buses`, `dc_buses`
- local voltage vectors
- per-island `converged` and `residual`

It calls the adaptive solver and then computes per-island residual checks.

## PowerFlowOptions

```julia
PowerFlowOptions(
    max_iter=200,
    tol=1e-6,
    enable_pv_pq_conversion=true,
    enable_auto_swing_selection=true,
    enable_converter_mode_switching=true,
    verbose=false,
)
```
