# JuliaPowerCase Integration

## 1. HybridPowerSystem Path (Preferred)

Supported directly:

- `solve_power_flow(hps::HybridPowerSystem)` — hybrid AC/DC
- `solve_dc_power_flow(hps::HybridPowerSystem)` — pure DC only
- `solve_power_flow_adaptive(hps::HybridPowerSystem)`
- `solve_power_flow_islanded(hps::HybridPowerSystem)`
- distributed-slack overloads for `HybridPowerSystem`

Use conversion explicitly when needed:

```julia
sys = to_solver_system(hps)
res = solve_power_flow(sys)
```

## 2. HybridPowerCaseData Path (Adapter)

Adapter APIs:

- `to_hybrid_system(h::HybridPowerCaseData) -> HybridSystem`
- `solve_power_flow(h::HybridPowerCaseData; update=true)`
- `update_results!(h, result)`

Example:

```julia
using JuliaPowerCase
using HybridACDCPowerFlow

h = case_hybrid_5ac3dc()
res = solve_power_flow(h; update=true)
```

With `update=true`, AC `VM/VA` and DC `VDC` are written back into `h`.

## 3. DC Bus Types

DC buses use the `bus_type` field from `Bus{DC}`:
- `DC_V` (= `REF_BUS`): Voltage reference bus — V fixed, P solved
- `DC_P` (= `PQ_BUS`): Power-specified bus — P specified, V solved

When creating DC-only systems:
```julia
# Bus 1 as voltage reference
dc_bus_1 = Bus{DC}(..., bus_type=DC_V, vm_pu=1.0, ...)
# Other buses as power-specified
dc_bus_2 = Bus{DC}(..., bus_type=DC_P, pd_mw=50.0, ...)
```

## 4. Notes

- Bus IDs are remapped to contiguous indices during adapter conversion.
- Converter control modes are mapped to Symbols: `:pq`, `:vdc_q`, `:vdc_vac`.
- Returned voltages always use solver conventions (`Va` in radians).
- If no `DC_V` bus is set, DC bus 1 defaults to reference (backward compatible).
