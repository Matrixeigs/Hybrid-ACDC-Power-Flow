# JuliaPowerCase Integration

## 1. HybridPowerSystem Path (Preferred)

Supported directly:

- `solve_power_flow(hps::HybridPowerSystem)`
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

## 3. Notes

- Bus IDs are remapped to contiguous indices during adapter conversion.
- Converter control modes are mapped to Symbols: `:pq`, `:vdc_q`, `:vdc_vac`.
- Returned voltages always use solver conventions (`Va` in radians).
