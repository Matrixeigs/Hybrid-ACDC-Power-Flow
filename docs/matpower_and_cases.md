# MATPOWER and Built-In Cases

## MATPOWER Parser

API:

- `parse_matpower(filepath::String) -> MatpowerData`

`MatpowerData` contains:

- `baseMVA`
- `bus`, `gen`, `branch` matrices
- `busid_to_idx` and original `bus_ids`

Parser behaviors:

- remaps non-sequential bus IDs to contiguous indices
- handles tap ratio `0.0 -> 1.0`
- detects distribution-style unit blocks and converts branch/load values as needed

## Built-In Case Builders

Hand-crafted hybrid benchmark systems:

- `build_ieee14_acdc()`
- `build_ieee24_3area_acdc()`
- `build_ieee118_acdc()`

MATPOWER-backed systems with HVDC augmentation:

- `build_case33bw_acdc()`
- `build_case33mg_acdc()`
- `build_case69_acdc()`
- `build_case300_acdc()`
- `build_case2000_acdc()`

## AC-Only Projection

- `build_ac_only_version(sys)` removes DC and converter layers from a hybrid case for comparative studies.

## Typical Workflow

```julia
using HybridACDCPowerFlow

sys = build_case300_acdc()
res = solve_power_flow(sys)
println((res.converged, res.iterations, res.residual))
```
