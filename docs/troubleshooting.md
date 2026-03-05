# Troubleshooting

## 1. Non-Convergence in `solve_power_flow`

Symptoms:

- `converged=false`
- high `residual`

Checks:

- verify slack bus exists in each AC island
- verify generator/load balance
- verify converter setpoints are physically consistent
- retry with warm-start `init`

## 2. Islanded/Adaptive Unexpected Zero Voltages

Adaptive solver sets dead or infeasible islands to zero voltage by design.

Action:

- inspect `result.islands`
- run `print_island_summary(detect_islands(sys), sys)`

## 3. Extension API Errors

If feasibility calls error with missing dependencies, load extension packages first:

```julia
using JuMP, Ipopt, NLsolve
```

## 4. Unit Mismatch Problems

Common issue: mixing MW/MVar and p.u. in manual system construction.

Rules:

- bus and generator powers are MW/MVar in system structs
- solver internally converts by `baseMVA`
- result angles `Va` are radians

## 5. MATPOWER Import Surprises

The parser remaps non-sequential bus IDs and may auto-convert units for distribution-style files.

Action:

- inspect returned `MatpowerData.busid_to_idx`
- confirm converted branch/load magnitudes before solving
