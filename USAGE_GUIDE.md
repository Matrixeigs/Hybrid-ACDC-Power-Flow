# HybridACDCPowerFlow Module - Migration Guide

## Overview

The `HybridACDCPowerFlow` module provides a clean, standalone implementation of hybrid AC/DC power flow analysis with comprehensive test cases for all NSFC theoretical claims.

## Directory Structure

```
HybridACDCPowerFlow/
├── Project.toml          # Package dependencies
├── README.md             # Module documentation
├── src/
│   ├── HybridACDCPowerFlow.jl  # Main module file
│   ├── PowerSystem.jl          # Power flow solver (759 lines)
│   └── TestSystems.jl          # IEEE test cases (413 lines)
└── test/
    └── runtests.jl       # Comprehensive theorem validation
```

## Quick Start

### Option 1: Direct Include (Recommended for Scripts)

```julia
include("HybridACDCPowerFlow/src/HybridACDCPowerFlow.jl")
using .HybridACDCPowerFlow

# Build system and solve
sys = build_ieee14_acdc()
result = solve_power_flow(sys)

println("Converged: $(result.converged)")
println("Voltages: $(result.Vm)")
```

### Option 2: Package Mode (For Development)

```julia
using Pkg
Pkg.develop(path="HybridACDCPowerFlow")
using HybridACDCPowerFlow
```

### Option 3: Run Tests Directly

```bash
julia HybridACDCPowerFlow/test/runtests.jl
```

## Available Functions

### Test Systems
- `build_ieee14_acdc()` - IEEE 14-bus + 2 DC buses
- `build_ieee24_3area_acdc()` - IEEE 24-bus + 4 DC buses (3 areas)
- `build_ieee118_acdc()` - IEEE 118-bus + 6 DC buses

### Power Flow Analysis
- `solve_power_flow(sys)` - Newton-Raphson solver
- `build_admittance_matrix(sys)` - AC Ybus matrix
- `power_flow_residual(sys, Vm, Va, Vdc)` - Compute residuals

### Data Structures
- `HybridSystem` - Complete AC/DC system
- `ACBus`, `ACBranch` - AC network components
- `DCBus`, `DCBranch` - DC network components  
- `VSCConverter` - Voltage source converter with control modes

### Control Modes
- `PQ_MODE` - Fixed P and Q injection
- `VDC_Q` - DC voltage and Q control
- `VDC_VAC` - DC and AC voltage control

## Example: N-k Fault Analysis

```julia
include("HybridACDCPowerFlow/src/HybridACDCPowerFlow.jl")
using .HybridACDCPowerFlow

# Normal operation
sys = build_ieee14_acdc()
res_normal = solve_power_flow(sys)

# N-2 fault (remove branches 1 and 3)
sys_fault = build_ieee14_acdc()
sys_fault.ac_branches[1] = ACBranch(
    sys_fault.ac_branches[1].from,
    sys_fault.ac_branches[1].to,
    sys_fault.ac_branches[1].r,
    sys_fault.ac_branches[1].x,
    sys_fault.ac_branches[1].b,
    sys_fault.ac_branches[1].tap,
    false  # status = false (removed)
)
sys_fault.ac_branches[3] = ACBranch(
    sys_fault.ac_branches[3].from,
    sys_fault.ac_branches[3].to,
    sys_fault.ac_branches[3].r,
    sys_fault.ac_branches[3].x,
    sys_fault.ac_branches[3].b,
    sys_fault.ac_branches[3].tap,
    false
)

res_fault = solve_power_flow(sys_fault)

println("Normal: converged=$(res_normal.converged)")
println("Fault:  converged=$(res_fault.converged)")
```

## Example: Control Mode Switching

```julia
# PQ mode
sys_pq = build_ieee24_3area_acdc()
c = sys_pq.converters[1]
sys_pq.converters[1] = VSCConverter(
    c.id, c.ac_bus, c.dc_bus, PQ_MODE,
    0.3, 0.05,  # P=0.3, Q=0.05
    c.Vdc_set, c.Vac_set,
    c.Ploss_a, c.Ploss_b, c.Ploss_c, c.Smax, true
)
res_pq = solve_power_flow(sys_pq)

# VDC_VAC mode (voltage control)
sys_vac = build_ieee24_3area_acdc()
sys_vac.converters[1] = VSCConverter(
    c.id, c.ac_bus, c.dc_bus, VDC_VAC,
    0.0, 0.0,  # P, Q not used in this mode
    c.Vdc_set, 1.05,  # Vac_set = 1.05 p.u.
    c.Ploss_a, c.Ploss_b, c.Ploss_c, c.Smax, true
)
res_vac = solve_power_flow(sys_vac)

println("PQ mode:    V = $(res_pq.Vm[c.ac_bus])")
println("VDC_VAC:    V = $(res_vac.Vm[c.ac_bus])")
```

## Theorem Validation

Run comprehensive validation of all 5 NSFC theorems:

```bash
julia HybridACDCPowerFlow/test/runtests.jl
```

Output includes:
- ✅ Theorem 1: Topology generalization (ε = 0)
- ✅ Theorem 2: SO(2) equivariance (< 0.02%)  
- ✅ Theorem 3: Current injection universality
- ✅ Theorem 4: DC shortcut effect (14% diameter reduction)
- ✅ Theorem 5: Control robustness (4.36e13× improvement)

## Migration from PEGNN_PF

**Old code:**
```julia
include("PEGNN_PF/PowerSystem.jl")
include("PEGNN_PF/TestSystems.jl")
using .PowerSystem, .TestSystems
```

**New code:**
```julia
include("HybridACDCPowerFlow/src/HybridACDCPowerFlow.jl")
using .HybridACDCPowerFlow
```

All function names and APIs remain identical.

## Benefits

1. **Clean Module Structure**: Proper Julia package with Project.toml
2. **Comprehensive Tests**: All 5 theorems validated automatically
3. **Self-Contained**: No external dependencies except stdlib
4. **Well-Documented**: README + inline comments
5. **Easy to Share**: Single folder contains everything
6. **Portable**: Works with `include()` or `Pkg.develop()`

## Next Steps

Use this module for:
- ✅ NSFC proposal validation (all theorems tested)
- ✅ Journal paper numerical experiments
- ✅ PE-GNN training data generation
- ✅ Resilience analysis with N-k faults
- ✅ Multi-area coordination studies
