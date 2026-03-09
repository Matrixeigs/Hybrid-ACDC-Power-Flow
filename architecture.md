# Architecture

## Purpose

`HybridACDCPowerFlow` is a **pure application layer** providing deterministic hybrid AC/DC power-flow solvers that:

- Accept `HybridPowerSystem` from **JuliaPowerCase** (information/data layer)
- Perform Newton-Raphson based power flow computations
- Return results as immutable NamedTuples
- Support AC/DC + converter coupling, island detection, adaptive solving

**Key Design Principle**: **JuliaPowerCase** = data structures, **HybridACDCPowerFlow** = algorithms

## Module Structure

- `src/PowerSystem.jl`: core solver and matrix assembly (Newton-Raphson)
- `src/PowerSystemEnhanced.jl`: islanding, adaptive workflow, distributed slack
- `src/MatpowerParser.jl`: MATPOWER parser into `MatpowerData`
- `src/TestSystems.jl`: IEEE and MATPOWER-based system builders (returns `HybridPowerSystem`)
- `src/JuliaPowerCaseAdapter.jl`: Legacy `HybridPowerCaseData` conversion/update bridge
- `src/FeasibilityChecker.jl`: extension API stubs + `FeasibilityResult`
- `ext/FeasibilityExt.jl`: JuMP/Ipopt/NLsolve implementations

## Data Model

### Public Interface (from JuliaPowerCase)

All data structures are imported from **JuliaPowerCase v1.0**:

- `ACBus`, `ACBranch`, `DCBus`, `DCBranch`, `Generator`, `VSCConverter`
- `HybridPowerSystem`: main container with `ac::ACSystem`, `dc::DCSystem`, and connecting components
- `IslandInfo`: island detection results

### Internal Solver Container (Not Exported)

- `SolverData`: Internal workspace that wraps `HybridPowerSystem` with cached matrices (`Ybus`, `Gdc`) and aggregated generation (`Pg`, `Qg`). 
  - **Users never interact with `SolverData` directly**
  - Automatically created via `to_solver_data(::HybridPowerSystem)` when needed
  - Provides efficient matrix operations for Newton-Raphson solver

## Solver Data Flow

1. **Input**: User provides `HybridPowerSystem` (from JuliaPowerCase) or legacy `HybridPowerCaseData`
2. **Internal Conversion**: System internally converts to `SolverData` via `to_solver_data()`
3. **Matrix Assembly**: `rebuild_matrices!` refreshes aggregated generation and sparse matrices
4. **Newton-Raphson**: Solves AC and DC equations with converter coupling
5. **Output**: Results returned as immutable `NamedTuple` with voltage magnitudes, angles, convergence status

```
User Code:                   HybridACDCPowerFlow Internal:
----------                   -----------------------------
HybridPowerSystem  ──────►  to_solver_data()  ──────►  SolverData
     (from JPC)                                          (cached Ybus, Gdc)
                                                               │
                                                               ▼
                                                       Newton-Raphson Solver
                                                               │
                                                               ▼
                                                         NamedTuple Result
                                                    (Vm, Va, Vdc, converged, ...)
```

## Design Notes

- **Sparse Jacobian** pattern is pre-built once per solve via `SolverWorkspace`
- **Jacobian numeric values** are updated in-place each iteration for efficiency
- **PowerSystemEnhanced** functions are functional at system level where possible (e.g., local PV→PQ conversion within island sub-systems)
- **No mutation** of input `HybridPowerSystem` - all algorithms are functionally pure from user perspective
