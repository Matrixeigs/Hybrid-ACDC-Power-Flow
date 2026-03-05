# Architecture

## Purpose

`HybridACDCPowerFlow` provides deterministic hybrid AC/DC power-flow solvers with:

- direct use of JuliaPowerCase types
- AC/DC + converter coupling
- island-aware and adaptive solve paths
- optional feasibility extension

## Module Structure

- `src/PowerSystem.jl`: core solver and matrix assembly
- `src/PowerSystemEnhanced.jl`: islanding, adaptive workflow, distributed slack
- `src/MatpowerParser.jl`: MATPOWER parser into `MatpowerData`
- `src/TestSystems.jl`: IEEE and MATPOWER-based system builders
- `src/JuliaPowerCaseAdapter.jl`: `HybridPowerCaseData` conversion/update bridge
- `src/FeasibilityChecker.jl`: extension API stubs + `FeasibilityResult`
- `ext/FeasibilityExt.jl`: JuMP/Ipopt/NLsolve implementations

## Data Model

Core types are imported from JuliaPowerCase:

- `ACBus`, `ACBranch`, `DCBus`, `DCBranch`, `Generator`, `VSCConverter`
- `HybridPowerSystem`, `IslandInfo`

Local solver container:

- `HybridSystem`: wraps component vectors, cached `Ybus`, cached `Gdc`, and aggregated `Pg/Qg`

## Solver Data Flow

1. Input can be `HybridSystem`, `HybridPowerSystem`, or `HybridPowerCaseData`.
2. Inputs are converted (if needed) into `HybridSystem`.
3. `rebuild_matrices!` refreshes aggregated generation and sparse matrices.
4. Newton-Raphson solves AC and DC equations with converter coupling.
5. Results are returned as a `NamedTuple`; adapter path can write back to case data.

## Design Notes

- Sparse Jacobian pattern is pre-built once per solve via `SolverWorkspace`.
- Jacobian numeric values are updated in-place each iteration.
- `PowerSystemEnhanced` functions are functional at system level where possible (for example, local PV->PQ conversion within island sub-systems).
