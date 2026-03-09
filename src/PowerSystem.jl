"""
    PowerSystem.jl

Hybrid AC/DC power flow solver module for JuliaPowerCase data model.

# Architecture (v0.6.0 - Pure Application Layer)
This module is a pure *application* layer that operates on data structures 
from JuliaPowerCase (the *information* layer).

## Data Flow
```
JuliaPowerCase.HybridPowerSystem (data)
         ↓ (to_solver_system)
   SolverData (internal workspace)
         ↓ (solve_power_flow)
    Results stored back in HybridPowerSystem
```

## Key Design Principles
1. **No data type definitions** - All types imported from JuliaPowerCase
2. **Input**: `HybridPowerSystem` from JuliaPowerCase 
3. **Internal**: `SolverData` workspace (not exported)
4. **Output**: Results written back to input data structure

## Supported Features
- AC power flow with Newton-Raphson solver
- DC network with linear analysis  
- VSC converters with multiple control modes (PQ, VDC_Q, VDC_VAC)
- Sparse Jacobian with pre-allocated workspace
- Island detection and adaptive solving

v0.6.0: Refactored as pure application layer (data from JuliaPowerCase)
"""
module PowerSystem

using LinearAlgebra, SparseArrays
using LoopVectorization: @turbo
using StaticArrays: SVector, SMatrix, @SVector, @SMatrix
using JuliaPowerCase: Bus, Branch, Generator,
                      VSCConverter, PowerSystem as JPCPowerSystem,
                      HybridPowerSystem,
                      BusType, PQ_BUS, PV_BUS, REF_BUS, ISOLATED_BUS,
                      AC, DC,
                      ACBus, ACBranch, DCBus, DCBranch,
                      # Topology functions
                      detect_islands as jpc_detect_islands,
                      extract_island_subsystem as jpc_extract_island_subsystem,
                      IslandInfo

# Re-export JuliaPowerCase types for convenient access
export ACBus, ACBranch, DCBus, DCBranch, VSCConverter, Generator,
       HybridPowerSystem, IslandInfo,
       BusType, PQ, PV, SLACK, DC_V, DC_P, PQ_MODE, VDC_Q, VDC_VAC,
       # Main solver functions
       solve_power_flow, solve_dc_power_flow, power_flow_residual,
       get_bus_voltages, get_branch_flows,
       # Conversion (limited - internal use mostly)
       to_solver_data,
       # Utility functions that may be needed externally
       rebuild_matrices!, full_residual_simple,
       # Topology functions from JuliaPowerCase
       jpc_detect_islands, jpc_extract_island_subsystem,
       # Loss model types
       LossModelType, LINEAR_LOSS, CURRENT_BASED_LOSS

# NOTE: SolverData (formerly HybridSystem) and SolverWorkspace are INTERNAL and not exported
# Users should interact with HybridPowerSystem from JuliaPowerCase

# ═══════════════════════════════════════════════════════════════════════════════
#  BUS TYPE ALIASES (from JuliaPowerCase)
# ═══════════════════════════════════════════════════════════════════════════════

const PQ = PQ_BUS
const PV = PV_BUS
const SLACK = REF_BUS

# DC bus type aliases (using BusType enum for DC buses)
# - DC_V: Voltage-controlled DC bus (reference/slack)
# - DC_P: Power-controlled DC bus
"""DC voltage reference bus (uses REF_BUS from BusType)."""
const DC_V = REF_BUS

"""DC power-specified bus (uses PQ_BUS from BusType)."""
const DC_P = PQ_BUS

# ═══════════════════════════════════════════════════════════════════════════════
#  CONVERTER CONTROL MODE (Symbol constants matching JuliaPowerCase)
# ═══════════════════════════════════════════════════════════════════════════════

"""Converter control modes (Symbol constants matching VSCConverter.control_mode)."""
const PQ_MODE  = :pq
const VDC_Q    = :vdc_q
const VDC_VAC  = :vdc_vac

# ═══════════════════════════════════════════════════════════════════════════════
#  CONVERTER LOSS MODEL TYPES
# ═══════════════════════════════════════════════════════════════════════════════

"""
    LossModelType

Enumeration for converter loss model selection.
- `LINEAR_LOSS`: Ploss = (1-η)|P| — simple efficiency-based model
- `CURRENT_BASED_LOSS`: Ploss = a·I² + b·I + c where I = |P|/Vdc
  - a = loss_mw/baseMVA, b = loss_percent/100, c = 1-eta
"""
@enum LossModelType begin
    LINEAR_LOSS         # Ploss = (1-η)|P|
    CURRENT_BASED_LOSS  # Ploss = a·I² + b·I + c, I = |P|/Vdc
end

# ═══════════════════════════════════════════════════════════════════════════════
#  SOLVER DATA (Internal Workspace - uses JuliaPowerCase types)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    SolverData

Internal solver workspace for hybrid AC/DC power flow analysis.
NOT EXPORTED - users should work with `HybridPowerSystem` from JuliaPowerCase.

This struct wraps component data from JuliaPowerCase with solver-specific cached values:
- Component references (busesConverters from HybridPowerSystem)
- Computed matrices (Ybus, Gdc)
- Aggregated data (Pg, Qg per bus)

# Fields
- `ac_buses::Vector{ACBus}` — references to AC buses from JuliaPowerCase
- `ac_branches::Vector{ACBranch}` — references to AC branches from JuliaPowerCase
- `dc_buses::Vector{DCBus}` — references to DC buses from JuliaPowerCase
- `dc_branches::Vector{DCBranch}` — references to DC branches from JuliaPowerCase  
- `converters::Vector{VSCConverter}` — references to converters from JuliaPowerCase
- `generators::Vector{Generator}` — references to generators from JuliaPowerCase
- `baseMVA::Float64` — base power for per-unit conversion
- `Pg, Qg::Vector{Float64}` — aggregated generation per AC bus (p.u.)
- `Ybus::SparseMatrixCSC{ComplexF64}` — AC bus admittance matrix
- `Gdc::SparseMatrixCSC{Float64}` — DC bus conductance matrix
- `loss_model::LossModelType` — converter loss model selection

# Constructor
    SolverData(hps::HybridPowerSystem; loss_model=LINEAR_LOSS)

Create solver workspace from a JuliaPowerCase `HybridPowerSystem`.
"""
mutable struct SolverData
    # Component references (from JuliaPowerCase, not owned)
    ac_buses::Vector{ACBus}
    ac_branches::Vector{ACBranch}
    baseMVA::Float64
    dc_buses::Vector{DCBus}
    dc_branches::Vector{DCBranch}
    converters::Vector{VSCConverter}
    generators::Vector{Generator}
    # Cached solver data
    Pg::Vector{Float64}
    Qg::Vector{Float64}
    Ybus::Union{Nothing, SparseMatrixCSC{ComplexF64, Int}}
    Gdc::Union{Nothing, SparseMatrixCSC{Float64, Int}}
    loss_model::LossModelType
end

"""
    aggregate_generation!(sys::SolverData)

Aggregate generator active/reactive power into per-bus Pg/Qg vectors (p.u.).
"""
function aggregate_generation!(sys::SolverData)
    nac = length(sys.ac_buses)
    fill!(sys.Pg, 0.0)
    fill!(sys.Qg, 0.0)
    for gen in sys.generators
        gen.in_service || continue
        bus_idx = gen.bus
        if 1 <= bus_idx <= nac
            sys.Pg[bus_idx] += gen.pg_mw / sys.baseMVA
            sys.Qg[bus_idx] += gen.qg_mvar / sys.baseMVA
        end
    end
end

function SolverData(ac_buses, ac_branches, dc_buses, dc_branches, converters;
                    generators=Generator[], baseMVA=100.0, 
                    loss_model::LossModelType=LINEAR_LOSS)
    nac = length(ac_buses)
    Pg = zeros(nac)
    Qg = zeros(nac)
    sys = SolverData(ac_buses, ac_branches, baseMVA, dc_buses, dc_branches,
                     converters, generators, Pg, Qg, nothing, nothing, loss_model)
    aggregate_generation!(sys)
    sys.Ybus = build_admittance_matrix(sys)
    sys.Gdc = build_dc_conductance(sys)
    return sys
end

function rebuild_matrices!(sys::SolverData)
    aggregate_generation!(sys)
    sys.Ybus = build_admittance_matrix(sys)
    sys.Gdc = build_dc_conductance(sys)
    return sys
end

# ═══════════════════════════════════════════════════════════════════════════════
#  CONVERSION: HybridPowerSystem → SolverData (solver-ready)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    to_solver_data(hps::HybridPowerSystem; loss_model=LINEAR_LOSS) -> SolverData

Convert a JuliaPowerCase `HybridPowerSystem` to internal solver workspace
with pre-computed admittance matrices.

This is called automatically by `solve_power_flow(hps)`, so users typically
don't need to call this directly.

# Example
```julia
using JuliaPowerCase, HybridACDCPowerFlow
hps = HybridPowerSystem(...)  # Define your system using JuliaPowerCase types
result = solve_power_flow(hps)  # Automatically creates solver workspace
```
"""
function to_solver_data(hps::HybridPowerSystem; loss_model::LossModelType=LINEAR_LOSS)
    # Extract components from nested structure
    ac_buses = hps.ac.buses
    ac_branches = hps.ac.branches
    dc_buses = hps.dc.buses
    dc_branches = hps.dc.branches
    converters = hps.vsc_converters
    generators = hps.ac.generators
    baseMVA = hps.base_mva
    
    return SolverData(ac_buses, ac_branches, dc_buses, dc_branches, converters;
                      generators=generators, baseMVA=baseMVA, loss_model=loss_model)
end

"""
    SolverData(hps::HybridPowerSystem; loss_model=LINEAR_LOSS)

Construct a solver-ready SolverData from a JuliaPowerCase HybridPowerSystem.
Equivalent to `to_solver_data(hps)`.
"""
SolverData(hps::HybridPowerSystem; loss_model::LossModelType=LINEAR_LOSS) = 
    to_solver_data(hps; loss_model=loss_model)

# ═══════════════════════════════════════════════════════════════════════════════
#  ADMITTANCE AND CONDUCTANCE MATRICES (SPARSE)
# ═══════════════════════════════════════════════════════════════════════════════

function build_admittance_matrix(sys::SolverData)
    n = length(sys.ac_buses)
    # Pre-count nonzeros for efficient sparse construction
    nnz_est = n + 2 * count(br -> br.in_service, sys.ac_branches)
    I = Vector{Int}(undef, nnz_est)
    J = Vector{Int}(undef, nnz_est)
    V = Vector{ComplexF64}(undef, nnz_est)
    
    # Diagonal entries first
    diag_vals = zeros(ComplexF64, n)
    idx = 0
    
    for br in sys.ac_branches
        br.in_service || continue
        i, j = br.from_bus, br.to_bus
        ys = 1.0 / complex(br.r_pu, br.x_pu)
        yc = im * br.b_pu / 2.0
        tap = br.tap == 0.0 ? 1.0 : br.tap
        
        # Diagonal contributions
        diag_vals[i] += ys / tap^2 + yc
        diag_vals[j] += ys + yc
        
        # Off-diagonal entries
        off_diag = -ys / tap
        idx += 1; I[idx] = i; J[idx] = j; V[idx] = off_diag
        idx += 1; I[idx] = j; J[idx] = i; V[idx] = off_diag
    end
    
    # Add diagonal entries
    for i in 1:n
        idx += 1; I[idx] = i; J[idx] = i; V[idx] = diag_vals[i]
    end
    
    resize!(I, idx); resize!(J, idx); resize!(V, idx)
    return sparse(I, J, V, n, n)
end

function build_dc_conductance(sys::SolverData)
    n = length(sys.dc_buses)
    n == 0 && return spzeros(Float64, 0, 0)
    
    # Pre-count nonzeros
    nnz_est = n + 2 * count(br -> br.in_service, sys.dc_branches)
    I = Vector{Int}(undef, nnz_est)
    J = Vector{Int}(undef, nnz_est)
    V = Vector{Float64}(undef, nnz_est)
    
    diag_vals = zeros(n)
    idx = 0
    
    for br in sys.dc_branches
        br.in_service || continue
        i, j = br.from_bus, br.to_bus
        g = 1.0 / br.r_pu
        diag_vals[i] += g
        diag_vals[j] += g
        idx += 1; I[idx] = i; J[idx] = j; V[idx] = -g
        idx += 1; I[idx] = j; J[idx] = i; V[idx] = -g
    end
    
    for i in 1:n
        idx += 1; I[idx] = i; J[idx] = i; V[idx] = diag_vals[i]
    end
    
    resize!(I, idx); resize!(J, idx); resize!(V, idx)
    return sparse(I, J, V, n, n)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  SOLVER WORKSPACE (PRE-ALLOCATED BUFFERS)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    SolverWorkspace

Pre-allocated workspace for Newton-Raphson solver. Eliminates allocations
during the NR iteration loop.
"""
mutable struct SolverWorkspace
    # Dimensions
    nac::Int
    ndc::Int
    nf::Int                    # total equations
    np::Int                    # P equations (non-slack)
    nq::Int                    # Q equations (PQ only)
    ndc_eq::Int                # DC equations (non-slack DC)
    
    # AC bus index sets (computed once)
    pq_idx::Vector{Int}
    pv_idx::Vector{Int}
    slack_idx::Vector{Int}
    non_slack::Vector{Int}
    vac_control_buses::Vector{Int}
    ac_isolated::Set{Int}
    
    # DC bus index sets
    dc_slack_idx::Int          # DC reference bus (0 if no DC or no ref)
    dc_p_idx::Vector{Int}      # DC P-type buses (power specified)
    dc_non_slack::Vector{Int}  # DC buses with equations (all except ref)
    dc_isolated::Set{Int}      # Isolated DC buses
    
    # AC index mappings for Jacobian
    bus_to_p_row::Vector{Int}  # bus → P equation row (0 if slack)
    bus_to_q_row::Vector{Int}  # bus → Q equation row (0 if not PQ)
    bus_to_va_col::Vector{Int} # bus → Va column (0 if slack/isolated)
    bus_to_vm_col::Vector{Int} # bus → Vm column (0 if not PQ)
    
    # DC index mappings for Jacobian
    dc_bus_to_pdc_row::Vector{Int}  # DC bus → Pdc equation row (0 if ref)
    dc_bus_to_vdc_col::Vector{Int}  # DC bus → Vdc column (0 if ref)
    
    # State vectors (reused)
    Vm::Vector{Float64}
    Va::Vector{Float64}
    Vdc::Vector{Float64}
    Vm_old::Vector{Float64}
    Va_old::Vector{Float64}
    Vdc_old::Vector{Float64}
    
    # Power injection buffers
    Pcalc::Vector{Float64}
    Qcalc::Vector{Float64}
    Psch::Vector{Float64}
    Qsch::Vector{Float64}
    V_complex::Vector{ComplexF64}
    I_complex::Vector{ComplexF64}
    S_complex::Vector{ComplexF64}
    
    # DC power buffers
    Pdc_calc::Vector{Float64}
    Pdc_sch::Vector{Float64}
    
    # Residual and update vectors
    F::Vector{Float64}
    dx::Vector{Float64}
    
    # Sparse Jacobian (COO triplets for value updates)
    J_I::Vector{Int}           # row indices (fixed)
    J_J::Vector{Int}           # col indices (fixed)
    J_V::Vector{Float64}       # values (updated each iter)
    J_nnz::Int                 # number of nonzeros
    J_sparse::SparseMatrixCSC{Float64, Int}  # assembled sparse matrix
    
    # LU factorization (reusable symbolic structure)
    lu_factor::Union{Nothing, SparseArrays.UMFPACK.UmfpackLU{Float64, Int}}
    
    # G and B sparse matrices (views into Ybus)
    G_sparse::SparseMatrixCSC{Float64, Int}
    B_sparse::SparseMatrixCSC{Float64, Int}
    
    # Pre-computed (row,col) -> nzval index map (avoids Dict allocation in NR loop)
    entry_map::Dict{Tuple{Int,Int}, Int}
    
    # Cached Ybus nonzero indices (avoids findnz allocation in NR loop)
    Ybus_rows::Vector{Int}
    Ybus_cols::Vector{Int}
    
    # COO to CSC nzval index mapping (for in-place Jacobian updates)
    coo_to_csc::Vector{Int}
    
    # Sin/cos cache for AC bus angles (avoids repeated trig calls in Jacobian)
    sin_Va::Vector{Float64}    # sin(Va[i]) for each AC bus
    cos_Va::Vector{Float64}    # cos(Va[i]) for each AC bus
end

"""
    create_solver_workspace(sys::SolverData) -> SolverWorkspace

Create pre-allocated workspace for Newton-Raphson solver.
"""
function create_solver_workspace(sys::SolverData)
    nac = length(sys.ac_buses)
    ndc = length(sys.dc_buses)
    
    # Identify AC bus types
    pq_idx = Int[]
    pv_idx = Int[]
    slack_idx = Int[]
    for (i, b) in enumerate(sys.ac_buses)
        if b.bus_type == PQ
            push!(pq_idx, i)
        elseif b.bus_type == PV
            push!(pv_idx, i)
        else
            push!(slack_idx, i)
        end
    end
    
    # Identify DC bus types (DC_V = REF_BUS, DC_P = PQ_BUS)
    dc_slack_idx = 0  # DC reference bus (voltage controlled)
    dc_p_idx = Int[]  # DC P-type buses (power specified)
    for (i, b) in enumerate(sys.dc_buses)
        if b.bus_type == DC_V  # REF_BUS
            if dc_slack_idx == 0
                dc_slack_idx = i
            end
            # Additional DC_V buses treated as P-type (only one reference allowed)
        else  # DC_P (PQ_BUS) or other
            push!(dc_p_idx, i)
        end
    end
    # If no explicit DC reference and we have DC buses, use bus 1 as default
    if dc_slack_idx == 0 && ndc > 0
        dc_slack_idx = 1
    end
    
    # VDC_VAC converters (treat as PV-like)
    vac_control_buses = Int[]
    for conv in sys.converters
        if conv.in_service && conv.control_mode == VDC_VAC
            push!(vac_control_buses, conv.bus_ac)
        end
    end
    
    pq_idx_modified = setdiff(pq_idx, vac_control_buses)
    pv_idx_modified = union(pv_idx, vac_control_buses)
    
    # Detect AC-isolated buses
    ac_isolated = Set{Int}()
    if nac > 0
        for i in 1:nac
            has_branch = any(br.in_service && (br.from_bus == i || br.to_bus == i) for br in sys.ac_branches)
            if !has_branch
                push!(ac_isolated, i)
            end
        end
    end
    
    # Detect DC-isolated buses
    dc_isolated = Set{Int}()
    if ndc > 0
        for i in 1:ndc
            has_branch = any(br.in_service && (br.from_bus == i || br.to_bus == i) for br in sys.dc_branches)
            has_conv = any(c.in_service && c.bus_dc == i for c in sys.converters)
            if !has_branch && !has_conv
                push!(dc_isolated, i)
            end
        end
    end
    
    pq_idx_final = sort(setdiff(pq_idx_modified, ac_isolated))
    pv_idx_final = sort(setdiff(pv_idx_modified, ac_isolated))
    non_slack = sort(union(pq_idx_final, pv_idx_final))
    
    # DC non-slack buses (all except reference and isolated)
    dc_non_slack = Int[]
    for i in 1:ndc
        if i != dc_slack_idx && i ∉ dc_isolated
            push!(dc_non_slack, i)
        end
    end
    sort!(dc_non_slack)
    
    np = length(non_slack)
    nq = length(pq_idx_final)
    ndc_eq = length(dc_non_slack)
    nf = np + nq + ndc_eq
    
    # Build AC index mappings
    bus_to_p_row = zeros(Int, max(nac, 1))
    bus_to_q_row = zeros(Int, max(nac, 1))
    bus_to_va_col = zeros(Int, max(nac, 1))
    bus_to_vm_col = zeros(Int, max(nac, 1))
    
    for (k, i) in enumerate(non_slack)
        bus_to_p_row[i] = k
        bus_to_va_col[i] = k
    end
    for (k, i) in enumerate(pq_idx_final)
        bus_to_q_row[i] = np + k
        bus_to_vm_col[i] = np + k
    end
    
    # Build DC index mappings
    dc_bus_to_pdc_row = zeros(Int, max(ndc, 1))
    dc_bus_to_vdc_col = zeros(Int, max(ndc, 1))
    for (k, i) in enumerate(dc_non_slack)
        dc_bus_to_pdc_row[i] = np + nq + k
        dc_bus_to_vdc_col[i] = np + nq + k
    end
    
    # Pre-allocate state vectors
    Vm = ones(nac)
    Va = zeros(nac)
    Vdc = ndc > 0 ? ones(ndc) : Float64[]
    Vm_old = similar(Vm)
    Va_old = similar(Va)
    Vdc_old = similar(Vdc)
    
    # Sin/cos cache for angle differences (populated during Jacobian build)
    sin_Va = zeros(nac)
    cos_Va = zeros(nac)
    
    # Power buffers
    Pcalc = zeros(nac)
    Qcalc = zeros(nac)
    Psch = zeros(nac)
    Qsch = zeros(nac)
    V_complex = zeros(ComplexF64, nac)
    I_complex = zeros(ComplexF64, nac)
    S_complex = zeros(ComplexF64, nac)
    
    # DC buffers
    Pdc_calc = zeros(ndc)
    Pdc_sch = zeros(ndc)
    
    # Residual and update
    F = zeros(nf)
    dx = zeros(nf)
    
    # Build Jacobian sparsity pattern (passing DC bus type info)
    J_I, J_J, J_V, J_nnz = build_jacobian_sparsity(sys, non_slack, pq_idx_final, 
                                                    slack_idx, ac_isolated, nac, ndc, np, nq, ndc_eq,
                                                    dc_slack_idx, dc_non_slack, dc_isolated)
    J_sparse = sparse(J_I, J_J, J_V, nf, nf)
    
    # Extract G and B from Ybus
    Y = sys.Ybus === nothing ? build_admittance_matrix(sys) : sys.Ybus
    G_sparse = real.(Y)
    B_sparse = imag.(Y)
    
    # Pre-compute entry_map: (row, col) -> index in J_V
    entry_map = Dict{Tuple{Int,Int}, Int}()
    sizehint!(entry_map, J_nnz)
    for k in 1:J_nnz
        entry_map[(J_I[k], J_J[k])] = k
    end
    
    # Pre-compute Ybus nonzero indices (avoid findnz allocation in NR loop)
    Ybus_rows, Ybus_cols, _ = findnz(Y)
    
    # Build COO to CSC nzval index mapping for in-place Jacobian updates
    # This avoids calling sparse() on every iteration
    coo_to_csc = Vector{Int}(undef, J_nnz)
    for k in 1:J_nnz
        row, col = J_I[k], J_J[k]
        found = false
        # Find position in CSC: look in column col's range for row
        for idx in J_sparse.colptr[col]:(J_sparse.colptr[col+1]-1)
            if J_sparse.rowval[idx] == row
                coo_to_csc[k] = idx
                found = true
                break
            end
        end
        @assert found "Sparsity mismatch: COO entry ($row, $col) not found in CSC matrix"
    end
    
    return SolverWorkspace(
        nac, ndc, nf, np, nq, ndc_eq,
        pq_idx_final, pv_idx_final, slack_idx, non_slack, vac_control_buses, ac_isolated,
        dc_slack_idx, dc_p_idx, dc_non_slack, dc_isolated,
        bus_to_p_row, bus_to_q_row, bus_to_va_col, bus_to_vm_col,
        dc_bus_to_pdc_row, dc_bus_to_vdc_col,
        Vm, Va, Vdc, Vm_old, Va_old, Vdc_old,
        Pcalc, Qcalc, Psch, Qsch, V_complex, I_complex, S_complex,
        Pdc_calc, Pdc_sch,
        F, dx,
        J_I, J_J, J_V, J_nnz, J_sparse,
        nothing,
        G_sparse, B_sparse,
        entry_map,
        Ybus_rows, Ybus_cols,
        coo_to_csc,
        sin_Va, cos_Va
    )
end

"""
Build the sparsity pattern for the Jacobian. Returns (I, J, V, nnz).
The pattern is determined by Ybus sparsity plus DC coupling.
"""
function build_jacobian_sparsity(sys, non_slack, pq_idx, slack_idx, ac_isolated, 
                                  nac, ndc, np, nq, ndc_eq,
                                  dc_slack_idx=1, dc_non_slack=Int[], dc_isolated=Set{Int}())
    Y = sys.Ybus
    nf = np + nq + ndc_eq
    
    # Handle pure DC case (no AC)
    if nac == 0 || Y === nothing
        # For pure DC, only DC Jacobian
        if ndc_eq == 0
            return Int[], Int[], Float64[], 0
        end
    end
    
    # Count nonzeros: each Y[i,j] ≠ 0 contributes to J
    # For each (i,j) with Y[i,j] ≠ 0:
    #   - If i ∉ slack and j ∉ slack|isolated: H[i,j] entry
    #   - If i ∉ slack and j ∈ pq: N[i,j] entry
    #   - If i ∈ pq and j ∉ slack|isolated: J[i,j] entry
    #   - If i ∈ pq and j ∈ pq: L[i,j] entry
    
    non_slack_set = Set(non_slack)
    pq_set = Set(pq_idx)
    slack_set = Set(slack_idx)
    
    # Estimate nnz (upper bound)
    Y_nnz = Y === nothing ? 0 : nnz(Y)
    Gdc_nnz = ndc > 0 && sys.Gdc !== nothing ? nnz(sys.Gdc) : 0
    nnz_est = 4 * Y_nnz + 4 * Gdc_nnz + 10 * length(sys.converters) + 1
    
    I = Vector{Int}(undef, nnz_est)
    J = Vector{Int}(undef, nnz_est)
    V = Vector{Float64}(undef, nnz_est)
    idx = 0
    
    # Build column maps for AC
    col_va = Dict(non_slack[k] => k for k in 1:np)
    col_vm = Dict(pq_idx[k] => np + k for k in 1:nq)
    
    # Build column maps for DC
    dc_non_slack_set = Set(dc_non_slack)
    col_vdc = Dict{Int, Int}()
    for (k, dc_bus) in enumerate(dc_non_slack)
        col_vdc[dc_bus] = np + nq + k
    end
    
    # AC Jacobian blocks (H, N, J, L)
    if Y !== nothing && nac > 0
        rows_Y, cols_Y, _ = findnz(Y)
        for k in 1:length(rows_Y)
            i, j = rows_Y[k], cols_Y[k]
            
            # P equation row for bus i
            if i ∉ slack_set && i ∉ ac_isolated && haskey(col_va, i)
                row_p = col_va[i]
                
                # ∂P_i/∂θ_j (H block)
                if haskey(col_va, j)
                    idx += 1; I[idx] = row_p; J[idx] = col_va[j]; V[idx] = 0.0
                end
                # ∂P_i/∂V_j (N block)
                if haskey(col_vm, j)
                    idx += 1; I[idx] = row_p; J[idx] = col_vm[j]; V[idx] = 0.0
                end
            end
            
            # Q equation row for bus i (only if i ∈ PQ)
            if i ∈ pq_set && haskey(col_vm, i)
                row_q = col_vm[i]  # Q row = np + pq_position
                
                # ∂Q_i/∂θ_j (J block)
                if haskey(col_va, j)
                    idx += 1; I[idx] = row_q; J[idx] = col_va[j]; V[idx] = 0.0
                end
                # ∂Q_i/∂V_j (L block)
                if haskey(col_vm, j)
                    idx += 1; I[idx] = row_q; J[idx] = col_vm[j]; V[idx] = 0.0
                end
            end
        end
    end
    
    # DC Jacobian block (using actual DC bus types)
    if ndc_eq > 0 && sys.Gdc !== nothing
        Gdc = sys.Gdc
        rows_G, cols_G, _ = findnz(Gdc)
        for k in 1:length(rows_G)
            ki, li = rows_G[k], cols_G[k]
            # Skip DC slack bus (reference) and isolated DC buses
            ki == dc_slack_idx && continue
            li == dc_slack_idx && continue
            ki ∈ dc_isolated && continue
            li ∈ dc_isolated && continue
            
            haskey(col_vdc, ki) || continue
            haskey(col_vdc, li) || continue
            row = col_vdc[ki]
            col = col_vdc[li]
            idx += 1; I[idx] = row; J[idx] = col; V[idx] = 0.0
        end
        
        # Converter DC contributions (droop modes)
        for conv in sys.converters
            conv.in_service || continue
            (conv.control_mode == VDC_Q || conv.control_mode == VDC_VAC) || continue
            k = conv.bus_dc
            k == dc_slack_idx && continue
            k ∈ dc_isolated && continue
            haskey(col_vdc, k) || continue
            row = col_vdc[k]
            col = col_vdc[k]
            # This entry may already exist from Gdc, but adding duplicate is okay
            idx += 1; I[idx] = row; J[idx] = col; V[idx] = 0.0
        end
        
        # AC-DC cross-coupling for droop converters
        for conv in sys.converters
            conv.in_service || continue
            (conv.control_mode == VDC_Q || conv.control_mode == VDC_VAC) || continue
            k = conv.bus_dc
            k == dc_slack_idx && continue
            k ∈ dc_isolated && continue
            haskey(col_va, conv.bus_ac) || continue
            haskey(col_vdc, k) || continue
            row_p = col_va[conv.bus_ac]
            col_dc = col_vdc[k]
            idx += 1; I[idx] = row_p; J[idx] = col_dc; V[idx] = 0.0
        end
    end
    
    resize!(I, idx); resize!(J, idx); resize!(V, idx)
    return I, J, V, idx
end

# ═══════════════════════════════════════════════════════════════════════════════
#  IN-PLACE POWER INJECTION COMPUTATION
# ═══════════════════════════════════════════════════════════════════════════════

"""
Compute AC power injections in-place into workspace buffers.
Also caches sin(Va) and cos(Va) for use in Jacobian building.
"""
function compute_power_injections!(ws::SolverWorkspace, Y::SparseMatrixCSC{ComplexF64, Int})
    nac = ws.nac
    Vm = ws.Vm
    Va = ws.Va
    sin_Va = ws.sin_Va
    cos_Va = ws.cos_Va
    
    # Cache sin/cos values with @turbo for SIMD acceleration
    @turbo for i in 1:nac
        sin_Va[i] = sin(Va[i])
        cos_Va[i] = cos(Va[i])
    end
    
    # Build complex voltage vector (requires sequential due to complex ops)
    @inbounds for i in 1:nac
        ws.V_complex[i] = complex(Vm[i] * cos_Va[i], Vm[i] * sin_Va[i])
    end
    
    # Sparse matrix-vector product (uses BLAS)
    mul!(ws.I_complex, Y, ws.V_complex)
    
    # Compute power injections
    @inbounds for i in 1:nac
        ws.S_complex[i] = ws.V_complex[i] * conj(ws.I_complex[i])
        ws.Pcalc[i] = real(ws.S_complex[i])
        ws.Qcalc[i] = imag(ws.S_complex[i])
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
#  CONVERTER LOSS MODELS (using JuliaPowerCase VSCConverter fields)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    converter_loss(conv, P, Vdc, baseMVA, loss_model) -> Float64

Compute converter power loss based on selected model.

# Models
- `LINEAR_LOSS`: Ploss = (1-η)|P|
- `CURRENT_BASED_LOSS`: Ploss = a·I² + b·I + c where I = |P|/Vdc
  - a = loss_mw/baseMVA, b = loss_percent/100, c = 1-eta

# Arguments
- `conv`: VSCConverter from JuliaPowerCase
- `P`: Power flow through converter (p.u.)
- `Vdc`: DC voltage at converter bus (p.u.)
- `baseMVA`: System base MVA
- `loss_model`: LossModelType enum
"""
@inline function converter_loss(conv::VSCConverter, P::Float64, Vdc_val::Float64, 
                                baseMVA::Float64, loss_model::LossModelType)
    if loss_model == LINEAR_LOSS
        # Ploss = (1 - η) × |P|
        return (1.0 - conv.eta) * abs(P)
    else  # CURRENT_BASED_LOSS
        # Ploss = a·I² + b·I + c, where I = |P|/Vdc
        a = conv.loss_mw / baseMVA    # quadratic coeff
        b = conv.loss_percent / 100.0  # linear coeff
        c = 1.0 - conv.eta             # constant term (repurposed)
        I = abs(P) / max(Vdc_val, 0.1) # current magnitude, with floor to avoid div/0
        return a * I^2 + b * I + c
    end
end

"""
    converter_loss_jacobian(conv, P, Vdc, baseMVA, loss_model) -> (dPloss_dP, dPloss_dVdc)

Compute Jacobian entries for converter loss: ∂Ploss/∂P and ∂Ploss/∂Vdc.

# Models
- `LINEAR_LOSS`: dPloss/dP = (1-η)·sign(P), dPloss/dVdc = 0
- `CURRENT_BASED_LOSS`: 
  - dPloss/dP = (2a·|P|/Vdc² + b/Vdc)·sign(P)
  - dPloss/dVdc = -2a·P²/Vdc³ - b·|P|/Vdc²

Returns an `SVector{2}` for stack allocation: [dPloss_dP, dPloss_dVdc].
"""
@inline function converter_loss_jacobian(conv::VSCConverter, P::Float64, Vdc_val::Float64,
                                         baseMVA::Float64, loss_model::LossModelType)
    if loss_model == LINEAR_LOSS
        # dPloss/dP = (1-η)·sign(P)
        dPloss_dP = (1.0 - conv.eta) * sign(P + 1e-30)
        return @SVector [dPloss_dP, 0.0]
    else  # CURRENT_BASED_LOSS
        a = conv.loss_mw / baseMVA
        b = conv.loss_percent / 100.0
        Vdc_safe = max(Vdc_val, 0.1)
        absP = abs(P)
        sgn = sign(P + 1e-30)
        
        # I = |P|/Vdc, Ploss = a·I² + b·I + c
        # dPloss/dP = d/dP[a·P²/Vdc² + b·|P|/Vdc + c]
        #           = 2a·P/Vdc² + b·sign(P)/Vdc
        dPloss_dP = 2.0 * a * P / Vdc_safe^2 + b * sgn / Vdc_safe
        
        # dPloss/dVdc = d/dVdc[a·P²/Vdc² + b·|P|/Vdc + c]
        #             = -2a·P²/Vdc³ - b·|P|/Vdc²
        dPloss_dVdc = -2.0 * a * P^2 / Vdc_safe^3 - b * absP / Vdc_safe^2
        
        return @SVector [dPloss_dP, dPloss_dVdc]
    end
end

"""
    compute_converter_jacobian_block(conv, Vdc, baseMVA, loss_model) -> SVector{2}

Compute converter Jacobian contributions as static vector for stack allocation.
Returns [∂Pdc_diag, ∂Pac_∂Vdc] where:
- ∂Pdc_diag: Contribution to DC diagonal (for VDC modes)
- ∂Pac_∂Vdc: AC-DC cross-coupling term

This uses StaticArrays for zero-allocation Jacobian block computation.
"""
@inline function compute_converter_jacobian_block(conv::VSCConverter, Vdc::Vector{Float64},
                                                  baseMVA::Float64, loss_model::LossModelType)
    Vdc_val = Vdc[conv.bus_dc]
    k_vdc = conv.k_vdc
    
    # DC diagonal contribution: 2·k_vdc·Vdc
    dPdc_diag = 2.0 * k_vdc * Vdc_val
    
    # AC-DC coupling through loss model
    P_transfer = k_vdc * (Vdc_val^2 - conv.v_dc_set_pu^2)
    loss_jac = converter_loss_jacobian(conv, P_transfer, Vdc_val, baseMVA, loss_model)
    dPloss_dP = loss_jac[1]
    dPloss_dVdc = loss_jac[2]
    
    # ∂Pac/∂Vdc = ∂Pdc/∂Vdc × (1 - ∂Ploss/∂P) - ∂Ploss/∂Vdc
    dPac_dVdc = dPdc_diag * (1.0 - dPloss_dP) - dPloss_dVdc
    
    return @SVector [dPdc_diag, dPac_dVdc]
end

# Legacy compatibility: single-argument version defaults to LINEAR_LOSS with Vdc=1.0
@inline function converter_loss(conv::VSCConverter, P::Float64, baseMVA::Float64)
    return converter_loss(conv, P, 1.0, baseMVA, LINEAR_LOSS)
end

@inline function conv_dc_power(conv::VSCConverter, Vdc::Vector{Float64}, baseMVA::Float64,
                               loss_model::LossModelType=LINEAR_LOSS)
    if conv.control_mode == PQ_MODE
        Pset_pu = conv.p_set_mw / baseMVA
        Vdc_val = Vdc[conv.bus_dc]
        Ploss = converter_loss(conv, Pset_pu, Vdc_val, baseMVA, loss_model)
        return -(Pset_pu + Ploss)
    else  # VDC_Q or VDC_VAC
        return conv.k_vdc * (Vdc[conv.bus_dc]^2 - conv.v_dc_set_pu^2)
    end
end

function converter_ac_injection(conv::VSCConverter, Vm, Va, Vdc, baseMVA::Float64,
                                loss_model::LossModelType=LINEAR_LOSS)
    Pset_pu = conv.p_set_mw / baseMVA
    Qset_pu = conv.q_set_mvar / baseMVA
    Vdc_val = Vdc[conv.bus_dc]
    if conv.control_mode == PQ_MODE
        return Pset_pu, Qset_pu
    elseif conv.control_mode == VDC_Q
        P_transfer = conv_dc_power(conv, Vdc, baseMVA, loss_model)
        Ploss = converter_loss(conv, P_transfer, Vdc_val, baseMVA, loss_model)
        Pac = P_transfer - Ploss
        return Pac, Qset_pu
    elseif conv.control_mode == VDC_VAC
        P_transfer = conv_dc_power(conv, Vdc, baseMVA, loss_model)
        Ploss = converter_loss(conv, P_transfer, Vdc_val, baseMVA, loss_model)
        Pac = P_transfer - Ploss
        return Pac, 0.0
    end
    return 0.0, 0.0
end

function converter_dc_injection(conv::VSCConverter, Vm, Va, Vdc, baseMVA::Float64,
                                loss_model::LossModelType=LINEAR_LOSS)
    if conv.control_mode == PQ_MODE
        Pset_pu = conv.p_set_mw / baseMVA
        Vdc_val = Vdc[conv.bus_dc]
        Ploss = converter_loss(conv, Pset_pu, Vdc_val, baseMVA, loss_model)
        return -(Pset_pu + Ploss)
    else
        return -conv_dc_power(conv, Vdc, baseMVA, loss_model)
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
#  JACOBIAN BUILDER (SHARED IMPLEMENTATION)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    build_jacobian_triplets!(ws, sys)

Update Jacobian values in COO format. Sparsity pattern is fixed;
only values are updated. Uses @inbounds and cached sin/cos with
angle-difference identities for SIMD-friendly computation.
"""
function build_jacobian_triplets!(ws::SolverWorkspace, sys::SolverData)
    G = ws.G_sparse
    B = ws.B_sparse
    Vm = ws.Vm
    Va = ws.Va
    Vdc = ws.Vdc
    Pcalc = ws.Pcalc
    Qcalc = ws.Qcalc
    nac = ws.nac
    ndc = ws.ndc
    np = ws.np
    nq = ws.nq
    ndc_eq = ws.ndc_eq
    
    # Cached sin/cos values (computed in compute_power_injections!)
    sin_Va = ws.sin_Va
    cos_Va = ws.cos_Va
    
    # Reset Jacobian values
    fill!(ws.J_V, 0.0)
    
    # Use pre-computed entry map (no allocation)
    J_I = ws.J_I
    J_J = ws.J_J
    J_V = ws.J_V
    entry_map = ws.entry_map
    
    bus_to_p_row = ws.bus_to_p_row
    bus_to_q_row = ws.bus_to_q_row
    bus_to_va_col = ws.bus_to_va_col
    bus_to_vm_col = ws.bus_to_vm_col
    pq_idx = ws.pq_idx
    slack_idx = ws.slack_idx
    ac_isolated = ws.ac_isolated
    
    # Use cached Ybus nonzero indices (no allocation)
    rows = ws.Ybus_rows
    cols = ws.Ybus_cols
    
    @inbounds for k in 1:length(rows)
        i, j = rows[k], cols[k]
        
        i in slack_idx && continue
        i in ac_isolated && continue
        
        Gij = G[i, j]
        Bij = B[i, j]
        (Gij == 0.0 && Bij == 0.0) && continue
        
        # Use angle-difference identities with cached sin/cos:
        # sin(θi - θj) = sin(θi)cos(θj) - cos(θi)sin(θj)
        # cos(θi - θj) = cos(θi)cos(θj) + sin(θi)sin(θj)
        sin_i, cos_i = sin_Va[i], cos_Va[i]
        sin_j, cos_j = sin_Va[j], cos_Va[j]
        sinθ = sin_i * cos_j - cos_i * sin_j
        cosθ = cos_i * cos_j + sin_i * sin_j
        
        row_p = bus_to_p_row[i]
        row_q = bus_to_q_row[i]
        col_va_j = bus_to_va_col[j]
        col_vm_j = bus_to_vm_col[j]
        
        if i == j
            # Diagonal elements
            H_ii = -Qcalc[i] - Bij * Vm[i]^2
            N_ii = Pcalc[i] / Vm[i] + Gij * Vm[i]
            J_Qθ_ii = Pcalc[i] - Gij * Vm[i]^2
            L_ii = Qcalc[i] / Vm[i] - Bij * Vm[i]
            
            # P row
            if row_p > 0 && col_va_j > 0
                key = (row_p, col_va_j)
                haskey(entry_map, key) && (J_V[entry_map[key]] = H_ii)
            end
            if row_p > 0 && col_vm_j > 0
                key = (row_p, col_vm_j)
                haskey(entry_map, key) && (J_V[entry_map[key]] = N_ii)
            end
            
            # Q row
            if row_q > 0 && col_va_j > 0
                key = (row_q, col_va_j)
                haskey(entry_map, key) && (J_V[entry_map[key]] = J_Qθ_ii)
            end
            if row_q > 0 && col_vm_j > 0
                key = (row_q, col_vm_j)
                haskey(entry_map, key) && (J_V[entry_map[key]] = L_ii)
            end
        else
            # Off-diagonal elements
            H_ij = Vm[i] * Vm[j] * (Gij * sinθ - Bij * cosθ)
            N_ij = Vm[i] * (Gij * cosθ + Bij * sinθ)
            J_Qθ_ij = -Vm[i] * Vm[j] * (Gij * cosθ + Bij * sinθ)
            L_ij = Vm[i] * (Gij * sinθ - Bij * cosθ)
            
            # P row
            if row_p > 0 && col_va_j > 0
                key = (row_p, col_va_j)
                haskey(entry_map, key) && (J_V[entry_map[key]] = H_ij)
            end
            if row_p > 0 && col_vm_j > 0
                key = (row_p, col_vm_j)
                haskey(entry_map, key) && (J_V[entry_map[key]] = N_ij)
            end
            
            # Q row
            if row_q > 0 && col_va_j > 0
                key = (row_q, col_va_j)
                haskey(entry_map, key) && (J_V[entry_map[key]] = J_Qθ_ij)
            end
            if row_q > 0 && col_vm_j > 0
                key = (row_q, col_vm_j)
                haskey(entry_map, key) && (J_V[entry_map[key]] = L_ij)
            end
        end
    end
    
    # DC Jacobian block (using actual DC bus types from workspace)
    dc_non_slack = ws.dc_non_slack
    dc_isolated = ws.dc_isolated
    dc_slack_idx = ws.dc_slack_idx
    dc_bus_to_vdc_col = ws.dc_bus_to_vdc_col
    
    if ndc_eq > 0 && sys.Gdc !== nothing
        Gdc = sys.Gdc
        # Iterate over DC non-slack buses
        for (ki, k) in enumerate(dc_non_slack)
            for (li, l) in enumerate(dc_non_slack)
                row = np + nq + ki
                col = np + nq + li
                
                if k == l
                    # Diagonal: ∂Pdc_k/∂Vdc_k = 2·G_kk·Vdc_k + Σ_{m≠k} G_km·Vdc_m
                    dPk = 2.0 * Gdc[k,k] * Vdc[k]
                    for m in 1:ndc
                        m == k && continue
                        dPk += Gdc[k,m] * Vdc[m]
                    end
                    key = (row, col)
                    haskey(entry_map, key) && (J_V[entry_map[key]] = dPk)
                else
                    # Off-diagonal: ∂Pdc_k/∂Vdc_l = G_kl·Vdc_k
                    key = (row, col)
                    haskey(entry_map, key) && (J_V[entry_map[key]] = Gdc[k,l] * Vdc[k])
                end
            end
        end
        
        # Converter Jacobian contributions using static arrays (droop modes)
        loss_model = sys.loss_model
        for conv in sys.converters
            conv.in_service || continue
            (conv.control_mode == VDC_Q || conv.control_mode == VDC_VAC) || continue
            k = conv.bus_dc
            k == dc_slack_idx && continue
            k ∈ dc_isolated && continue
            col_idx = dc_bus_to_vdc_col[k]
            col_idx == 0 && continue
            
            # Compute both DC diagonal and AC-DC coupling with static arrays
            jac_block = compute_converter_jacobian_block(conv, Vdc, sys.baseMVA, loss_model)
            dPdc_diag = jac_block[1]   # DC diagonal contribution
            dPac_dVdc = jac_block[2]   # AC-DC cross-coupling
            
            # Add DC diagonal contribution
            key = (col_idx, col_idx)
            if haskey(entry_map, key)
                J_V[entry_map[key]] += dPdc_diag
            end
            
            # Add AC-DC cross-coupling (subtract because injected power convention)
            row_p = bus_to_p_row[conv.bus_ac]
            if row_p > 0
                key = (row_p, col_idx)
                if haskey(entry_map, key)
                    J_V[entry_map[key]] -= dPac_dVdc
                end
            end
        end
    end
    
    # Update CSC nzval in-place using pre-computed mapping (no allocation)
    nzval = ws.J_sparse.nzval
    coo_to_csc = ws.coo_to_csc
    @inbounds for k in 1:ws.J_nnz
        nzval[coo_to_csc[k]] = J_V[k]
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
#  NEWTON-RAPHSON SOLVER (OPTIMIZED)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    solve_power_flow(sys::SolverData; max_iter=50, tol=1e-8, init=nothing)

Internal method - users should call `solve_power_flow(hps::HybridPowerSystem)` instead.

Solve hybrid AC/DC power flow using Newton-Raphson with:
- Sparse Jacobian and sparse LU factorization
- Pre-allocated workspace (zero allocation in main loop)
- Symbolic factorization reuse across iterations

Keyword arguments:
- `max_iter`: Maximum NR iterations
- `tol`: Infinity-norm convergence tolerance on mismatch vector
- `init`: Optional warm-start state as `(Vm=..., Va=..., Vdc=...)`
"""
function solve_power_flow(sys::SolverData; max_iter::Int=50, tol::Float64=1e-8,
                          init::Union{Nothing,NamedTuple}=nothing)
    rebuild_matrices!(sys)
    ws = create_solver_workspace(sys)
    
    nac = ws.nac
    ndc = ws.ndc
    nf = ws.nf
    np = ws.np
    nq = ws.nq
    ndc_eq = ws.ndc_eq
    
    # Initialize from bus data or optional warm-start state
    if init === nothing
        @inbounds for i in 1:nac
            ws.Vm[i] = sys.ac_buses[i].vm_pu
            ws.Va[i] = deg2rad(sys.ac_buses[i].va_deg)
        end
        @inbounds for i in 1:ndc
            ws.Vdc[i] = sys.dc_buses[i].vm_pu
        end
    else
        hasproperty(init, :Vm) || error("init must contain field Vm")
        hasproperty(init, :Va) || error("init must contain field Va")
        hasproperty(init, :Vdc) || error("init must contain field Vdc")
        length(init.Vm) == nac || error("length(init.Vm) must be $nac, got $(length(init.Vm))")
        length(init.Va) == nac || error("length(init.Va) must be $nac, got $(length(init.Va))")
        length(init.Vdc) == ndc || error("length(init.Vdc) must be $ndc, got $(length(init.Vdc))")
        copyto!(ws.Vm, init.Vm)
        copyto!(ws.Va, init.Va)
        copyto!(ws.Vdc, init.Vdc)
    end
    
    # Set VDC_VAC controlled bus voltages
    for conv in sys.converters
        if conv.in_service && conv.control_mode == VDC_VAC
            ws.Vm[conv.bus_ac] = conv.v_ac_set_pu
        end
    end
    
    Y = sys.Ybus
    baseMVA = sys.baseMVA
    
    for iter in 1:max_iter
        # Compute power injections
        compute_power_injections!(ws, Y)
        
        # Compute scheduled injections (p.u.)
        @inbounds for i in 1:nac
            ws.Psch[i] = sys.Pg[i] - sys.ac_buses[i].pd_mw / baseMVA
            ws.Qsch[i] = sys.Qg[i] - sys.ac_buses[i].qd_mvar / baseMVA
        end
        
        # Add converter contributions
        for conv in sys.converters
            conv.in_service || continue
            if conv.control_mode == VDC_VAC
                ws.Vm[conv.bus_ac] = conv.v_ac_set_pu
            end
            Pac, Qac = converter_ac_injection(conv, ws.Vm, ws.Va, ws.Vdc, baseMVA, sys.loss_model)
            ws.Psch[conv.bus_ac] += Pac
            ws.Qsch[conv.bus_ac] += Qac
        end
        
        # Build residual vector
        fill!(ws.F, 0.0)
        @inbounds for (k, i) in enumerate(ws.non_slack)
            ws.F[k] = ws.Pcalc[i] - ws.Psch[i]
        end
        @inbounds for (k, i) in enumerate(ws.pq_idx)
            ws.F[np + k] = ws.Qcalc[i] - ws.Qsch[i]
        end
        
        # DC residual (using actual DC bus types)
        if ndc_eq > 0 && ndc > 0 && sys.Gdc !== nothing
            mul!(ws.Pdc_calc, sys.Gdc, ws.Vdc)
            @inbounds for i in 1:ndc
                ws.Pdc_calc[i] *= ws.Vdc[i]
            end
            @inbounds for i in 1:ndc
                ws.Pdc_sch[i] = sys.dc_buses[i].pd_mw / baseMVA
            end
            for conv in sys.converters
                conv.in_service || continue
                ws.Pdc_sch[conv.bus_dc] += converter_dc_injection(conv, ws.Vm, ws.Va, ws.Vdc, baseMVA, sys.loss_model)
            end
            @inbounds for (k, dc_bus) in enumerate(ws.dc_non_slack)
                ws.F[np + nq + k] = ws.Pdc_calc[dc_bus] - ws.Pdc_sch[dc_bus]
            end
        end
        
        res_norm = norm(ws.F, Inf)
        if res_norm < tol
            return (Vm=copy(ws.Vm), Va=copy(ws.Va), Vdc=copy(ws.Vdc), 
                    converged=true, iterations=iter, residual=res_norm)
        end
        
        # Build Jacobian
        build_jacobian_triplets!(ws, sys)
        
        # Solve linear system with sparse LU
        try
            if ws.lu_factor === nothing
                ws.lu_factor = lu(ws.J_sparse)
            else
                # Try to reuse symbolic factorization
                try
                    lu!(ws.lu_factor, ws.J_sparse)
                catch
                    ws.lu_factor = lu(ws.J_sparse)
                end
            end
            ldiv!(ws.dx, ws.lu_factor, ws.F)
            ws.dx .*= -1.0
        catch e
            return (Vm=copy(ws.Vm), Va=copy(ws.Va), Vdc=copy(ws.Vdc),
                    converged=false, iterations=iter, residual=res_norm)
        end
        
        # Check for NaN/Inf
        if any(!isfinite, ws.dx)
            return (Vm=copy(ws.Vm), Va=copy(ws.Va), Vdc=copy(ws.Vdc),
                    converged=false, iterations=iter, residual=res_norm)
        end
        
        # Backtracking line search
        copyto!(ws.Vm_old, ws.Vm)
        copyto!(ws.Va_old, ws.Va)
        copyto!(ws.Vdc_old, ws.Vdc)
        
        α = 1.0
        best_α = α
        best_res = Inf
        
        for _ls in 1:8
            @inbounds for (k, i) in enumerate(ws.non_slack)
                ws.Va[i] = ws.Va_old[i] + α * ws.dx[k]
            end
            @inbounds for (k, i) in enumerate(ws.pq_idx)
                ws.Vm[i] = ws.Vm_old[i] + α * ws.dx[np + k]
            end
            if ndc_eq > 0
                @inbounds for (k, dc_bus) in enumerate(ws.dc_non_slack)
                    ws.Vdc[dc_bus] = ws.Vdc_old[dc_bus] + α * ws.dx[np + nq + k]
                end
            end
            
            @inbounds for i in 1:nac
                ws.Vm[i] = max(ws.Vm[i], 0.05)
            end
            @inbounds for i in 1:ndc
                ws.Vdc[i] = max(ws.Vdc[i], 0.05)
            end
            
            # Evaluate residual
            compute_power_injections!(ws, Y)
            @inbounds for i in 1:nac
                ws.Psch[i] = sys.Pg[i] - sys.ac_buses[i].pd_mw / baseMVA
                ws.Qsch[i] = sys.Qg[i] - sys.ac_buses[i].qd_mvar / baseMVA
            end
            for conv in sys.converters
                conv.in_service || continue
                Pac, Qac = converter_ac_injection(conv, ws.Vm, ws.Va, ws.Vdc, baseMVA, sys.loss_model)
                ws.Psch[conv.bus_ac] += Pac
                ws.Qsch[conv.bus_ac] += Qac
            end
            
            new_res = 0.0
            @inbounds for (k, i) in enumerate(ws.non_slack)
                new_res = max(new_res, abs(ws.Pcalc[i] - ws.Psch[i]))
            end
            @inbounds for (k, i) in enumerate(ws.pq_idx)
                new_res = max(new_res, abs(ws.Qcalc[i] - ws.Qsch[i]))
            end
            if ndc_eq > 0 && sys.Gdc !== nothing
                mul!(ws.Pdc_calc, sys.Gdc, ws.Vdc)
                @inbounds for i in 1:ndc
                    ws.Pdc_calc[i] *= ws.Vdc[i]
                end
                @inbounds for i in 1:ndc
                    ws.Pdc_sch[i] = sys.dc_buses[i].pd_mw / baseMVA
                end
                for conv in sys.converters
                    conv.in_service || continue
                    ws.Pdc_sch[conv.bus_dc] += converter_dc_injection(conv, ws.Vm, ws.Va, ws.Vdc, baseMVA, sys.loss_model)
                end
                @inbounds for (k, dc_bus) in enumerate(ws.dc_non_slack)
                    new_res = max(new_res, abs(ws.Pdc_calc[dc_bus] - ws.Pdc_sch[dc_bus]))
                end
            end
            
            if new_res < best_res
                best_res = new_res
                best_α = α
            end
            new_res < res_norm && break
            α *= 0.5
        end
        
        # Apply best step if no α gave strict decrease
        if best_res >= res_norm
            α = best_α
            @inbounds for (k, i) in enumerate(ws.non_slack)
                ws.Va[i] = ws.Va_old[i] + α * ws.dx[k]
            end
            @inbounds for (k, i) in enumerate(ws.pq_idx)
                ws.Vm[i] = ws.Vm_old[i] + α * ws.dx[np + k]
            end
            if ndc_eq > 0
                @inbounds for (k, dc_bus) in enumerate(ws.dc_non_slack)
                    ws.Vdc[dc_bus] = ws.Vdc_old[dc_bus] + α * ws.dx[np + nq + k]
                end
            end
            @inbounds for i in 1:nac
                ws.Vm[i] = max(ws.Vm[i], 0.05)
            end
            @inbounds for i in 1:ndc
                ws.Vdc[i] = max(ws.Vdc[i], 0.05)
            end
        end
    end
    
    # Compute final residual
    compute_power_injections!(ws, Y)
    @inbounds for i in 1:nac
        ws.Psch[i] = sys.Pg[i] - sys.ac_buses[i].pd_mw / baseMVA
        ws.Qsch[i] = sys.Qg[i] - sys.ac_buses[i].qd_mvar / baseMVA
    end
    for conv in sys.converters
        conv.in_service || continue
        Pac, Qac = converter_ac_injection(conv, ws.Vm, ws.Va, ws.Vdc, baseMVA, sys.loss_model)
        ws.Psch[conv.bus_ac] += Pac
        ws.Qsch[conv.bus_ac] += Qac
    end
    final_res = 0.0
    @inbounds for (k, i) in enumerate(ws.non_slack)
        final_res = max(final_res, abs(ws.Pcalc[i] - ws.Psch[i]))
    end
    @inbounds for (k, i) in enumerate(ws.pq_idx)
        final_res = max(final_res, abs(ws.Qcalc[i] - ws.Qsch[i]))
    end
    
    return (Vm=copy(ws.Vm), Va=copy(ws.Va), Vdc=copy(ws.Vdc),
            converged=false, iterations=max_iter, residual=final_res)
end

"""
    solve_power_flow(hps::HybridPowerSystem; kwargs...)

Solve hybrid AC/DC power flow using Newton-Raphson method.

This is the main user-facing API. Accepts a `HybridPowerSystem` from JuliaPowerCase,
internally creates solver workspace, solves, and writes results back to the input.

# Example
```julia
using JuliaPowerCase, HybridACDCPowerFlow
hps = case_hybrid_5ac3dc()  # Load a test case
result = solve_power_flow(hps)
```

See `solve_power_flow(::SolverData)` for implementation details.
"""
function solve_power_flow(hps::HybridPowerSystem; kwargs...)
    sys = to_solver_data(hps)
    return solve_power_flow(sys; kwargs...)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  PURE DC POWER FLOW SOLVER  
# ═══════════════════════════════════════════════════════════════════════════════

"""
    solve_dc_power_flow(sys; max_iter=50, tol=1e-8, init=nothing)

Solve pure DC power flow using Newton-Raphson. This is a specialized solver
for systems with only DC buses and branches (no AC components).

# DC Bus Types
- `DC_V` (3): Voltage reference bus - voltage is fixed, power is solved
- `DC_P` (1): Power-specified bus - power is specified, voltage is solved

# Mathematical Formulation
The DC power flow equations are:
    Pdc_k = Vdc_k × Σⱼ G_kj × Vdc_j  for all buses k

For each DC_P bus, we solve for Vdc such that:
    Pdc_calc_k - Pdc_sch_k = 0

The DC_V bus serves as the voltage reference (like a slack bus in AC).

# Arguments
- `sys::SolverData`: System with DC buses and DC lines (AC should be empty)
- `max_iter::Int=50`: Maximum Newton-Raphson iterations  
- `tol::Float64=1e-8`: Convergence tolerance (infinity norm)
- `init::NamedTuple=nothing`: Optional initial guess as `(Vdc=...,)`

# Returns
Named tuple with:
- `Vdc::Vector{Float64}`: DC bus voltages (p.u.)
- `converged::Bool`: Whether solution converged
- `iterations::Int`: Number of iterations taken
- `residual::Float64`: Final mismatch magnitude

# Example
```julia
# Build a DC-only system using JuliaPowerCase
using JuliaPowerCase, HybridACDCPowerFlow
dc_sys = DCPowerSystem()
# Add DC buses and solve
result = solve_dc_power_flow(to_solver_data(dc_sys))
```
"""
function solve_dc_power_flow(sys::SolverData; max_iter::Int=50, tol::Float64=1e-8,
                              init::Union{Nothing,NamedTuple}=nothing)
    ndc = length(sys.dc_buses)
    ndc == 0 && error("No DC buses in the system")
    
    # Ensure DC admittance matrix is built
    if sys.Gdc === nothing
        rebuild_matrices!(sys)
    end
    sys.Gdc === nothing && error("Failed to build DC conductance matrix")
    
    Gdc = sys.Gdc
    baseMVA = sys.baseMVA
    
    # Identify DC bus types
    dc_slack_idx = 0
    dc_p_idx = Int[]
    
    for (i, bus) in enumerate(sys.dc_buses)
        bt = bus.bus_type
        if bt == DC_V || bt == REF_BUS
            if dc_slack_idx == 0
                dc_slack_idx = i
            end
        else
            push!(dc_p_idx, i)
        end
    end
    
    # Default: if no DC_V bus specified, use bus 1 as reference
    if dc_slack_idx == 0
        dc_slack_idx = 1
        dc_p_idx = collect(2:ndc)
    end
    
    ndc_eq = length(dc_p_idx)  # Number of DC equations (non-slack buses)
    
    if ndc_eq == 0
        # Only one bus or all are references - trivial case
        Vdc = [bus.vm_pu for bus in sys.dc_buses]
        return (Vdc=Vdc, converged=true, iterations=0, residual=0.0)
    end
    
    # Initialize DC voltages
    Vdc = zeros(ndc)
    if init === nothing
        @inbounds for i in 1:ndc
            Vdc[i] = sys.dc_buses[i].vm_pu
        end
    else
        hasproperty(init, :Vdc) || error("init must contain field Vdc")
        copyto!(Vdc, init.Vdc)
    end
    
    # Allocate workspace
    F = zeros(ndc_eq)                    # Mismatch vector
    Pdc_calc = zeros(ndc)                # Calculated DC power
    Pdc_sch = zeros(ndc)                 # Scheduled DC power
    dx = zeros(ndc_eq)                   # Newton step
    Vdc_old = zeros(ndc)                 # For line search
    
    # Build DC Jacobian sparsity pattern
    # For DC power flow: ∂Pdc_k/∂Vdc_l
    # Row k-th equation (for k in dc_p_idx) -> row index in Jacobian
    # Col l-th variable (for l in dc_p_idx) -> col index in Jacobian
    
    dc_bus_to_row = zeros(Int, ndc)
    dc_bus_to_col = zeros(Int, ndc)
    for (k, dc_bus) in enumerate(dc_p_idx)
        dc_bus_to_row[dc_bus] = k
        dc_bus_to_col[dc_bus] = k
    end
    
    # Build Jacobian structure
    rows_G, cols_G, _ = findnz(Gdc)
    I_jac = Int[]
    J_jac = Int[]
    for k in 1:length(rows_G)
        i, j = rows_G[k], cols_G[k]
        i == dc_slack_idx && continue
        j == dc_slack_idx && continue
        row = dc_bus_to_row[i]
        col = dc_bus_to_col[j]
        row == 0 && continue
        col == 0 && continue
        push!(I_jac, row)
        push!(J_jac, col)
    end
    
    V_jac = zeros(length(I_jac))
    J = sparse(I_jac, J_jac, V_jac, ndc_eq, ndc_eq)
    
    # Build entry map for Jacobian updates
    entry_map = Dict{Tuple{Int,Int}, Int}()
    for k in 1:length(I_jac)
        entry_map[(I_jac[k], J_jac[k])] = k
    end
    
    lu_factor = nothing
    
    for iter in 1:max_iter
        # Compute DC power injections: Pdc_k = Vdc_k × Σⱼ G_kj × Vdc_j
        mul!(Pdc_calc, Gdc, Vdc)
        @inbounds for i in 1:ndc
            Pdc_calc[i] *= Vdc[i]
        end
        
        # Scheduled power (load - generation)
        @inbounds for i in 1:ndc
            # Positive Pd means consumption, so scheduled injection = -Pd
            Pdc_sch[i] = -sys.dc_buses[i].pd_mw / baseMVA
        end
        
        # Build mismatch vector
        @inbounds for (k, dc_bus) in enumerate(dc_p_idx)
            F[k] = Pdc_calc[dc_bus] - Pdc_sch[dc_bus]
        end
        
        res_norm = norm(F, Inf)
        if res_norm < tol
            return (Vdc=copy(Vdc), converged=true, iterations=iter, residual=res_norm)
        end
        
        # Build Jacobian
        fill!(V_jac, 0.0)
        for (k_eq, k) in enumerate(dc_p_idx)
            for (l_eq, l) in enumerate(dc_p_idx)
                key = (k_eq, l_eq)
                haskey(entry_map, key) || continue
                
                if k == l
                    # Diagonal: ∂Pdc_k/∂Vdc_k = 2·G_kk·Vdc_k + Σ_{m≠k} G_km·Vdc_m
                    dPk = 2.0 * Gdc[k,k] * Vdc[k]
                    for m in 1:ndc
                        m == k && continue
                        dPk += Gdc[k,m] * Vdc[m]
                    end
                    V_jac[entry_map[key]] = dPk
                else
                    # Off-diagonal: ∂Pdc_k/∂Vdc_l = G_kl·Vdc_k
                    V_jac[entry_map[key]] = Gdc[k,l] * Vdc[k]
                end
            end
        end
        
        # Update sparse matrix values
        for (k, (i, j)) in enumerate(zip(I_jac, J_jac))
            J[i, j] = V_jac[k]
        end
        
        # Solve linear system
        try
            if lu_factor === nothing
                lu_factor = lu(J)
            else
                try
                    lu!(lu_factor, J)
                catch
                    lu_factor = lu(J)
                end
            end
            ldiv!(dx, lu_factor, F)
            dx .*= -1.0
        catch e
            return (Vdc=copy(Vdc), converged=false, iterations=iter, residual=res_norm)
        end
        
        # Check for NaN/Inf
        if any(!isfinite, dx)
            return (Vdc=copy(Vdc), converged=false, iterations=iter, residual=res_norm)
        end
        
        # Simple line search
        copyto!(Vdc_old, Vdc)
        α = 1.0
        
        for _ls in 1:8
            @inbounds for (k, dc_bus) in enumerate(dc_p_idx)
                Vdc[dc_bus] = Vdc_old[dc_bus] + α * dx[k]
            end
            
            # Clamp voltages
            @inbounds for i in 1:ndc
                Vdc[i] = max(Vdc[i], 0.05)
            end
            
            # Check new residual
            mul!(Pdc_calc, Gdc, Vdc)
            @inbounds for i in 1:ndc
                Pdc_calc[i] *= Vdc[i]
            end
            @inbounds for i in 1:ndc
                Pdc_sch[i] = -sys.dc_buses[i].pd_mw / baseMVA
            end
            
            new_res = 0.0
            @inbounds for (k, dc_bus) in enumerate(dc_p_idx)
                new_res = max(new_res, abs(Pdc_calc[dc_bus] - Pdc_sch[dc_bus]))
            end
            
            new_res < res_norm && break
            α *= 0.5
        end
    end
    
    # Did not converge
    return (Vdc=copy(Vdc), converged=false, iterations=max_iter, residual=norm(F, Inf))
end

"""
    solve_dc_power_flow(hps::HybridPowerSystem; kwargs...)

Solve pure DC power flow for a JuliaPowerCase `HybridPowerSystem`.
See `solve_dc_power_flow(::SolverData)` for full documentation.
"""
function solve_dc_power_flow(hps::HybridPowerSystem; kwargs...)
    sys = to_solver_system(hps)
    return solve_dc_power_flow(sys; kwargs...)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  LEGACY INTERFACE (FOR BACKWARD COMPATIBILITY)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    ac_power_injections(sys, Vm, Va)

Compute AC power injections (legacy interface, allocates).
"""
function ac_power_injections(sys::SolverData, Vm::Vector{Float64}, Va::Vector{Float64})
    Y = sys.Ybus === nothing ? build_admittance_matrix(sys) : sys.Ybus
    V = Vm .* cis.(Va)
    I = Y * V
    S = V .* conj.(I)
    return real.(S), imag.(S)
end

"""
    power_flow_residual(sys, Vm, Va, Vdc)

Compute full AC/DC power flow mismatch vector (legacy interface).
"""
function power_flow_residual(sys::SolverData, Vm::Vector{Float64}, Va::Vector{Float64},
                              Vdc::Vector{Float64})
    nac = length(sys.ac_buses)
    ndc = length(sys.dc_buses)

    Pcalc, Qcalc = ac_power_injections(sys, Vm, Va)

    Psch = [sys.Pg[i] - sys.ac_buses[i].pd_mw / sys.baseMVA for i in 1:nac]
    Qsch = [sys.Qg[i] - sys.ac_buses[i].qd_mvar / sys.baseMVA for i in 1:nac]

    for conv in sys.converters
        conv.in_service || continue
        Pac, Qac = converter_ac_injection(conv, Vm, Va, Vdc, sys.baseMVA, sys.loss_model)
        Psch[conv.bus_ac] += Pac
        Qsch[conv.bus_ac] += Qac
    end

    F_P = Float64[]
    F_Q = Float64[]
    pq_idx = Int[]
    pv_idx = Int[]

    for (i, bus) in enumerate(sys.ac_buses)
        if bus.bus_type == PQ
            push!(pq_idx, i)
        elseif bus.bus_type == PV
            push!(pv_idx, i)
        end
    end

    non_slack = sort(union(pq_idx, pv_idx))
    for i in non_slack
        push!(F_P, Pcalc[i] - Psch[i])
    end
    for i in pq_idx
        push!(F_Q, Qcalc[i] - Qsch[i])
    end

    F_dc = Float64[]
    if ndc > 0
        Gdc = sys.Gdc
        Pdc_calc2 = Vdc .* (Gdc * Vdc)
        Pdc_sch = [b.pd_mw / sys.baseMVA for b in sys.dc_buses]
        for conv in sys.converters
            conv.in_service || continue
            Pdc_inj = converter_dc_injection(conv, Vm, Va, Vdc, sys.baseMVA, sys.loss_model)
            Pdc_sch[conv.bus_dc] += Pdc_inj
        end
        dc_non_slack = collect(2:ndc)
        for k in dc_non_slack
            push!(F_dc, Pdc_calc2[k] - Pdc_sch[k])
        end
    end

    return vcat(F_P, F_Q), F_dc
end

"""
Full residual as a single vector for NR solver (legacy).
"""
function full_residual(sys::SolverData, x::Vector{Float64})
    nac = length(sys.ac_buses)
    ndc = length(sys.dc_buses)

    pq_idx = findall(b -> b.bus_type == PQ, sys.ac_buses)
    pv_idx = findall(b -> b.bus_type == PV, sys.ac_buses)
    non_slack = sort(union(pq_idx, pv_idx))
    npq = length(pq_idx)

    n_va = length(non_slack)
    n_vm = npq
    n_vdc = max(0, ndc - 1)

    Va_full = [deg2rad(b.va_deg) for b in sys.ac_buses]
    Vm_full = [b.vm_pu for b in sys.ac_buses]
    Vdc_full = [b.vm_pu for b in sys.dc_buses]

    idx = 0
    for (k, i) in enumerate(non_slack)
        Va_full[i] = x[idx + k]
    end
    idx += n_va
    for (k, i) in enumerate(pq_idx)
        Vm_full[i] = x[idx + k]
    end
    idx += n_vm
    for k in 1:n_vdc
        Vdc_full[k+1] = x[idx + k]
    end

    F_ac, F_dc = power_flow_residual(sys, Vm_full, Va_full, Vdc_full)
    return vcat(F_ac, F_dc)
end

function full_residual_simple(sys, Vm, Va, Vdc)
    F_ac, F_dc = power_flow_residual(sys, Vm, Va, Vdc)
    return vcat(F_ac, F_dc)
end

"""
    full_residual_simple(hps::HybridPowerSystem, Vm, Va, Vdc)

Convenience wrapper for computing residual from HybridPowerSystem.
Internally converts to SolverData.
"""
function full_residual_simple(hps::HybridPowerSystem, Vm, Va, Vdc)
    sys = to_solver_data(hps)
    return full_residual_simple(sys, Vm, Va, Vdc)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  TOPOLOGY MODIFICATION
# ═══════════════════════════════════════════════════════════════════════════════

function remove_ac_branch(sys::SolverData, branch_idx::Int)
    new_branches = copy(sys.ac_branches)
    br = new_branches[branch_idx]
    # Create a copy with in_service=false
    new_branches[branch_idx] = Branch{AC}(
        index=br.index, name=br.name, from_bus=br.from_bus, to_bus=br.to_bus,
        in_service=false, branch_type=br.branch_type, length_km=br.length_km,
        n_parallel=br.n_parallel, v_rated_kv=br.v_rated_kv, s_rated_mva=br.s_rated_mva,
        s_max_mva=br.s_max_mva, r_pu=br.r_pu, x_pu=br.x_pu, b_pu=br.b_pu,
        r_ohm_km=br.r_ohm_km, x_ohm_km=br.x_ohm_km, c_nf_km=br.c_nf_km,
        b_us_km=br.b_us_km, r0_pu=br.r0_pu, x0_pu=br.x0_pu, b0_pu=br.b0_pu,
        c0_nf_km=br.c0_nf_km, rate_a_mva=br.rate_a_mva, rate_b_mva=br.rate_b_mva,
        rate_c_mva=br.rate_c_mva, tap=br.tap, shift_deg=br.shift_deg,
        angmin_deg=br.angmin_deg, angmax_deg=br.angmax_deg,
        mtbf_hours=br.mtbf_hours, mttr_hours=br.mttr_hours,
        t_scheduled_h=br.t_scheduled_h, sw_hours=br.sw_hours, rp_hours=br.rp_hours
    )
    return SolverData(sys.ac_buses, new_branches, sys.dc_buses, sys.dc_branches,
                        sys.converters; generators=sys.generators, baseMVA=sys.baseMVA)
end

function remove_dc_branch(sys::SolverData, branch_idx::Int)
    new_branches = copy(sys.dc_branches)
    br = new_branches[branch_idx]
    new_branches[branch_idx] = Branch{DC}(
        index=br.index, name=br.name, from_bus=br.from_bus, to_bus=br.to_bus,
        in_service=false, branch_type=br.branch_type, length_km=br.length_km,
        n_parallel=br.n_parallel, v_rated_kv=br.v_rated_kv, s_rated_mva=br.s_rated_mva,
        s_max_mva=br.s_max_mva, r_pu=br.r_pu, x_pu=br.x_pu, b_pu=br.b_pu,
        r_ohm_km=br.r_ohm_km, x_ohm_km=br.x_ohm_km, c_nf_km=br.c_nf_km,
        b_us_km=br.b_us_km, r0_pu=br.r0_pu, x0_pu=br.x0_pu, b0_pu=br.b0_pu,
        c0_nf_km=br.c0_nf_km, rate_a_mva=br.rate_a_mva, rate_b_mva=br.rate_b_mva,
        rate_c_mva=br.rate_c_mva, tap=br.tap, shift_deg=br.shift_deg,
        angmin_deg=br.angmin_deg, angmax_deg=br.angmax_deg,
        mtbf_hours=br.mtbf_hours, mttr_hours=br.mttr_hours,
        t_scheduled_h=br.t_scheduled_h, sw_hours=br.sw_hours, rp_hours=br.rp_hours
    )
    return SolverData(sys.ac_buses, sys.ac_branches, sys.dc_buses, new_branches,
                        sys.converters; generators=sys.generators, baseMVA=sys.baseMVA)
end

function remove_converter(sys::SolverData, conv_idx::Int)
    new_convs = copy(sys.converters)
    c = new_convs[conv_idx]
    new_convs[conv_idx] = VSCConverter(
        index=c.index, name=c.name, bus_ac=c.bus_ac, bus_dc=c.bus_dc,
        in_service=false, vsc_type=c.vsc_type, p_rated_mw=c.p_rated_mw,
        vn_ac_kv=c.vn_ac_kv, vn_dc_kv=c.vn_dc_kv, p_mw=c.p_mw, q_mvar=c.q_mvar,
        vm_ac_pu=c.vm_ac_pu, vm_dc_pu=c.vm_dc_pu, pmax_mw=c.pmax_mw, pmin_mw=c.pmin_mw,
        qmax_mvar=c.qmax_mvar, qmin_mvar=c.qmin_mvar, eta=c.eta,
        loss_percent=c.loss_percent, loss_mw=c.loss_mw, controllable=c.controllable,
        control_mode=c.control_mode, p_set_mw=c.p_set_mw, q_set_mvar=c.q_set_mvar,
        v_ac_set_pu=c.v_ac_set_pu, v_dc_set_pu=c.v_dc_set_pu,
        k_vdc=c.k_vdc, k_p=c.k_p, k_q=c.k_q, v_ref_pu=c.v_ref_pu, f_ref_hz=c.f_ref_hz,
        mtbf_hours=c.mtbf_hours, mttr_hours=c.mttr_hours, t_scheduled_h=c.t_scheduled_h
    )
    return SolverData(sys.ac_buses, sys.ac_branches, sys.dc_buses, sys.dc_branches,
                        new_convs; generators=sys.generators, baseMVA=sys.baseMVA)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  RESULT EXTRACTION
# ═══════════════════════════════════════════════════════════════════════════════

function get_bus_voltages(result)
    return result.Vm, result.Va, result.Vdc
end

function get_branch_flows(sys::SolverData, result)
    Vm, Va, _ = get_bus_voltages(result)
    return get_branch_flows(sys, Vm, Va)
end

function get_branch_flows(sys::SolverData, Vm::Vector{Float64}, Va::Vector{Float64})
    n_br = length(sys.ac_branches)
    Pij = zeros(n_br)
    Qij = zeros(n_br)
    @inbounds for (idx, br) in enumerate(sys.ac_branches)
        br.in_service || continue
        i = br.from_bus
        j = br.to_bus
        ys = 1.0 / complex(br.r_pu, br.x_pu)
        yc = im * br.b_pu / 2.0
        tap = br.tap == 0.0 ? 1.0 : br.tap
        Vi = Vm[i] * exp(im * Va[i])
        Vj = Vm[j] * exp(im * Va[j])
        Iij = (ys / tap^2 + yc) * Vi - (ys / tap) * Vj
        Sij = Vi * conj(Iij)
        Pij[idx] = real(Sij)
        Qij[idx] = imag(Sij)
    end
    return Pij, Qij
end

# ═══════════════════════════════════════════════════════════════════════════════
#  GRAPH DATA EXTRACTION FOR GNN
# ═══════════════════════════════════════════════════════════════════════════════

function extract_graph_data(sys::SolverData)
    nac = length(sys.ac_buses)
    ndc = length(sys.dc_buses)
    baseMVA = sys.baseMVA

    ac_degree = zeros(Int, nac)
    for br in sys.ac_branches
        br.in_service || continue
        ac_degree[br.from_bus] += 1
        ac_degree[br.to_bus] += 1
    end

    dc_degree = zeros(Int, ndc)
    for br in sys.dc_branches
        br.in_service || continue
        dc_degree[br.from_bus] += 1
        dc_degree[br.to_bus] += 1
    end

    ac_adm_sum = zeros(nac)
    for br in sys.ac_branches
        br.in_service || continue
        ys = 1.0 / complex(br.r_pu, br.x_pu)
        adm_mag = abs(ys)
        ac_adm_sum[br.from_bus] += adm_mag
        ac_adm_sum[br.to_bus]   += adm_mag
    end
    ac_adm_sum = max.(ac_adm_sum, 1e-10)

    dc_adm_sum = zeros(ndc)
    for br in sys.dc_branches
        br.in_service || continue
        g = 1.0 / br.r_pu
        dc_adm_sum[br.from_bus] += g
        dc_adm_sum[br.to_bus]   += g
    end
    dc_adm_sum = max.(dc_adm_sum, 1e-10)

    ac_node_feats = zeros(nac, 8)
    for (i, b) in enumerate(sys.ac_buses)
        ac_node_feats[i, Int(b.bus_type)] = 1.0
        ac_node_feats[i, 4] = b.pd_mw / baseMVA
        ac_node_feats[i, 5] = b.qd_mvar / baseMVA
        ac_node_feats[i, 6] = sys.Pg[i]
        ac_node_feats[i, 7] = b.vm_pu
        ac_node_feats[i, 8] = Float64(b.area)
    end

    dc_node_feats = zeros(ndc, 3)
    for (i, b) in enumerate(sys.dc_buses)
        dc_node_feats[i, 1] = b.vm_pu
        dc_node_feats[i, 2] = b.pd_mw / baseMVA
        dc_node_feats[i, 3] = 1.0
    end

    ac_edges_src = Int[]
    ac_edges_dst = Int[]
    for br in sys.ac_branches
        br.in_service || continue
        push!(ac_edges_src, br.from_bus)
        push!(ac_edges_dst, br.to_bus)
        push!(ac_edges_src, br.to_bus)
        push!(ac_edges_dst, br.from_bus)
    end

    n_ac_edges = length(ac_edges_src)
    ac_edge_attr = zeros(n_ac_edges, 5)
    ac_edge_status = ones(n_ac_edges)

    eidx = 0
    for br in sys.ac_branches
        br.in_service || continue
        ys = 1.0 / complex(br.r_pu, br.x_pu)
        adm_mag = abs(ys)
        tap = br.tap == 0.0 ? 1.0 : br.tap
        for dir in 1:2
            eidx += 1
            ac_edge_attr[eidx, 1] = real(ys)
            ac_edge_attr[eidx, 2] = imag(ys)
            ac_edge_attr[eidx, 3] = br.b_pu
            ac_edge_attr[eidx, 4] = tap
            src_node = dir == 1 ? br.from_bus : br.to_bus
            ac_edge_attr[eidx, 5] = adm_mag / ac_adm_sum[src_node]
        end
    end

    dc_edges_src = Int[]
    dc_edges_dst = Int[]
    for br in sys.dc_branches
        br.in_service || continue
        push!(dc_edges_src, br.from_bus)
        push!(dc_edges_dst, br.to_bus)
        push!(dc_edges_src, br.to_bus)
        push!(dc_edges_dst, br.from_bus)
    end

    n_dc_edges = length(dc_edges_src)
    dc_edge_attr = zeros(n_dc_edges, 2)
    dc_edge_status = ones(n_dc_edges)

    deidx = 0
    for br in sys.dc_branches
        br.in_service || continue
        g = 1.0 / br.r_pu
        for dir in 1:2
            deidx += 1
            dc_edge_attr[deidx, 1] = g
            src_node = dir == 1 ? br.from_bus : br.to_bus
            dc_edge_attr[deidx, 2] = g / dc_adm_sum[src_node]
        end
    end

    conv_data = NamedTuple{(:ac_bus, :dc_bus, :mode, :Pset, :Qset, :loss_mw, :loss_percent, :eta, :Smax, :status),
                            Tuple{Int, Int, Symbol, Float64, Float64, Float64, Float64, Float64, Float64, Bool}}[]
    for conv in sys.converters
        conv.in_service || continue
        push!(conv_data, (ac_bus=conv.bus_ac, dc_bus=conv.bus_dc,
                          mode=conv.control_mode, 
                          Pset=conv.p_set_mw / baseMVA,
                          Qset=conv.q_set_mvar / baseMVA,
                          loss_mw=conv.loss_mw, loss_percent=conv.loss_percent,
                          eta=conv.eta, Smax=conv.p_rated_mw / baseMVA,
                          status=conv.in_service))
    end

    return (
        ac_nodes = ac_node_feats,
        dc_nodes = dc_node_feats,
        ac_edges_src = ac_edges_src,
        ac_edges_dst = ac_edges_dst,
        ac_edge_attr = ac_edge_attr,
        ac_edge_status = ac_edge_status,
        dc_edges_src = dc_edges_src,
        dc_edges_dst = dc_edges_dst,
        dc_edge_attr = dc_edge_attr,
        dc_edge_status = dc_edge_status,
        converters = conv_data,
        ac_degree = ac_degree,
        dc_degree = dc_degree,
        nac = nac,
        ndc = ndc
    )
end

end  # module PowerSystem
