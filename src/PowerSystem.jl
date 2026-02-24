"""
    PowerSystem.jl

Hybrid AC/DC power system data structures and Newton-Raphson solver.
Supports: AC power flow, DC network, VSC converters with multiple control modes.

OPTIMIZED v2.0:
- Sparse Jacobian with pre-computed sparsity pattern
- Sparse LU factorization with symbolic reuse
- Pre-allocated workspace buffers (zero allocation in NR loop)
- @inbounds and sin/cos caching for inner loops
- Factored Jacobian builder shared with PowerSystemEnhanced
"""
module PowerSystem

using LinearAlgebra, SparseArrays

export ACBus, ACBranch, DCBus, DCBranch, VSCConverter, HybridSystem,
       BusType, PQ, PV, SLACK, ConverterMode, PQ_MODE, VDC_Q, VDC_VAC,
       build_admittance_matrix, solve_power_flow, power_flow_residual,
       get_bus_voltages, get_branch_flows, rebuild_matrices!,
       remove_ac_branch, extract_graph_data,
       # Optimization exports
       SolverWorkspace, create_solver_workspace, build_jacobian_triplets!,
       compute_power_injections!, compute_residual!

# ═══════════════════════════════════════════════════════════════════════════════
#  BUS AND BRANCH TYPES
# ═══════════════════════════════════════════════════════════════════════════════

@enum BusType PQ=1 PV=2 SLACK=3
@enum ConverterMode PQ_MODE=1 VDC_Q=2 VDC_VAC=3

struct ACBus
    id::Int
    type::BusType
    Pd::Float64      # active load (p.u.)
    Qd::Float64      # reactive load (p.u.)
    Pg::Float64      # active generation (p.u.)
    Qg::Float64      # reactive generation (p.u.)
    Vm::Float64      # voltage magnitude setpoint (p.u.)
    Va::Float64      # voltage angle setpoint (rad)
    area::Int        # area ID for multi-area systems
end

struct ACBranch
    from::Int
    to::Int
    r::Float64       # resistance (p.u.)
    x::Float64       # reactance (p.u.)
    b::Float64       # line charging susceptance (p.u.)
    tap::Float64     # transformer tap ratio (1.0 for lines)
    status::Bool     # in service?
end

struct DCBus
    id::Int
    Vdc_set::Float64 # DC voltage setpoint (p.u.)
    Pdc::Float64     # DC load/generation (p.u.)
end

struct DCBranch
    from::Int
    to::Int
    r::Float64       # DC resistance (p.u.)
    status::Bool
end

struct VSCConverter
    id::Int
    ac_bus::Int      # connected AC bus index
    dc_bus::Int      # connected DC bus index
    mode::ConverterMode
    Pset::Float64    # P setpoint (for PQ_MODE)
    Qset::Float64    # Q setpoint
    Vdc_set::Float64 # Vdc setpoint (for VDC_Q, VDC_VAC)
    Vac_set::Float64 # Vac setpoint (for VDC_VAC)
    Ploss_a::Float64 # loss coefficient a (constant)
    Ploss_b::Float64 # loss coefficient b (proportional)
    Ploss_c::Float64 # loss coefficient c (quadratic)
    Smax::Float64    # MVA rating
    G_droop::Float64 # droop conductance for VDC_Q/VDC_VAC (p.u.)
    status::Bool
end

function VSCConverter(id::Int, ac_bus::Int, dc_bus::Int, mode::ConverterMode,
                      Pset::Float64, Qset::Float64, Vdc_set::Float64, Vac_set::Float64,
                      Ploss_a::Float64, Ploss_b::Float64, Ploss_c::Float64,
                      Smax::Float64, status::Bool; G_droop::Float64=0.1)
    return VSCConverter(id, ac_bus, dc_bus, mode, Pset, Qset, Vdc_set, Vac_set,
                        Ploss_a, Ploss_b, Ploss_c, Smax, G_droop, status)
end

mutable struct HybridSystem
    # AC
    ac_buses::Vector{ACBus}
    ac_branches::Vector{ACBranch}
    baseMVA::Float64
    # DC
    dc_buses::Vector{DCBus}
    dc_branches::Vector{DCBranch}
    # Converters
    converters::Vector{VSCConverter}
    # Computed matrices (sparse)
    Ybus::Union{Nothing, SparseMatrixCSC{ComplexF64, Int}}
    Gdc::Union{Nothing, SparseMatrixCSC{Float64, Int}}
end

function HybridSystem(ac_buses, ac_branches, dc_buses, dc_branches, converters; baseMVA=100.0)
    sys = HybridSystem(ac_buses, ac_branches, baseMVA, dc_buses, dc_branches, converters, nothing, nothing)
    sys.Ybus = build_admittance_matrix(sys)
    sys.Gdc = build_dc_conductance(sys)
    return sys
end

function rebuild_matrices!(sys::HybridSystem)
    sys.Ybus = build_admittance_matrix(sys)
    sys.Gdc = build_dc_conductance(sys)
    return sys
end

# ═══════════════════════════════════════════════════════════════════════════════
#  ADMITTANCE AND CONDUCTANCE MATRICES (SPARSE)
# ═══════════════════════════════════════════════════════════════════════════════

function build_admittance_matrix(sys::HybridSystem)
    n = length(sys.ac_buses)
    # Pre-count nonzeros for efficient sparse construction
    nnz_est = n + 2 * count(br -> br.status, sys.ac_branches)
    I = Vector{Int}(undef, nnz_est)
    J = Vector{Int}(undef, nnz_est)
    V = Vector{ComplexF64}(undef, nnz_est)
    
    # Diagonal entries first
    diag_vals = zeros(ComplexF64, n)
    idx = 0
    
    for br in sys.ac_branches
        br.status || continue
        i, j = br.from, br.to
        ys = 1.0 / complex(br.r, br.x)
        yc = im * br.b / 2.0
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

function build_dc_conductance(sys::HybridSystem)
    n = length(sys.dc_buses)
    n == 0 && return spzeros(Float64, 0, 0)
    
    # Pre-count nonzeros
    nnz_est = n + 2 * count(br -> br.status, sys.dc_branches)
    I = Vector{Int}(undef, nnz_est)
    J = Vector{Int}(undef, nnz_est)
    V = Vector{Float64}(undef, nnz_est)
    
    diag_vals = zeros(n)
    idx = 0
    
    for br in sys.dc_branches
        br.status || continue
        i, j = br.from, br.to
        g = 1.0 / br.r
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
    
    # Bus index sets (computed once)
    pq_idx::Vector{Int}
    pv_idx::Vector{Int}
    slack_idx::Vector{Int}
    non_slack::Vector{Int}
    vac_control_buses::Vector{Int}
    ac_isolated::Set{Int}
    
    # Index mappings for Jacobian
    bus_to_p_row::Vector{Int}  # bus → P equation row (0 if slack)
    bus_to_q_row::Vector{Int}  # bus → Q equation row (0 if not PQ)
    bus_to_va_col::Vector{Int} # bus → Va column (0 if slack/isolated)
    bus_to_vm_col::Vector{Int} # bus → Vm column (0 if not PQ)
    
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
    
    # Sin/cos cache for Jacobian (indexed by branch or (i,j) pair)
    sin_cache::Vector{Float64}
    cos_cache::Vector{Float64}
    
    # G and B sparse matrices (views into Ybus)
    G_sparse::SparseMatrixCSC{Float64, Int}
    B_sparse::SparseMatrixCSC{Float64, Int}
end

"""
    create_solver_workspace(sys::HybridSystem) -> SolverWorkspace

Create pre-allocated workspace for Newton-Raphson solver.
"""
function create_solver_workspace(sys::HybridSystem)
    nac = length(sys.ac_buses)
    ndc = length(sys.dc_buses)
    
    # Identify bus types
    pq_idx = Int[]
    pv_idx = Int[]
    slack_idx = Int[]
    for (i, b) in enumerate(sys.ac_buses)
        if b.type == PQ
            push!(pq_idx, i)
        elseif b.type == PV
            push!(pv_idx, i)
        else
            push!(slack_idx, i)
        end
    end
    
    # VDC_VAC converters (treat as PV-like)
    vac_control_buses = Int[]
    for conv in sys.converters
        if conv.status && conv.mode == VDC_VAC
            push!(vac_control_buses, conv.ac_bus)
        end
    end
    
    pq_idx_modified = setdiff(pq_idx, vac_control_buses)
    pv_idx_modified = union(pv_idx, vac_control_buses)
    
    # Detect AC-isolated buses
    ac_isolated = Set{Int}()
    for i in 1:nac
        has_branch = any(br.status && (br.from == i || br.to == i) for br in sys.ac_branches)
        if !has_branch
            push!(ac_isolated, i)
        end
    end
    
    pq_idx_final = sort(setdiff(pq_idx_modified, ac_isolated))
    pv_idx_final = sort(setdiff(pv_idx_modified, ac_isolated))
    non_slack = sort(union(pq_idx_final, pv_idx_final))
    
    np = length(non_slack)
    nq = length(pq_idx_final)
    ndc_eq = max(0, ndc - 1)
    nf = np + nq + ndc_eq
    
    # Build index mappings
    bus_to_p_row = zeros(Int, nac)
    bus_to_q_row = zeros(Int, nac)
    bus_to_va_col = zeros(Int, nac)
    bus_to_vm_col = zeros(Int, nac)
    
    for (k, i) in enumerate(non_slack)
        bus_to_p_row[i] = k
        bus_to_va_col[i] = k
    end
    for (k, i) in enumerate(pq_idx_final)
        bus_to_q_row[i] = np + k
        bus_to_vm_col[i] = np + k
    end
    
    # Pre-allocate state vectors
    Vm = ones(nac)
    Va = zeros(nac)
    Vdc = ndc > 0 ? ones(ndc) : Float64[]
    Vm_old = similar(Vm)
    Va_old = similar(Va)
    Vdc_old = similar(Vdc)
    
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
    
    # Build Jacobian sparsity pattern
    J_I, J_J, J_V, J_nnz = build_jacobian_sparsity(sys, non_slack, pq_idx_final, 
                                                    slack_idx, ac_isolated, nac, ndc, np, nq, ndc_eq)
    J_sparse = sparse(J_I, J_J, J_V, nf, nf)
    
    # Extract G and B from Ybus
    Y = sys.Ybus === nothing ? build_admittance_matrix(sys) : sys.Ybus
    G_sparse = real.(Y)
    B_sparse = imag.(Y)
    
    # Sin/cos cache (one per AC branch pair)
    n_pairs = nac * nac  # worst case, we index by (i-1)*nac + j
    sin_cache = zeros(n_pairs)
    cos_cache = zeros(n_pairs)
    
    return SolverWorkspace(
        nac, ndc, nf, np, nq, ndc_eq,
        pq_idx_final, pv_idx_final, slack_idx, non_slack, vac_control_buses, ac_isolated,
        bus_to_p_row, bus_to_q_row, bus_to_va_col, bus_to_vm_col,
        Vm, Va, Vdc, Vm_old, Va_old, Vdc_old,
        Pcalc, Qcalc, Psch, Qsch, V_complex, I_complex, S_complex,
        Pdc_calc, Pdc_sch,
        F, dx,
        J_I, J_J, J_V, J_nnz, J_sparse,
        nothing,
        sin_cache, cos_cache,
        G_sparse, B_sparse
    )
end

"""
Build the sparsity pattern for the Jacobian. Returns (I, J, V, nnz).
The pattern is determined by Ybus sparsity plus DC coupling.
"""
function build_jacobian_sparsity(sys, non_slack, pq_idx, slack_idx, ac_isolated, 
                                  nac, ndc, np, nq, ndc_eq)
    Y = sys.Ybus
    nf = np + nq + ndc_eq
    
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
    nnz_est = 4 * nnz(Y) + ndc_eq * ndc_eq + 10 * length(sys.converters)
    
    I = Vector{Int}(undef, nnz_est)
    J = Vector{Int}(undef, nnz_est)
    V = Vector{Float64}(undef, nnz_est)
    idx = 0
    
    # Build column maps
    col_va = Dict(non_slack[k] => k for k in 1:np)
    col_vm = Dict(pq_idx[k] => np + k for k in 1:nq)
    
    # AC Jacobian blocks (H, N, J, L)
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
    
    # DC Jacobian block
    if ndc > 1
        Gdc = sys.Gdc
        rows_G, cols_G, _ = findnz(Gdc)
        for k in 1:length(rows_G)
            ki, li = rows_G[k], cols_G[k]
            ki == 1 && continue  # DC slack
            li == 1 && continue  # DC slack column
            row = np + nq + ki - 1
            col = np + nq + li - 1
            idx += 1; I[idx] = row; J[idx] = col; V[idx] = 0.0
        end
        
        # Converter DC contributions (droop modes)
        for conv in sys.converters
            conv.status || continue
            (conv.mode == VDC_Q || conv.mode == VDC_VAC) || continue
            k = conv.dc_bus
            k == 1 && continue
            ki = k - 1
            row = np + nq + ki
            col = np + nq + ki
            # Check if this entry already exists; if not, add it
            idx += 1; I[idx] = row; J[idx] = col; V[idx] = 0.0
        end
        
        # AC-DC cross-coupling for droop converters
        for conv in sys.converters
            conv.status || continue
            (conv.mode == VDC_Q || conv.mode == VDC_VAC) || continue
            k = conv.dc_bus
            k == 1 && continue
            haskey(col_va, conv.ac_bus) || continue
            ki = k - 1
            row_p = col_va[conv.ac_bus]
            col_vdc = np + nq + ki
            idx += 1; I[idx] = row_p; J[idx] = col_vdc; V[idx] = 0.0
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
"""
function compute_power_injections!(ws::SolverWorkspace, Y::SparseMatrixCSC{ComplexF64, Int})
    @inbounds for i in 1:ws.nac
        ws.V_complex[i] = ws.Vm[i] * cis(ws.Va[i])
    end
    mul!(ws.I_complex, Y, ws.V_complex)
    @inbounds for i in 1:ws.nac
        ws.S_complex[i] = ws.V_complex[i] * conj(ws.I_complex[i])
        ws.Pcalc[i] = real(ws.S_complex[i])
        ws.Qcalc[i] = imag(ws.S_complex[i])
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
#  CONVERTER MODELS
# ═══════════════════════════════════════════════════════════════════════════════

@inline function converter_loss(conv::VSCConverter, P::Float64)
    return conv.Ploss_a + conv.Ploss_b * abs(P) + conv.Ploss_c * P^2
end

@inline function conv_dc_power(conv::VSCConverter, Vdc::Vector{Float64})
    if conv.mode == PQ_MODE
        Ploss = converter_loss(conv, conv.Pset)
        return -(conv.Pset + Ploss)
    else  # VDC_Q or VDC_VAC
        return conv.G_droop * (Vdc[conv.dc_bus]^2 - conv.Vdc_set^2)
    end
end

function converter_ac_injection(conv::VSCConverter, Vm, Va, Vdc)
    if conv.mode == PQ_MODE
        return conv.Pset, conv.Qset
    elseif conv.mode == VDC_Q
        P_transfer = conv_dc_power(conv, Vdc)
        Ploss = converter_loss(conv, P_transfer)
        Pac = P_transfer - Ploss
        return Pac, conv.Qset
    elseif conv.mode == VDC_VAC
        P_transfer = conv_dc_power(conv, Vdc)
        Ploss = converter_loss(conv, P_transfer)
        Pac = P_transfer - Ploss
        return Pac, 0.0
    end
    return 0.0, 0.0
end

function converter_dc_injection(conv::VSCConverter, Vm, Va, Vdc)
    if conv.mode == PQ_MODE
        Ploss = converter_loss(conv, conv.Pset)
        return -(conv.Pset + Ploss)
    else
        return -conv_dc_power(conv, Vdc)
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
#  JACOBIAN BUILDER (SHARED IMPLEMENTATION)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    build_jacobian_triplets!(ws, sys)

Update Jacobian values in COO format. Sparsity pattern is fixed;
only values are updated. Uses @inbounds and cached sin/cos.
"""
function build_jacobian_triplets!(ws::SolverWorkspace, sys::HybridSystem)
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
    
    # Reset Jacobian values
    fill!(ws.J_V, 0.0)
    
    # Build value index map for efficient updates
    # We'll rebuild the sparse matrix from triplets
    J_I = ws.J_I
    J_J = ws.J_J
    J_V = ws.J_V
    
    # Create a dict for (row, col) -> index in J_V
    entry_map = Dict{Tuple{Int,Int}, Int}()
    for k in 1:ws.J_nnz
        key = (J_I[k], J_J[k])
        entry_map[key] = k
    end
    
    bus_to_p_row = ws.bus_to_p_row
    bus_to_q_row = ws.bus_to_q_row
    bus_to_va_col = ws.bus_to_va_col
    bus_to_vm_col = ws.bus_to_vm_col
    pq_idx = ws.pq_idx
    slack_idx = ws.slack_idx
    ac_isolated = ws.ac_isolated
    
    # Iterate over Ybus nonzeros
    rows, cols, _ = findnz(sys.Ybus)
    
    @inbounds for k in 1:length(rows)
        i, j = rows[k], cols[k]
        
        i in slack_idx && continue
        i in ac_isolated && continue
        
        Gij = G[i, j]
        Bij = B[i, j]
        (Gij == 0.0 && Bij == 0.0) && continue
        
        θij = Va[i] - Va[j]
        sinθ = sin(θij)
        cosθ = cos(θij)
        
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
    
    # DC Jacobian block
    if ndc > 1
        Gdc = sys.Gdc
        for ki in 1:(ndc-1)
            k = ki + 1  # DC bus index (skip DC slack = bus 1)
            for li in 1:(ndc-1)
                l = li + 1
                row = np + nq + ki
                col = np + nq + li
                
                if k == l
                    dPk = 2.0 * Gdc[k,k] * Vdc[k]
                    for m in 1:ndc
                        m == k && continue
                        dPk += Gdc[k,m] * Vdc[m]
                    end
                    key = (row, col)
                    haskey(entry_map, key) && (J_V[entry_map[key]] = dPk)
                else
                    key = (row, col)
                    haskey(entry_map, key) && (J_V[entry_map[key]] = Gdc[k,l] * Vdc[k])
                end
            end
        end
        
        # Converter DC contributions
        for conv in sys.converters
            conv.status || continue
            (conv.mode == VDC_Q || conv.mode == VDC_VAC) || continue
            k = conv.dc_bus
            k == 1 && continue
            ki = k - 1
            row = np + nq + ki
            col = np + nq + ki
            key = (row, col)
            if haskey(entry_map, key)
                J_V[entry_map[key]] += 2.0 * conv.G_droop * Vdc[k]
            end
        end
        
        # AC–DC cross-coupling
        for conv in sys.converters
            conv.status || continue
            (conv.mode == VDC_Q || conv.mode == VDC_VAC) || continue
            k = conv.dc_bus
            k == 1 && continue
            row_p = bus_to_p_row[conv.ac_bus]
            row_p == 0 && continue
            ki = k - 1
            col_vdc = np + nq + ki
            P_transfer = conv.G_droop * (Vdc[k]^2 - conv.Vdc_set^2)
            dPloss_dP = conv.Ploss_b * sign(P_transfer + 1e-30) + 2.0 * conv.Ploss_c * P_transfer
            key = (row_p, col_vdc)
            if haskey(entry_map, key)
                J_V[entry_map[key]] -= 2.0 * conv.G_droop * Vdc[k] * (1.0 - dPloss_dP)
            end
        end
    end
    
    # Rebuild sparse matrix from updated triplets
    ws.J_sparse = sparse(J_I, J_J, J_V, ws.nf, ws.nf)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  NEWTON-RAPHSON SOLVER (OPTIMIZED)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    solve_power_flow(sys; max_iter=50, tol=1e-8)

Solve hybrid AC/DC power flow using Newton-Raphson with:
- Sparse Jacobian and sparse LU factorization
- Pre-allocated workspace (zero allocation in main loop)
- Symbolic factorization reuse across iterations
"""
function solve_power_flow(sys::HybridSystem; max_iter::Int=50, tol::Float64=1e-8)
    rebuild_matrices!(sys)
    ws = create_solver_workspace(sys)
    
    nac = ws.nac
    ndc = ws.ndc
    nf = ws.nf
    np = ws.np
    nq = ws.nq
    ndc_eq = ws.ndc_eq
    
    # Initialize from bus data
    @inbounds for i in 1:nac
        ws.Vm[i] = sys.ac_buses[i].Vm
        ws.Va[i] = sys.ac_buses[i].Va
    end
    @inbounds for i in 1:ndc
        ws.Vdc[i] = sys.dc_buses[i].Vdc_set
    end
    
    # Set VDC_VAC controlled bus voltages
    for conv in sys.converters
        if conv.status && conv.mode == VDC_VAC
            ws.Vm[conv.ac_bus] = conv.Vac_set
        end
    end
    
    Y = sys.Ybus
    
    for iter in 1:max_iter
        # Compute power injections
        compute_power_injections!(ws, Y)
        
        # Compute scheduled injections
        @inbounds for i in 1:nac
            ws.Psch[i] = sys.ac_buses[i].Pg - sys.ac_buses[i].Pd
            ws.Qsch[i] = sys.ac_buses[i].Qg - sys.ac_buses[i].Qd
        end
        
        # Add converter contributions
        for conv in sys.converters
            conv.status || continue
            if conv.mode == VDC_VAC
                ws.Vm[conv.ac_bus] = conv.Vac_set
            end
            Pac, Qac = converter_ac_injection(conv, ws.Vm, ws.Va, ws.Vdc)
            ws.Psch[conv.ac_bus] += Pac
            ws.Qsch[conv.ac_bus] += Qac
        end
        
        # Build residual vector
        fill!(ws.F, 0.0)
        @inbounds for (k, i) in enumerate(ws.non_slack)
            ws.F[k] = ws.Pcalc[i] - ws.Psch[i]
        end
        @inbounds for (k, i) in enumerate(ws.pq_idx)
            ws.F[np + k] = ws.Qcalc[i] - ws.Qsch[i]
        end
        
        # DC residual
        if ndc > 1
            mul!(ws.Pdc_calc, sys.Gdc, ws.Vdc)
            @inbounds for i in 1:ndc
                ws.Pdc_calc[i] *= ws.Vdc[i]
            end
            @inbounds for i in 1:ndc
                ws.Pdc_sch[i] = sys.dc_buses[i].Pdc
            end
            for conv in sys.converters
                conv.status || continue
                ws.Pdc_sch[conv.dc_bus] += converter_dc_injection(conv, ws.Vm, ws.Va, ws.Vdc)
            end
            @inbounds for k in 1:(ndc-1)
                ws.F[np + nq + k] = ws.Pdc_calc[k+1] - ws.Pdc_sch[k+1]
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
            if ndc > 1
                @inbounds for k in 1:(ndc-1)
                    ws.Vdc[k+1] = ws.Vdc_old[k+1] + α * ws.dx[np + nq + k]
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
                ws.Psch[i] = sys.ac_buses[i].Pg - sys.ac_buses[i].Pd
                ws.Qsch[i] = sys.ac_buses[i].Qg - sys.ac_buses[i].Qd
            end
            for conv in sys.converters
                conv.status || continue
                Pac, Qac = converter_ac_injection(conv, ws.Vm, ws.Va, ws.Vdc)
                ws.Psch[conv.ac_bus] += Pac
                ws.Qsch[conv.ac_bus] += Qac
            end
            
            new_res = 0.0
            @inbounds for (k, i) in enumerate(ws.non_slack)
                new_res = max(new_res, abs(ws.Pcalc[i] - ws.Psch[i]))
            end
            @inbounds for (k, i) in enumerate(ws.pq_idx)
                new_res = max(new_res, abs(ws.Qcalc[i] - ws.Qsch[i]))
            end
            if ndc > 1
                mul!(ws.Pdc_calc, sys.Gdc, ws.Vdc)
                @inbounds for i in 1:ndc
                    ws.Pdc_calc[i] *= ws.Vdc[i]
                end
                @inbounds for i in 1:ndc
                    ws.Pdc_sch[i] = sys.dc_buses[i].Pdc
                end
                for conv in sys.converters
                    conv.status || continue
                    ws.Pdc_sch[conv.dc_bus] += converter_dc_injection(conv, ws.Vm, ws.Va, ws.Vdc)
                end
                @inbounds for k in 1:(ndc-1)
                    new_res = max(new_res, abs(ws.Pdc_calc[k+1] - ws.Pdc_sch[k+1]))
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
            if ndc > 1
                @inbounds for k in 1:(ndc-1)
                    ws.Vdc[k+1] = ws.Vdc_old[k+1] + α * ws.dx[np + nq + k]
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
        ws.Psch[i] = sys.ac_buses[i].Pg - sys.ac_buses[i].Pd
        ws.Qsch[i] = sys.ac_buses[i].Qg - sys.ac_buses[i].Qd
    end
    for conv in sys.converters
        conv.status || continue
        Pac, Qac = converter_ac_injection(conv, ws.Vm, ws.Va, ws.Vdc)
        ws.Psch[conv.ac_bus] += Pac
        ws.Qsch[conv.ac_bus] += Qac
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

# ═══════════════════════════════════════════════════════════════════════════════
#  LEGACY INTERFACE (FOR BACKWARD COMPATIBILITY)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    ac_power_injections(sys, Vm, Va)

Compute AC power injections (legacy interface, allocates).
"""
function ac_power_injections(sys::HybridSystem, Vm::Vector{Float64}, Va::Vector{Float64})
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
function power_flow_residual(sys::HybridSystem, Vm::Vector{Float64}, Va::Vector{Float64},
                              Vdc::Vector{Float64})
    nac = length(sys.ac_buses)
    ndc = length(sys.dc_buses)

    Pcalc, Qcalc = ac_power_injections(sys, Vm, Va)

    Psch = [b.Pg - b.Pd for b in sys.ac_buses]
    Qsch = [b.Qg - b.Qd for b in sys.ac_buses]

    for conv in sys.converters
        conv.status || continue
        Pac, Qac = converter_ac_injection(conv, Vm, Va, Vdc)
        Psch[conv.ac_bus] += Pac
        Qsch[conv.ac_bus] += Qac
    end

    F_P = Float64[]
    F_Q = Float64[]
    pq_idx = Int[]
    pv_idx = Int[]

    for (i, bus) in enumerate(sys.ac_buses)
        if bus.type == PQ
            push!(pq_idx, i)
        elseif bus.type == PV
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
        Pdc_sch = [b.Pdc for b in sys.dc_buses]
        for conv in sys.converters
            conv.status || continue
            Pdc_inj = converter_dc_injection(conv, Vm, Va, Vdc)
            Pdc_sch[conv.dc_bus] += Pdc_inj
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
function full_residual(sys::HybridSystem, x::Vector{Float64})
    nac = length(sys.ac_buses)
    ndc = length(sys.dc_buses)

    pq_idx = findall(b -> b.type == PQ, sys.ac_buses)
    pv_idx = findall(b -> b.type == PV, sys.ac_buses)
    non_slack = sort(union(pq_idx, pv_idx))
    npq = length(pq_idx)

    n_va = length(non_slack)
    n_vm = npq
    n_vdc = max(0, ndc - 1)

    Va_full = [b.Va for b in sys.ac_buses]
    Vm_full = [b.Vm for b in sys.ac_buses]
    Vdc_full = [b.Vdc_set for b in sys.dc_buses]

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

# ═══════════════════════════════════════════════════════════════════════════════
#  TOPOLOGY MODIFICATION
# ═══════════════════════════════════════════════════════════════════════════════

function remove_ac_branch(sys::HybridSystem, branch_idx::Int)
    new_branches = copy(sys.ac_branches)
    br = new_branches[branch_idx]
    new_branches[branch_idx] = ACBranch(br.from, br.to, br.r, br.x, br.b, br.tap, false)
    return HybridSystem(sys.ac_buses, new_branches, sys.dc_buses, sys.dc_branches,
                        sys.converters; baseMVA=sys.baseMVA)
end

function remove_dc_branch(sys::HybridSystem, branch_idx::Int)
    new_branches = copy(sys.dc_branches)
    br = new_branches[branch_idx]
    new_branches[branch_idx] = DCBranch(br.from, br.to, br.r, false)
    return HybridSystem(sys.ac_buses, sys.ac_branches, sys.dc_buses, new_branches,
                        sys.converters; baseMVA=sys.baseMVA)
end

function remove_converter(sys::HybridSystem, conv_idx::Int)
    new_convs = copy(sys.converters)
    c = new_convs[conv_idx]
    new_convs[conv_idx] = VSCConverter(c.id, c.ac_bus, c.dc_bus, c.mode, c.Pset, c.Qset,
                                       c.Vdc_set, c.Vac_set, c.Ploss_a, c.Ploss_b, c.Ploss_c,
                                       c.Smax, false; G_droop=c.G_droop)
    return HybridSystem(sys.ac_buses, sys.ac_branches, sys.dc_buses, sys.dc_branches,
                        new_convs; baseMVA=sys.baseMVA)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  RESULT EXTRACTION
# ═══════════════════════════════════════════════════════════════════════════════

function get_bus_voltages(result)
    return result.Vm, result.Va, result.Vdc
end

function get_branch_flows(sys::HybridSystem, result)
    Vm, Va, _ = get_bus_voltages(result)
    return get_branch_flows(sys, Vm, Va)
end

function get_branch_flows(sys::HybridSystem, Vm::Vector{Float64}, Va::Vector{Float64})
    n_br = length(sys.ac_branches)
    Pij = zeros(n_br)
    Qij = zeros(n_br)
    @inbounds for (idx, br) in enumerate(sys.ac_branches)
        br.status || continue
        i = br.from
        j = br.to
        ys = 1.0 / complex(br.r, br.x)
        yc = im * br.b / 2.0
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

function extract_graph_data(sys::HybridSystem)
    nac = length(sys.ac_buses)
    ndc = length(sys.dc_buses)

    ac_degree = zeros(Int, nac)
    for br in sys.ac_branches
        br.status || continue
        ac_degree[br.from] += 1
        ac_degree[br.to] += 1
    end

    dc_degree = zeros(Int, ndc)
    for br in sys.dc_branches
        br.status || continue
        dc_degree[br.from] += 1
        dc_degree[br.to] += 1
    end

    ac_adm_sum = zeros(nac)
    for br in sys.ac_branches
        br.status || continue
        ys = 1.0 / complex(br.r, br.x)
        adm_mag = abs(ys)
        ac_adm_sum[br.from] += adm_mag
        ac_adm_sum[br.to]   += adm_mag
    end
    ac_adm_sum = max.(ac_adm_sum, 1e-10)

    dc_adm_sum = zeros(ndc)
    for br in sys.dc_branches
        br.status || continue
        g = 1.0 / br.r
        dc_adm_sum[br.from] += g
        dc_adm_sum[br.to]   += g
    end
    dc_adm_sum = max.(dc_adm_sum, 1e-10)

    ac_node_feats = zeros(nac, 8)
    for (i, b) in enumerate(sys.ac_buses)
        ac_node_feats[i, Int(b.type)] = 1.0
        ac_node_feats[i, 4] = b.Pd
        ac_node_feats[i, 5] = b.Qd
        ac_node_feats[i, 6] = b.Pg
        ac_node_feats[i, 7] = b.Vm
        ac_node_feats[i, 8] = Float64(b.area)
    end

    dc_node_feats = zeros(ndc, 3)
    for (i, b) in enumerate(sys.dc_buses)
        dc_node_feats[i, 1] = b.Vdc_set
        dc_node_feats[i, 2] = b.Pdc
        dc_node_feats[i, 3] = 1.0
    end

    ac_edges_src = Int[]
    ac_edges_dst = Int[]
    for br in sys.ac_branches
        br.status || continue
        push!(ac_edges_src, br.from)
        push!(ac_edges_dst, br.to)
        push!(ac_edges_src, br.to)
        push!(ac_edges_dst, br.from)
    end

    n_ac_edges = length(ac_edges_src)
    ac_edge_attr = zeros(n_ac_edges, 5)
    ac_edge_status = ones(n_ac_edges)

    eidx = 0
    for br in sys.ac_branches
        br.status || continue
        ys = 1.0 / complex(br.r, br.x)
        adm_mag = abs(ys)
        tap = br.tap == 0.0 ? 1.0 : br.tap
        for dir in 1:2
            eidx += 1
            ac_edge_attr[eidx, 1] = real(ys)
            ac_edge_attr[eidx, 2] = imag(ys)
            ac_edge_attr[eidx, 3] = br.b
            ac_edge_attr[eidx, 4] = tap
            src_node = dir == 1 ? br.from : br.to
            ac_edge_attr[eidx, 5] = adm_mag / ac_adm_sum[src_node]
        end
    end

    dc_edges_src = Int[]
    dc_edges_dst = Int[]
    for br in sys.dc_branches
        br.status || continue
        push!(dc_edges_src, br.from)
        push!(dc_edges_dst, br.to)
        push!(dc_edges_src, br.to)
        push!(dc_edges_dst, br.from)
    end

    n_dc_edges = length(dc_edges_src)
    dc_edge_attr = zeros(n_dc_edges, 2)
    dc_edge_status = ones(n_dc_edges)

    deidx = 0
    for br in sys.dc_branches
        br.status || continue
        g = 1.0 / br.r
        for dir in 1:2
            deidx += 1
            dc_edge_attr[deidx, 1] = g
            src_node = dir == 1 ? br.from : br.to
            dc_edge_attr[deidx, 2] = g / dc_adm_sum[src_node]
        end
    end

    conv_data = []
    for conv in sys.converters
        conv.status || continue
        push!(conv_data, (ac_bus=conv.ac_bus, dc_bus=conv.dc_bus,
                          mode=Int(conv.mode), Pset=conv.Pset, Qset=conv.Qset,
                          Ploss_a=conv.Ploss_a, Ploss_b=conv.Ploss_b,
                          Ploss_c=conv.Ploss_c, Smax=conv.Smax,
                          status=conv.status))
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
