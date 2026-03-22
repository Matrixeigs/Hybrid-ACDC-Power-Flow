"""
    MatpowerParser

Parses MATPOWER .m case files into Julia data structures.
"""
module MatpowerParser

export MatpowerData, parse_matpower

"""
    MatpowerData

Parsed MATPOWER case data with sequential bus indexing.

All bus references in `gen` and `branch` matrices are remapped to sequential
1:N indices matching the row order of the `bus` matrix.
"""
struct MatpowerData
    bus::Matrix{Float64}          # n_bus × ncols;  column 1 = original bus ID
    gen::Matrix{Float64}          # n_gen × ncols;  column 1 = sequential bus index
    branch::Matrix{Float64}       # n_br  × ncols;  columns 1,2 = sequential bus indices
    baseMVA::Float64
    busid_to_idx::Dict{Int, Int}  # original bus ID → sequential (1-based) index
end

# ─────────────────────────────────────────────────────────────────────────────
#  INTERNAL HELPERS
# ─────────────────────────────────────────────────────────────────────────────

"""Strip MATLAB-style line comment (everything from % to end of line)."""
_strip_comment(s::AbstractString) = let i = findfirst('%', s); i === nothing ? s : s[1:i-1] end

"""
Evaluate simple MATLAB math expressions like `135/sqrt(3)`, `50/3`, `sqrt(3)`.
Returns the Float64 value or `nothing` if the expression is not recognized.
"""
function _eval_simple_expr(s::AbstractString)::Union{Float64, Nothing}
    # Direct parse
    v = tryparse(Float64, s)
    v !== nothing && return v
    # a/sqrt(b)
    m = match(r"^(-?[\d.eE+\-]+)/sqrt\((-?[\d.eE+\-]+)\)$", s)
    if m !== nothing
        a = tryparse(Float64, m.captures[1])
        b = tryparse(Float64, m.captures[2])
        a !== nothing && b !== nothing && b > 0 && return a / sqrt(b)
    end
    # sqrt(a)/b
    m = match(r"^sqrt\((-?[\d.eE+\-]+)\)/(-?[\d.eE+\-]+)$", s)
    if m !== nothing
        a = tryparse(Float64, m.captures[1])
        b = tryparse(Float64, m.captures[2])
        a !== nothing && b !== nothing && b != 0 && return sqrt(a) / b
    end
    # sqrt(a)
    m = match(r"^sqrt\((-?[\d.eE+\-]+)\)$", s)
    if m !== nothing
        a = tryparse(Float64, m.captures[1])
        a !== nothing && return sqrt(a)
    end
    # a/b
    m = match(r"^(-?[\d.eE+\-]+)/(-?[\d.eE+\-]+)$", s)
    if m !== nothing
        a = tryparse(Float64, m.captures[1])
        b = tryparse(Float64, m.captures[2])
        a !== nothing && b !== nothing && b != 0 && return a / b
    end
    return nothing
end

"""
Collect all content between the first '[' and the matching ']' in `lines`
starting from `start_line`. Returns the raw content string with newlines preserved.
"""
function _collect_block_content(lines::Vector{<:AbstractString}, start_line::Int)
    content = IOBuffer()
    found_open = false
    for i in start_line:length(lines)
        raw = _strip_comment(lines[i])
        if !found_open
            pos = findfirst('[', raw)
            pos === nothing && continue
            raw = raw[pos+1:end]
            found_open = true
        end
        close_pos = findfirst(']', raw)
        if close_pos !== nothing
            write(content, raw[1:close_pos-1])
            write(content, "\n")  # ensure last line is terminated
            break
        else
            write(content, raw, "\n")  # preserve newlines as row delimiters
        end
    end
    return String(take!(content))
end

"""
Parse a block content string (between '[' and ']') into a Matrix{Float64}.
Both semicolons and newlines delimit rows; whitespace/commas delimit values
within a row. Simple MATLAB math expressions (e.g. `135/sqrt(3)`) are evaluated.
"""
function _parse_block_content(content::String)
    rows = Vector{Float64}[]
    # Split on BOTH semicolons and newlines so that rows without trailing
    # semicolons (valid in MATPOWER, e.g. bus 1 in some files) are handled.
    for row_str in split(content, r"[;\n]")
        # Replace commas used as separators with spaces
        row_str = replace(row_str, ',' => ' ')
        row_str = strip(row_str)
        isempty(row_str) && continue
        # Handle MATLAB line continuation dots
        row_str = replace(row_str, "..." => " ")
        # Evaluate each token, including simple math expressions
        parsed = Float64[]
        for s in split(row_str)
            isempty(s) && continue
            v = _eval_simple_expr(s)
            v !== nothing && push!(parsed, v)
        end
        isempty(parsed) && continue
        push!(rows, parsed)
    end
    isempty(rows) && return Matrix{Float64}(undef, 0, 0)
    ncols = maximum(length(r) for r in rows)
    M = zeros(Float64, length(rows), ncols)
    for (ri, r) in enumerate(rows)
        M[ri, 1:length(r)] .= r
    end
    return M
end

# ─────────────────────────────────────────────────────────────────────────────
#  MAIN PARSER
# ─────────────────────────────────────────────────────────────────────────────

"""
    parse_matpower(filepath::String) -> MatpowerData

Parse a MATPOWER .m case file and return a `MatpowerData` struct.

Bus IDs are remapped to sequential 1:N indices. The `busid_to_idx` field
provides the mapping from original bus IDs to sequential indices.
"""
function parse_matpower(filepath::String)
    isfile(filepath) || error("File not found: $filepath")
    lines = readlines(filepath)

    baseMVA = 100.0
    bus_mat    = Matrix{Float64}(undef, 0, 0)
    gen_mat    = Matrix{Float64}(undef, 0, 0)
    branch_mat = Matrix{Float64}(undef, 0, 0)

    i = 1
    while i <= length(lines)
        raw = strip(_strip_comment(lines[i]))

        # baseMVA — supports plain numbers and simple a/b expressions (e.g. 50/3)
        if occursin(r"\.baseMVA\s*=", raw)
            m = match(r"=\s*([^\s;]+)", raw)
            if m !== nothing
                v = _eval_simple_expr(m.captures[1])
                v !== nothing && (baseMVA = v)
            end

        # bus matrix
        elseif occursin(r"\.bus\s*=", raw) && !occursin(r"busname", raw)
            content = _collect_block_content(lines, i)
            bus_mat = _parse_block_content(content)

        # branch matrix
        elseif occursin(r"\.branch\s*=", raw)
            content = _collect_block_content(lines, i)
            branch_mat = _parse_block_content(content)

        # gen matrix (not gencost, genfuel, gentype)
        elseif occursin(r"\.gen\s*=", raw) &&
               !occursin("gencost", raw) && !occursin("genfuel", raw) &&
               !occursin("gentype", raw)
            content = _collect_block_content(lines, i)
            gen_mat = _parse_block_content(content)
        end

        i += 1
    end

    size(bus_mat, 1) == 0 && error("No bus data found in $filepath")

    # ── Apply MATLAB post-processing conversions ──────────────────────────────
    # Some distribution MATPOWER files store loads in kW/kVAr and branch
    # impedances in Ohms, with explicit MATLAB conversion code.  We detect
    # those patterns and apply the same conversions in Julia.

    # Pattern 1: loads from kW/kVAr → MW/MVAR
    #   mpc.bus(:, [PD, QD]) = mpc.bus(:, [PD, QD]) / 1e3
    load_kw = any(l -> occursin(r"bus\(.*\[PD|bus\(:,.*PD|bus\s*\(:,\s*\[PD", l) &&
                       occursin(r"/\s*1[Ee]3|/\s*1000", l),
              lines)
    if load_kw && size(bus_mat, 2) >= 4
        bus_mat[:, 3] ./= 1000.0   # Pd: kW → MW
        bus_mat[:, 4] ./= 1000.0   # Qd: kVAr → MVAR
    end

    # Pattern 2: branch R, X from Ohms → per-unit
    #   Vbase = mpc.bus(1, BASE_KV) * 1e3     (V)
    #   Sbase = mpc.baseMVA * 1e6             (VA)
    #   mpc.branch(:, [BR_R BR_X]) = ... / (Vbase^2 / Sbase)
    #   ⟹  Z_base_pu = V_base_kV^2 / baseMVA
    branch_ohm = any(l -> occursin(r"branch\(.*\[BR_R|branch\(:.*BR_R|Vbase\^2\s*/\s*Sbase", l),
                 lines)
    if branch_ohm && size(bus_mat, 2) >= 10 && size(branch_mat, 2) >= 4
        bus1_baseKV = bus_mat[1, 10]  # column 10 = baseKV of first bus
        if bus1_baseKV > 0
            Z_base = bus1_baseKV^2 / baseMVA   # Ω base
            branch_mat[:, 3] ./= Z_base   # R: Ω → pu
            branch_mat[:, 4] ./= Z_base   # X: Ω → pu
        end
    end

    nbus = size(bus_mat, 1)

    # Build bus ID → sequential index mapping
    busid_to_idx = Dict{Int, Int}()
    for k in 1:nbus
        orig_id = Int(round(bus_mat[k, 1]))
        busid_to_idx[orig_id] = k
    end

    # If no generators found, create empty matrix
    if size(gen_mat, 1) == 0
        gen_mat = Matrix{Float64}(undef, 0, 10)
    end

    # Remap generator bus column from original ID to sequential index
    gen_remapped = copy(gen_mat)
    for g in 1:size(gen_mat, 1)
        orig_id = Int(round(gen_mat[g, 1]))
        idx = get(busid_to_idx, orig_id, 0)
        gen_remapped[g, 1] = Float64(idx)
    end

    # If no branches found, create empty matrix
    if size(branch_mat, 1) == 0
        branch_mat = Matrix{Float64}(undef, 0, 11)
    end

    # Remap branch from/to columns from original IDs to sequential indices
    branch_remapped = copy(branch_mat)
    for k in 1:size(branch_mat, 1)
        f_orig = Int(round(branch_mat[k, 1]))
        t_orig = Int(round(branch_mat[k, 2]))
        branch_remapped[k, 1] = Float64(get(busid_to_idx, f_orig, 0))
        branch_remapped[k, 2] = Float64(get(busid_to_idx, t_orig, 0))
    end

    return MatpowerData(bus_mat, gen_remapped, branch_remapped, baseMVA, busid_to_idx)
end

end # module

