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
Collect all content between the first '[' and the matching ']' in `lines`
starting from `start_line`. Returns the raw content string.
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
            break
        else
            write(content, raw, " ")
        end
    end
    return String(take!(content))
end

"""
Parse a block content string (between '[' and ']') into a Matrix{Float64}.
Semicolons delimit rows; whitespace/commas delimit values within a row.
"""
function _parse_block_content(content::String)
    rows = Vector{Float64}[]
    # Replace commas with spaces and split on semicolons
    for row_str in split(content, ';')
        # Replace commas used as separators with spaces
        row_str = replace(row_str, ',' => ' ')
        row_str = strip(row_str)
        isempty(row_str) && continue
        # Handle MATLAB line continuation dots
        row_str = replace(row_str, "..." => " ")
        vals = [tryparse(Float64, s) for s in split(row_str) if !isempty(s)]
        parsed = Float64[v for v in vals if v !== nothing]
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

        # baseMVA
        if occursin(r"\.baseMVA\s*=", raw)
            m = match(r"=\s*([0-9.eE+\-]+)", raw)
            m !== nothing && (baseMVA = parse(Float64, m.captures[1]))

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

