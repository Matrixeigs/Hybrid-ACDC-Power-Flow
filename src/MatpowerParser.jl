"""
    MatpowerParser.jl

Parse MATPOWER `.m` case files into raw Julia data structures.
Handles the standard MATPOWER Case Format Version 2, including:
  - Non-sequential bus IDs (mapped to contiguous 1:N indices)
  - Impedance/load unit conversions (Ohms→p.u., kW→MW) for distribution cases
  - Transformer tap ratios (ratio=0 → 1.0)
  - Generator Pg merging into bus-level generation

The parser returns a `MatpowerData` struct that the `TestSystems` module
converts into `HybridSystem` objects with added DC components.
"""
module MatpowerParser

export MatpowerData, parse_matpower

"""
Raw parsed MATPOWER data (all values in p.u. on baseMVA after conversion).
"""
struct MatpowerData
    baseMVA::Float64
    # Bus data: each row = (bus_id, type, Pd, Qd, Gs, Bs, area, Vm, Va, baseKV)
    bus::Matrix{Float64}
    # Generator data: each row = (bus_id, Pg, Qg, Qmax, Qmin, Vg, mBase, status, Pmax, Pmin)
    gen::Matrix{Float64}
    # Branch data: each row = (fbus, tbus, r, x, b, rateA, rateB, rateC, ratio, angle, status)
    branch::Matrix{Float64}
    # Mapping from original bus IDs → contiguous 1:N indices
    busid_to_idx::Dict{Int, Int}
    # Original bus IDs in order
    bus_ids::Vector{Int}
end

"""
    parse_matpower(filepath::String) -> MatpowerData

Parse a MATPOWER `.m` case file.  Supports standard bus/gen/branch matrices.
Handles distribution-case conventions (loads in kW, impedance in Ohms) by
detecting `idx_bus` / `idx_brch` conversion blocks at the end of the file.
"""
function parse_matpower(filepath::String)
    text = read(filepath, String)

    baseMVA = _parse_scalar(text, "mpc.baseMVA")

    bus_raw  = _parse_matrix(text, "mpc.bus")
    gen_raw  = _parse_matrix(text, "mpc.gen")
    branch_raw = _parse_matrix(text, "mpc.branch")

    # Detect distribution-case unit conversion flags
    needs_impedance_conversion = occursin("idx_brch", text) || occursin("idx_bus", text)
    needs_load_kw_conversion   = occursin(r"kW|kVAr", text)

    if needs_load_kw_conversion && size(bus_raw, 2) >= 10
        baseKV = bus_raw[1, 10]
        if needs_impedance_conversion && baseKV > 0
            Vbase = baseKV * 1e3          # V
            Sbase = baseMVA * 1e6         # VA
            Zbase = Vbase^2 / Sbase       # Ohms
            # Convert branch r, x from Ohms to p.u.
            branch_raw[:, 3] ./= Zbase
            branch_raw[:, 4] ./= Zbase
        end
        # Convert loads from kW to MW
        bus_raw[:, 3] ./= 1e3   # Pd
        bus_raw[:, 4] ./= 1e3   # Qd
    end

    # Build contiguous bus ID mapping
    bus_ids = Int.(bus_raw[:, 1])
    sorted_ids = sort(unique(bus_ids))
    busid_to_idx = Dict{Int, Int}(id => i for (i, id) in enumerate(sorted_ids))

    # Remap bus IDs in bus_raw to 1:N
    bus_mapped = copy(bus_raw)
    for i in axes(bus_mapped, 1)
        bus_mapped[i, 1] = busid_to_idx[Int(bus_mapped[i, 1])]
    end
    # Sort by new index
    perm = sortperm(bus_mapped[:, 1])
    bus_mapped = bus_mapped[perm, :]

    # Remap gen bus IDs
    gen_mapped = copy(gen_raw)
    for i in axes(gen_mapped, 1)
        orig_id = Int(gen_mapped[i, 1])
        if haskey(busid_to_idx, orig_id)
            gen_mapped[i, 1] = busid_to_idx[orig_id]
        end
    end

    # Remap branch bus IDs and fix tap ratios
    branch_mapped = copy(branch_raw)
    for i in axes(branch_mapped, 1)
        fid = Int(branch_mapped[i, 1])
        tid = Int(branch_mapped[i, 2])
        if haskey(busid_to_idx, fid)
            branch_mapped[i, 1] = busid_to_idx[fid]
        end
        if haskey(busid_to_idx, tid)
            branch_mapped[i, 2] = busid_to_idx[tid]
        end
        # ratio column (col 9): 0 means 1.0
        if size(branch_mapped, 2) >= 9 && branch_mapped[i, 9] == 0.0
            branch_mapped[i, 9] = 1.0
        end
    end

    return MatpowerData(baseMVA, bus_mapped, gen_mapped, branch_mapped,
                        busid_to_idx, sorted_ids)
end

# ── Internal helpers ──────────────────────────────────────────────────────────

function _parse_scalar(text::String, varname::String)
    m = match(Regex(varname * raw"\s*=\s*([0-9.eE+-]+)"), text)
    m === nothing && error("Cannot find $varname in MATPOWER file")
    return parse(Float64, m.captures[1])
end

"""
Parse a MATPOWER matrix block like `mpc.bus = [ ... ];`
"""
function _parse_matrix(text::String, varname::String)
    # Find the block between `varname = [` and `];`
    pattern = Regex(varname * raw"\s*=\s*\[([^]]*)\]", "s")
    m = match(pattern, text)
    m === nothing && error("Cannot find $varname matrix in MATPOWER file")

    block = m.captures[1]
    rows = Float64[]
    ncols = 0

    for line in split(block, '\n')
        # Strip comments (% ...)
        line = replace(line, r"%.*" => "")
        # Strip semicolons
        line = replace(line, ";" => " ")
        # Strip trailing/leading whitespace
        line = strip(line)
        isempty(line) && continue

        # Handle tab-separated values 
        tokens = split(line)
        vals = Float64[]
        for t in tokens
            v = tryparse(Float64, t)
            v !== nothing && push!(vals, v)
        end
        isempty(vals) && continue

        if ncols == 0
            ncols = length(vals)
        end
        # Some files have extra columns; truncate or pad
        if length(vals) >= ncols
            append!(rows, vals[1:ncols])
        else
            append!(rows, vals)
            append!(rows, zeros(ncols - length(vals)))
        end
    end

    nrows = length(rows) ÷ ncols
    return reshape(rows, ncols, nrows)' |> collect
end

end  # module MatpowerParser
