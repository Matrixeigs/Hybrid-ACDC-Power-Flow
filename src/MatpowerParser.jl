"""
    MatpowerParser

Minimal stub for MATPOWER .m file parsing.
"""
module MatpowerParser

export MatpowerData, parse_matpower

struct MatpowerData
    bus::Matrix{Float64}
    gen::Matrix{Float64}
    branch::Matrix{Float64}
    baseMVA::Float64
end

function parse_matpower(filepath::String)
    error("MatpowerParser: parse_matpower is not implemented. File: $filepath")
end

end # module
