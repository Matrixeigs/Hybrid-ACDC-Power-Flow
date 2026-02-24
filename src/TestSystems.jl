"""
    TestSystems.jl

Standard IEEE test cases modified with HVDC links for hybrid AC/DC experiments.
Systems range from small distribution (33-bus) to large transmission (2000-bus):
  IEEE14+DC, IEEE24+DC, IEEE118+DC  (hand-crafted)
  Case33bw+DC, Case33mg+DC, Case69+DC, Case300+DC, Case2000+DC  (parsed from MATPOWER)
"""
module TestSystems

using ..PowerSystem: ACBus, ACBranch, DCBus, DCBranch, VSCConverter, HybridSystem,
                     BusType, PQ, PV, SLACK, ConverterMode, PQ_MODE, VDC_Q, VDC_VAC
using ..MatpowerParser: parse_matpower, MatpowerData

export build_ieee14_acdc, build_ieee24_3area_acdc, build_ieee118_acdc, build_ac_only_version,
       build_case33bw_acdc, build_case33mg_acdc, build_case69_acdc,
       build_case300_acdc, build_case2000_acdc

# ─── IEEE 14-bus + 2 DC buses + 2 VSC converters ─────────────────────────────

"""
    build_ieee14_acdc()

IEEE 14-bus with a 2-terminal HVDC link between bus 2 (Area 1) and bus 9 (Area 2).
Small test system: 14 AC buses, 20 AC branches, 2 DC buses, 1 DC cable, 2 converters.
"""
function build_ieee14_acdc()
    # ── AC buses (IEEE 14 standard data, p.u. on 100 MVA) ──
    ac_buses = [
        ACBus(1,  SLACK, 0.0,    0.0,    2.324, 0.0,   1.060, 0.0, 1),
        ACBus(2,  PV,    0.217,  0.127,  0.40,  0.0,   1.045, 0.0, 1),
        ACBus(3,  PV,    0.942,  0.190,  0.0,   0.0,   1.010, 0.0, 1),
        ACBus(4,  PQ,    0.478,  0.039,  0.0,   0.0,   1.0,   0.0, 1),
        ACBus(5,  PQ,    0.076,  0.016,  0.0,   0.0,   1.0,   0.0, 1),
        ACBus(6,  PV,    0.112,  0.075,  0.0,   0.0,   1.070, 0.0, 2),
        ACBus(7,  PQ,    0.0,    0.0,    0.0,   0.0,   1.0,   0.0, 2),
        ACBus(8,  PV,    0.0,    0.0,    0.0,   0.0,   1.090, 0.0, 2),
        ACBus(9,  PQ,    0.295,  0.166,  0.0,   0.0,   1.0,   0.0, 2),
        ACBus(10, PQ,    0.090,  0.058,  0.0,   0.0,   1.0,   0.0, 2),
        ACBus(11, PQ,    0.035,  0.018,  0.0,   0.0,   1.0,   0.0, 2),
        ACBus(12, PQ,    0.061,  0.016,  0.0,   0.0,   1.0,   0.0, 2),
        ACBus(13, PQ,    0.135,  0.058,  0.0,   0.0,   1.0,   0.0, 2),
        ACBus(14, PQ,    0.149,  0.050,  0.0,   0.0,   1.0,   0.0, 2),
    ]

    # ── AC branches (IEEE 14 standard) ──
    ac_branches = [
        ACBranch(1,  2,  0.01938, 0.05917, 0.0528, 1.0, true),
        ACBranch(1,  5,  0.05403, 0.22304, 0.0492, 1.0, true),
        ACBranch(2,  3,  0.04699, 0.19797, 0.0438, 1.0, true),
        ACBranch(2,  4,  0.05811, 0.17632, 0.0374, 1.0, true),
        ACBranch(2,  5,  0.05695, 0.17388, 0.0340, 1.0, true),
        ACBranch(3,  4,  0.06701, 0.17103, 0.0346, 1.0, true),
        ACBranch(4,  5,  0.01335, 0.04211, 0.0128, 1.0, true),
        ACBranch(4,  7,  0.0,     0.20912, 0.0,    0.978, true),  # transformer
        ACBranch(4,  9,  0.0,     0.55618, 0.0,    0.969, true),  # transformer
        ACBranch(5,  6,  0.0,     0.25202, 0.0,    0.932, true),  # transformer
        ACBranch(6,  11, 0.09498, 0.19890, 0.0,    1.0, true),
        ACBranch(6,  12, 0.12291, 0.25581, 0.0,    1.0, true),
        ACBranch(6,  13, 0.06615, 0.13027, 0.0,    1.0, true),
        ACBranch(7,  8,  0.0,     0.17615, 0.0,    1.0, true),
        ACBranch(7,  9,  0.11001, 0.20640, 0.0,    1.0, true),
        ACBranch(9,  10, 0.03181, 0.08450, 0.0,    1.0, true),
        ACBranch(9,  14, 0.12711, 0.27038, 0.0,    1.0, true),
        ACBranch(10, 11, 0.08205, 0.19207, 0.0,    1.0, true),
        ACBranch(12, 13, 0.22092, 0.19988, 0.0,    1.0, true),
        ACBranch(13, 14, 0.17093, 0.34802, 0.0,    1.0, true),
    ]

    # ── DC buses (2-terminal HVDC) ──
    dc_buses = [
        DCBus(1, 1.0, 0.0),   # DC bus at converter 1
        DCBus(2, 1.0, 0.0),   # DC bus at converter 2
    ]

    # ── DC cable ──
    dc_branches = [
        DCBranch(1, 2, 0.01, true),   # Low resistance DC cable
    ]

    # ── VSC converters ──
    converters = [
        VSCConverter(1, 2, 1, PQ_MODE,  0.15, 0.0, 1.0, 1.045,
                     0.001, 0.01, 0.001, 1.0, true),  # Bus 2 → DC bus 1
        VSCConverter(2, 9, 2, VDC_Q,    0.0,  0.0, 1.0, 1.0,
                     0.001, 0.01, 0.001, 1.0, true),   # Bus 9 → DC bus 2
    ]

    return HybridSystem(ac_buses, ac_branches, dc_buses, dc_branches, converters)
end

# ─── IEEE 24-bus RTS, 3-area hybrid AC/DC ────────────────────────────────────

"""
    build_ieee24_3area_acdc()

Modified IEEE 24-bus RTS as a 3-area system connected by HVDC links.
- Area 1: buses 1-8  (generation-heavy)
- Area 2: buses 9-16  (mixed)
- Area 3: buses 17-24 (load-heavy)
DC interconnections: Area 1↔Area 2, Area 2↔Area 3.
"""
function build_ieee24_3area_acdc()
    # ── AC buses (simplified IEEE 24 RTS data) ──
    ac_buses = [
        # Area 1 (buses 1-8)
        ACBus(1,  PV,    1.08, 0.22, 1.92, 0.0, 1.035, 0.0, 1),
        ACBus(2,  PV,    0.97, 0.20, 1.92, 0.0, 1.035, 0.0, 1),
        ACBus(3,  PQ,    1.80, 0.37, 0.0,  0.0, 1.0,   0.0, 1),
        ACBus(4,  PQ,    0.74, 0.15, 0.0,  0.0, 1.0,   0.0, 1),
        ACBus(5,  PQ,    0.71, 0.14, 0.0,  0.0, 1.0,   0.0, 1),
        ACBus(6,  PQ,    1.36, 0.28, 0.0,  0.0, 1.0,   0.0, 1),
        ACBus(7,  PV,    1.25, 0.25, 2.40, 0.0, 1.025, 0.0, 1),
        ACBus(8,  PQ,    1.71, 0.35, 0.0,  0.0, 1.0,   0.0, 1),
        # Area 2 (buses 9-16)
        ACBus(9,  PQ,    1.75, 0.36, 0.0,  0.0, 1.0,   0.0, 2),
        ACBus(10, PQ,    1.95, 0.40, 0.0,  0.0, 1.0,   0.0, 2),
        ACBus(11, PQ,    0.0,  0.0,  0.0,  0.0, 1.0,   0.0, 2),
        ACBus(12, PQ,    0.0,  0.0,  0.0,  0.0, 1.0,   0.0, 2),
        ACBus(13, SLACK, 2.65, 0.54, 5.91, 0.0, 1.020, 0.0, 2),
        ACBus(14, PV,    1.94, 0.39, 0.0,  0.0, 1.0,   0.0, 2),
        ACBus(15, PV,    3.17, 0.64, 2.15, 0.0, 1.014, 0.0, 2),
        ACBus(16, PV,    1.00, 0.20, 1.55, 0.0, 1.017, 0.0, 2),
        # Area 3 (buses 17-24)
        ACBus(17, PQ,    0.0,  0.0,  0.0,  0.0, 1.0,   0.0, 3),
        ACBus(18, PV,    3.33, 0.68, 4.00, 0.0, 1.050, 0.0, 3),
        ACBus(19, PQ,    1.81, 0.37, 0.0,  0.0, 1.0,   0.0, 3),
        ACBus(20, PQ,    1.28, 0.26, 0.0,  0.0, 1.0,   0.0, 3),
        ACBus(21, PV,    0.0,  0.0,  4.00, 0.0, 1.050, 0.0, 3),
        ACBus(22, PV,    0.0,  0.0,  3.00, 0.0, 1.050, 0.0, 3),
        ACBus(23, PV,    0.0,  0.0,  6.60, 0.0, 1.050, 0.0, 3),
        ACBus(24, PQ,    0.0,  0.0,  0.0,  0.0, 1.0,   0.0, 3),
    ]

    # ── AC branches (simplified IEEE 24 RTS topology) ──
    ac_branches = [
        # Area 1 internal
        ACBranch(1, 2,  0.0026, 0.0139, 0.4611, 1.0, true),
        ACBranch(1, 3,  0.0546, 0.2112, 0.0572, 1.0, true),
        ACBranch(1, 5,  0.0218, 0.0845, 0.0229, 1.0, true),
        ACBranch(2, 4,  0.0328, 0.1267, 0.0343, 1.0, true),
        ACBranch(2, 6,  0.0497, 0.1920, 0.0520, 1.0, true),
        ACBranch(3, 9,  0.0308, 0.1190, 0.0322, 1.0, true),  # inter-area AC tie
        ACBranch(4, 9,  0.0268, 0.1037, 0.0281, 1.0, true),  # inter-area AC tie
        ACBranch(5, 10, 0.0228, 0.0883, 0.0239, 1.0, true),  # inter-area AC tie
        ACBranch(6, 10, 0.0139, 0.0605, 0.2459, 1.0, true),  # inter-area AC tie
        ACBranch(7, 8,  0.0159, 0.0614, 0.0166, 1.0, true),
        ACBranch(3, 24, 0.0023, 0.0839, 0.0,    1.015, true), # transformer
        ACBranch(8, 10, 0.0427, 0.1651, 0.0447, 1.0, true),
        # Area 2 internal
        ACBranch(9, 11,  0.0023, 0.0839, 0.0, 1.03, true),  # transformer
        ACBranch(9, 12,  0.0023, 0.0839, 0.0, 1.03, true),  # transformer
        ACBranch(10, 11, 0.0023, 0.0839, 0.0, 1.02, true),  # transformer
        ACBranch(10, 12, 0.0023, 0.0839, 0.0, 1.02, true),  # transformer
        ACBranch(11, 13, 0.0061, 0.0476, 0.0999, 1.0, true),
        ACBranch(11, 14, 0.0054, 0.0418, 0.0879, 1.0, true),
        ACBranch(12, 13, 0.0061, 0.0476, 0.0999, 1.0, true),
        ACBranch(12, 23, 0.0124, 0.0966, 0.2030, 1.0, true),
        ACBranch(13, 23, 0.0111, 0.0865, 0.1818, 1.0, true),
        # Area 3 internal
        ACBranch(14, 16, 0.0050, 0.0389, 0.0818, 1.0, true),
        ACBranch(15, 16, 0.0022, 0.0173, 0.0364, 1.0, true),
        ACBranch(15, 21, 0.0063, 0.0490, 0.1030, 1.0, true),
        ACBranch(15, 24, 0.0067, 0.0519, 0.1091, 1.0, true),
        ACBranch(16, 17, 0.0033, 0.0259, 0.0545, 1.0, true),
        ACBranch(16, 19, 0.0030, 0.0231, 0.0485, 1.0, true),
        ACBranch(17, 18, 0.0018, 0.0144, 0.0303, 1.0, true),
        ACBranch(17, 22, 0.0135, 0.1053, 0.2212, 1.0, true),
        ACBranch(18, 21, 0.0033, 0.0259, 0.0545, 1.0, true),
        ACBranch(19, 20, 0.0025, 0.0198, 0.0417, 1.0, true),
        ACBranch(20, 23, 0.0022, 0.0173, 0.0364, 1.0, true),
        ACBranch(21, 22, 0.0087, 0.0678, 0.1424, 1.0, true),
    ]

    # ── DC network (4 DC buses, 3 DC cables) ──
    dc_buses = [
        DCBus(1, 1.0, 0.0),   # near AC bus 7 (Area 1)
        DCBus(2, 1.0, 0.0),   # near AC bus 13 (Area 2)
        DCBus(3, 1.0, 0.0),   # near AC bus 15 (Area 2)
        DCBus(4, 1.0, 0.0),   # near AC bus 21 (Area 3)
    ]

    dc_branches = [
        DCBranch(1, 2, 0.005, true),   # Area 1 ↔ Area 2
        DCBranch(3, 4, 0.005, true),   # Area 2 ↔ Area 3
        DCBranch(2, 3, 0.003, true),   # DC backbone
    ]

    # ── VSC converters ──
    converters = [
        VSCConverter(1, 7,  1, PQ_MODE, 0.20, 0.0, 1.0, 1.025,
                     0.001, 0.01, 0.001, 2.0, true),    # Area 1
        VSCConverter(2, 13, 2, VDC_Q,   0.0,  0.0, 1.0, 1.020,
                     0.001, 0.01, 0.001, 2.0, true),    # Area 2 (slack side)
        VSCConverter(3, 15, 3, PQ_MODE, 0.15, 0.0, 1.0, 1.014,
                     0.001, 0.01, 0.001, 2.0, true),    # Area 2
        VSCConverter(4, 21, 4, VDC_Q,   0.0,  0.0, 1.0, 1.050,
                     0.001, 0.01, 0.001, 2.0, true),    # Area 3
    ]

    return HybridSystem(ac_buses, ac_branches, dc_buses, dc_branches, converters)
end

# ─── IEEE 118-bus hybrid AC/DC (large) ───────────────────────────────────────

const _CASE118_LOADED = Ref(false)

module _Case118Shim
    # Constants expected by PF/case118_hybrid_acdc.jl
    const PQ_BUS = 1
    const PV_BUS = 2
    const SLACK_BUS = 3

    const DC_NORMAL = 1
    const DC_SLACK = 2

    const PQ_MODE = 1
    const P_VDC_MODE = 2
    const DROOP_MODE = 3

    struct ACBus
        bus_type::Int
        Pd::Float64
        Qd::Float64
        Pg::Float64
        Qg::Float64
        Vset::Float64
        Qmin::Float64
        Qmax::Float64
    end

    struct ACBranch
        from::Int
        to::Int
        r::Float64
        x::Float64
        b::Float64
        tap::Float64
        shift::Float64
    end

    ACBranch(from::Int, to::Int, r::Real, x::Real, b::Real) =
        ACBranch(from, to, Float64(r), Float64(x), Float64(b), 1.0, 0.0)
    ACBranch(from::Int, to::Int, r::Real, x::Real, b::Real, tap::Real, shift::Real) =
        ACBranch(from, to, Float64(r), Float64(x), Float64(b), Float64(tap), Float64(shift))

    struct ACNetwork
        nb::Int
        nl::Int
        buses::Vector{ACBus}
        branches::Vector{ACBranch}
        baseMVA::Float64
    end

    struct DCBus
        dc_type::Int
        Pdc::Float64
        Vdc::Float64
    end

    struct DCBranch
        from::Int
        to::Int
        r::Float64
    end

    struct DCNetwork
        nb::Int
        nl::Int
        buses::Vector{DCBus}
        branches::Vector{DCBranch}
    end

    struct QuadraticLoss
        a::Float64
        b::Float64
        c::Float64
    end

    struct VSCConverter
        ac_bus::Int
        dc_bus::Int
        mode::Int
        P_sp::Float64
        Q_sp::Float64
        Vdc_sp::Float64
        loss_model::QuadraticLoss
    end

    function VSCConverter(ac_bus::Int, dc_bus::Int; mode::Int,
                          P_sp::Real=0.0, Q_sp::Real=0.0, Vdc_sp::Real=1.0,
                          loss_model::QuadraticLoss=QuadraticLoss(0.0, 0.0, 0.0),
                          kwargs...)
        return VSCConverter(ac_bus, dc_bus, mode, Float64(P_sp), Float64(Q_sp),
                            Float64(Vdc_sp), loss_model)
    end

    struct HybridNetwork
        ac_net::ACNetwork
        dc_net::DCNetwork
        converters::Vector{VSCConverter}
    end
end

"""
    build_ieee118_acdc()

IEEE 118-bus hybrid AC/DC system (large case) based on `PF/case118_hybrid_acdc.jl`:
- 118 AC buses, 186 AC branches
- 8 DC buses, 9 DC branches
- 6 VSC converters (two MTDC corridors)

Notes:
  - The original case file provides no official "area" partition. We assign a
    simple 4-area split by bus index ranges: 1--30, 31--60, 61--90, 91--118.
  - DC slack is handled by reindexing the DC network so the slack-type DC bus
    becomes DC bus 1 (matching the NR solver's convention).
"""
function build_ieee118_acdc()
    if !_CASE118_LOADED[]
        include(joinpath(@__DIR__, "..", "..", "PF", "case118_hybrid_acdc.jl"))
        _CASE118_LOADED[] = true
    end

    net = Base.invokelatest(create_case118_hybrid, _Case118Shim)

    # AC buses
    ac_buses = Vector{ACBus}(undef, length(net.ac_net.buses))
    for (i, b) in enumerate(net.ac_net.buses)
        bus_type = if b.bus_type == _Case118Shim.PQ_BUS
            PQ
        elseif b.bus_type == _Case118Shim.PV_BUS
            PV
        else
            SLACK
        end

        area = if i <= 30
            1
        elseif i <= 60
            2
        elseif i <= 90
            3
        else
            4
        end

        ac_buses[i] = ACBus(i, bus_type, b.Pd, b.Qd, b.Pg, b.Qg, b.Vset, 0.0, area)
    end

    # AC branches
    ac_branches = Vector{ACBranch}(undef, length(net.ac_net.branches))
    for (k, br) in enumerate(net.ac_net.branches)
        tap = br.tap == 0.0 ? 1.0 : br.tap
        ac_branches[k] = ACBranch(br.from, br.to, br.r, br.x, br.b, tap, true)
    end

    # DC buses: reorder so slack bus becomes index 1 (solver convention)
    ndc = length(net.dc_net.buses)
    slack_old = findfirst(b -> b.dc_type == _Case118Shim.DC_SLACK, net.dc_net.buses)
    slack_old === nothing && (slack_old = 1)
    perm = vcat([slack_old], [i for i in 1:ndc if i != slack_old])
    old2new = Dict(old => new for (new, old) in enumerate(perm))

    dc_buses = Vector{DCBus}(undef, ndc)
    for (new_i, old_i) in enumerate(perm)
        b = net.dc_net.buses[old_i]
        dc_buses[new_i] = DCBus(new_i, b.Vdc, b.Pdc)
    end

    dc_branches = DCBranch[]
    for br in net.dc_net.branches
        from_new = old2new[br.from]
        to_new = old2new[br.to]
        push!(dc_branches, DCBranch(from_new, to_new, br.r, true))
    end

    # Converters
    converters = VSCConverter[]
    for (cid, c) in enumerate(net.converters)
        dc_new = old2new[c.dc_bus]

        mode = if c.mode == _Case118Shim.PQ_MODE
            PQ_MODE
        else
            VDC_Q
        end

        # Map quadratic loss
        la = c.loss_model.a
        lb = c.loss_model.b
        lc = c.loss_model.c

        # Use AC bus setpoint as Vac_set (not enforced in simplified model)
        Vac_set = ac_buses[c.ac_bus].Vm

        push!(converters, VSCConverter(
            cid, c.ac_bus, dc_new, mode,
            c.P_sp, c.Q_sp, c.Vdc_sp, Vac_set,
            la, lb, lc, 5.0, true
        ))
    end

    return HybridSystem(ac_buses, ac_branches, dc_buses, dc_branches, converters)
end

# ─── AC-only version for comparison ──────────────────────────────────────────

"""
Strip DC components for AC-only baseline comparison.
"""
function build_ac_only_version(sys::HybridSystem)
    return HybridSystem(sys.ac_buses, sys.ac_branches,
                        DCBus[], DCBranch[], VSCConverter[];
                        baseMVA=sys.baseMVA)
end

# ═══════════════════════════════════════════════════════════════════════════════
#  MATPOWER-based hybrid AC/DC test systems
# ═══════════════════════════════════════════════════════════════════════════════

# ── Helper: convert MatpowerData → vectors of ACBus / ACBranch ────────────────

const _DATA_DIR = joinpath(@__DIR__, "..", "data")

"""
    _matpower_to_ac(mpd::MatpowerData) -> (Vector{ACBus}, Vector{ACBranch}, Float64)

Convert parsed MATPOWER data into ACBus/ACBranch vectors, merging generator
data into bus records.  Returns (buses, branches, baseMVA).
"""
function _matpower_to_ac(mpd::MatpowerData)
    nbus = size(mpd.bus, 1)
    nbranch = size(mpd.branch, 1)
    ngen = size(mpd.gen, 1)

    # Aggregate generation per bus: Pg, Qg, Vg (keep the last Vg encountered)
    bus_Pg = zeros(nbus)
    bus_Qg = zeros(nbus)
    bus_Vg = ones(nbus)
    bus_has_gen = falses(nbus)

    for g in 1:ngen
        # gen columns: bus, Pg, Qg, Qmax, Qmin, Vg, mBase, status, Pmax, Pmin
        status = size(mpd.gen, 2) >= 8 ? mpd.gen[g, 8] : 1.0
        status ≤ 0 && continue
        bi = Int(mpd.gen[g, 1])
        (bi < 1 || bi > nbus) && continue
        bus_Pg[bi] += mpd.gen[g, 2] / mpd.baseMVA   # MW → p.u.
        bus_Qg[bi] += mpd.gen[g, 3] / mpd.baseMVA
        bus_Vg[bi]  = mpd.gen[g, 6]
        bus_has_gen[bi] = true
    end

    ac_buses = Vector{ACBus}(undef, nbus)
    for i in 1:nbus
        # bus columns: id, type, Pd, Qd, Gs, Bs, area, Vm, Va, baseKV, ...
        orig_type = Int(mpd.bus[i, 2])
        btype = orig_type == 3 ? SLACK : (orig_type == 2 ? PV : PQ)
        Pd = mpd.bus[i, 3] / mpd.baseMVA   # MW → p.u.
        Qd = mpd.bus[i, 4] / mpd.baseMVA
        area = size(mpd.bus, 2) >= 7 ? Int(mpd.bus[i, 7]) : 1
        Vm = bus_has_gen[i] ? bus_Vg[i] : mpd.bus[i, 8]
        Va = size(mpd.bus, 2) >= 9 ? deg2rad(mpd.bus[i, 9]) : 0.0

        ac_buses[i] = ACBus(i, btype, Pd, Qd, bus_Pg[i], bus_Qg[i], Vm, Va, area)
    end

    ac_branches = ACBranch[]
    for k in 1:nbranch
        # branch columns: fbus, tbus, r, x, b, rateA, rateB, rateC, ratio, angle, status
        from = Int(mpd.branch[k, 1])
        to   = Int(mpd.branch[k, 2])
        r    = mpd.branch[k, 3]
        x    = mpd.branch[k, 4]
        b_ch = mpd.branch[k, 5]
        tap  = size(mpd.branch, 2) >= 9 ? mpd.branch[k, 9] : 1.0
        tap == 0.0 && (tap = 1.0)
        br_status = size(mpd.branch, 2) >= 11 ? (mpd.branch[k, 11] > 0) : true
        push!(ac_branches, ACBranch(from, to, r, x, b_ch, tap, br_status))
    end

    return ac_buses, ac_branches, mpd.baseMVA
end

# ─── 33-bus Baran & Wu (distribution) + HVDC/microgrid DC link ───────────────

"""
    build_case33bw_acdc()

33-bus distribution system (Baran & Wu, 1989) with a 3-terminal DC microgrid.
- 33 AC buses, 32 AC branches (radial, tie switches excluded)
- 3 DC buses, 2 DC branches
- 3 VSC converters connecting DC microgrids at buses 6, 18, 33

The DC microgrid represents a local area network for distributed generation
(PV + battery storage) typical of modern distribution systems.
"""
function build_case33bw_acdc()
    mpd = parse_matpower(joinpath(_DATA_DIR, "case33bw.m"))
    ac_buses, ac_branches_all, baseMVA = _matpower_to_ac(mpd)

    # Keep only active branches (status=true); the case33bw file has 5 normally-open
    # tie switches (status=0) at the end → filter them out for radial topology
    ac_branches = filter(b -> b.status, ac_branches_all)

    # DC microgrid: 3 nodes connected to ends of main feeders
    dc_buses = [
        DCBus(1, 1.0, 0.0),   # near bus 6  (main feeder junction)
        DCBus(2, 1.0, 0.0),   # near bus 18 (end of feeder 1)
        DCBus(3, 1.0, 0.0),   # near bus 33 (end of feeder 2)
    ]
    dc_branches = [
        DCBranch(1, 2, 0.02, true),
        DCBranch(1, 3, 0.02, true),
    ]

    converters = [
        VSCConverter(1, 6,  1, PQ_MODE, 0.05, 0.0, 1.0, 1.0,
                     0.001, 0.01, 0.001, 0.5, true),
        VSCConverter(2, 18, 2, VDC_Q,   0.0,  0.0, 1.0, 1.0,
                     0.001, 0.01, 0.001, 0.5, true),
        VSCConverter(3, 33, 3, PQ_MODE, 0.03, 0.0, 1.0, 1.0,
                     0.001, 0.01, 0.001, 0.5, true),
    ]

    return HybridSystem(ac_buses, ac_branches, dc_buses, dc_branches, converters;
                        baseMVA=baseMVA)
end

# ─── 33-bus (Kashem microgrid variant) + HVDC ────────────────────────────────

"""
    build_case33mg_acdc()

33-bus distribution system (Kashem microgrid variant) with a 3-terminal DC link.
Same topology as case33bw but with different branch 7 impedance and baseMVA=1.
- 33 AC buses, 32 AC branches (radial)
- 3 DC buses, 2 DC branches
- 3 VSC converters at buses 6, 18, 33
"""
function build_case33mg_acdc()
    mpd = parse_matpower(joinpath(_DATA_DIR, "case33mg.m"))
    ac_buses, ac_branches_all, baseMVA = _matpower_to_ac(mpd)
    ac_branches = filter(b -> b.status, ac_branches_all)

    dc_buses = [
        DCBus(1, 1.0, 0.0),
        DCBus(2, 1.0, 0.0),
        DCBus(3, 1.0, 0.0),
    ]
    dc_branches = [
        DCBranch(1, 2, 0.02, true),
        DCBranch(1, 3, 0.02, true),
    ]
    converters = [
        VSCConverter(1, 6,  1, PQ_MODE, 0.05, 0.0, 1.0, 1.0,
                     0.001, 0.01, 0.001, 0.5, true),
        VSCConverter(2, 18, 2, VDC_Q,   0.0,  0.0, 1.0, 1.0,
                     0.001, 0.01, 0.001, 0.5, true),
        VSCConverter(3, 33, 3, PQ_MODE, 0.03, 0.0, 1.0, 1.0,
                     0.001, 0.01, 0.001, 0.5, true),
    ]

    return HybridSystem(ac_buses, ac_branches, dc_buses, dc_branches, converters;
                        baseMVA=baseMVA)
end

# ─── 69-bus distribution system + HVDC ───────────────────────────────────────

"""
    build_case69_acdc()

69-bus PG&E distribution system (Baran & Wu, 1989) with a 4-terminal DC network.
- 69 AC buses, 69 AC branches
- 4 DC buses, 3 DC branches
- 4 VSC converters at load centers: buses 11, 27, 50, 61

The DC backbone interconnects the four main lateral feeders, enabling
power sharing between lightly-loaded and heavily-loaded segments.
"""
function build_case69_acdc()
    mpd = parse_matpower(joinpath(_DATA_DIR, "case69.m"))
    ac_buses, ac_branches, baseMVA = _matpower_to_ac(mpd)

    # DC backbone across four feeders
    dc_buses = [
        DCBus(1, 1.0, 0.0),   # near bus 11 (end of lateral 1)
        DCBus(2, 1.0, 0.0),   # near bus 27 (end of lateral 2)
        DCBus(3, 1.0, 0.0),   # near bus 50 (heavy-load lateral 3)
        DCBus(4, 1.0, 0.0),   # near bus 61 (heavy-load lateral 4)
    ]
    dc_branches = [
        DCBranch(1, 2, 0.015, true),
        DCBranch(2, 3, 0.015, true),
        DCBranch(3, 4, 0.015, true),
    ]

    converters = [
        VSCConverter(1, 11, 1, PQ_MODE, 0.04,  0.0, 1.0, 1.0,
                     0.001, 0.01, 0.001, 0.5, true),
        VSCConverter(2, 27, 2, VDC_Q,   0.0,   0.0, 1.0, 1.0,
                     0.001, 0.01, 0.001, 0.5, true),
        VSCConverter(3, 50, 3, PQ_MODE, 0.10,  0.0, 1.0, 1.0,
                     0.001, 0.01, 0.001, 1.0, true),
        VSCConverter(4, 61, 4, PQ_MODE, 0.20,  0.0, 1.0, 1.0,
                     0.001, 0.01, 0.001, 2.0, true),
    ]

    return HybridSystem(ac_buses, ac_branches, dc_buses, dc_branches, converters;
                        baseMVA=baseMVA)
end

# ─── IEEE 300-bus transmission system + MTDC overlay ─────────────────────────

"""
    build_case300_acdc()

IEEE 300-bus transmission system with a 6-terminal MTDC overlay.
- 300 AC buses, 411 AC branches, 69 generators
- 6 DC buses, 6 DC branches (meshed ring)
- 6 VSC converters at generation/load hubs

The MTDC grid connects key PV/generation buses across the network,
providing corridor reinforcement for inter-area power transfers.

Bus mapping: the original IEEE 300 case uses non-sequential bus IDs
(1–9533). These are mapped to contiguous 1:300 indices internally.
"""
function build_case300_acdc()
    mpd = parse_matpower(joinpath(_DATA_DIR, "case300.m"))
    ac_buses, ac_branches, baseMVA = _matpower_to_ac(mpd)

    # Identify key buses for converter placement by looking at PV generators
    # We pick 6 well-distributed large PV/generator buses from different parts
    # of the network. Original bus IDs → mapped sequential IDs via mpd.busid_to_idx
    # Slack bus (7049) maps to one index; pick major PV buses around the network.
    # We'll use the mapped indices for buses that correspond to
    # original IDs: 8 (gen), 76 (gen), 119 (gen), 149 (gen), 198 (gen), 243 (gen)
    conv_ac_buses = [
        mpd.busid_to_idx[8],    # generator bus, area 1
        mpd.busid_to_idx[76],   # generator bus
        mpd.busid_to_idx[119],  # generator bus
        mpd.busid_to_idx[149],  # generator bus
        mpd.busid_to_idx[198],  # generator bus
        mpd.busid_to_idx[243],  # generator bus
    ]

    dc_buses = [DCBus(i, 1.0, 0.0) for i in 1:6]

    # Meshed DC ring with cross-links
    dc_branches = [
        DCBranch(1, 2, 0.005, true),
        DCBranch(2, 3, 0.005, true),
        DCBranch(3, 4, 0.005, true),
        DCBranch(4, 5, 0.005, true),
        DCBranch(5, 6, 0.005, true),
        DCBranch(6, 1, 0.005, true),   # close the ring
    ]

    converters = [
        VSCConverter(1, conv_ac_buses[1], 1, PQ_MODE, 0.50, 0.0, 1.0,
                     ac_buses[conv_ac_buses[1]].Vm,
                     0.001, 0.01, 0.001, 3.0, true),
        VSCConverter(2, conv_ac_buses[2], 2, VDC_Q,   0.0,  0.0, 1.0,
                     ac_buses[conv_ac_buses[2]].Vm,
                     0.001, 0.01, 0.001, 3.0, true),
        VSCConverter(3, conv_ac_buses[3], 3, PQ_MODE, 0.30, 0.0, 1.0,
                     ac_buses[conv_ac_buses[3]].Vm,
                     0.001, 0.01, 0.001, 3.0, true),
        VSCConverter(4, conv_ac_buses[4], 4, PQ_MODE, 0.40, 0.0, 1.0,
                     ac_buses[conv_ac_buses[4]].Vm,
                     0.001, 0.01, 0.001, 3.0, true),
        VSCConverter(5, conv_ac_buses[5], 5, VDC_Q,   0.0,  0.0, 1.0,
                     ac_buses[conv_ac_buses[5]].Vm,
                     0.001, 0.01, 0.001, 3.0, true),
        VSCConverter(6, conv_ac_buses[6], 6, PQ_MODE, 0.35, 0.0, 1.0,
                     ac_buses[conv_ac_buses[6]].Vm,
                     0.001, 0.01, 0.001, 3.0, true),
    ]

    return HybridSystem(ac_buses, ac_branches, dc_buses, dc_branches, converters;
                        baseMVA=baseMVA)
end

# ─── ACTIVSg 2000-bus Texas synthetic + MTDC backbone ────────────────────────

"""
    build_case2000_acdc()

ACTIVSg 2000-bus synthetic Texas system with an 8-terminal MTDC backbone.
- 2000 AC buses, ~3206 AC branches, 544 generators, 8 areas
- 8 DC buses (one per area), 9 DC branches (inter-area ring + cross-links)
- 8 VSC converters (one per area at the largest PV/generator bus)

The MTDC backbone mirrors the 8-area structure, providing DC corridors
for inter-area power transfer reinforcement.  Bus IDs 1001–8160 are
mapped to contiguous 1:2000 indices internally.
"""
function build_case2000_acdc()
    mpd = parse_matpower(joinpath(_DATA_DIR, "case_ACTIVSg2000.m"))
    ac_buses, ac_branches, baseMVA = _matpower_to_ac(mpd)

    nbus = length(ac_buses)

    # Find the largest PV/generator bus in each area (by net Pg)
    area_best = Dict{Int, Tuple{Int, Float64}}()  # area → (bus_idx, Pg)
    for b in ac_buses
        if (b.type == PV || b.type == SLACK) && b.Pg > 0
            if !haskey(area_best, b.area) || b.Pg > area_best[b.area][2]
                area_best[b.area] = (b.id, b.Pg)
            end
        end
    end

    # Sort areas 1..8 and collect converter AC buses
    areas_sorted = sort(collect(keys(area_best)))
    n_areas = length(areas_sorted)
    conv_ac_buses = [area_best[a][1] for a in areas_sorted]

    # DC network: one DC bus per area
    dc_buses = [DCBus(i, 1.0, 0.0) for i in 1:n_areas]

    # DC branches: ring connecting adjacent areas + two cross-links
    dc_branches = DCBranch[]
    for i in 1:n_areas-1
        push!(dc_branches, DCBranch(i, i+1, 0.003, true))
    end
    push!(dc_branches, DCBranch(n_areas, 1, 0.003, true))    # close ring
    # Additional cross-links for redundancy (if ≥4 areas)
    if n_areas >= 4
        push!(dc_branches, DCBranch(1, div(n_areas, 2) + 1, 0.005, true))
    end

    # Converters: alternate PQ_MODE and VDC_Q so DC network has both setpoint
    # buses and voltage-regulating buses
    converters = VSCConverter[]
    for (i, ab) in enumerate(conv_ac_buses)
        mode = iseven(i) ? VDC_Q : PQ_MODE
        Pset = mode == PQ_MODE ? 0.50 : 0.0
        push!(converters, VSCConverter(
            i, ab, i, mode, Pset, 0.0, 1.0,
            ac_buses[ab].Vm,
            0.001, 0.01, 0.001, 5.0, true
        ))
    end

    return HybridSystem(ac_buses, ac_branches, dc_buses, dc_branches, converters;
                        baseMVA=baseMVA)
end

end  # module TestSystems
