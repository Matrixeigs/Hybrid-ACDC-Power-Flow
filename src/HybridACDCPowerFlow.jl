"""
    HybridACDCPowerFlow

A Julia module for hybrid AC/DC power flow analysis supporting:
- AC power flow with Newton-Raphson solver
- DC network with linear analysis
- VSC converters with multiple control modes (PQ, VDC_Q, VDC_VAC)
- IEEE test systems with HVDC extensions
- **ENHANCED v0.2.0**: Island detection, PV→PQ conversion, auto swing bus selection, converter mode switching
- **NEW v0.3.0**: Distributed slack bus model for realistic multi-generator control
- **NEW v0.4.0**: Optimized sparse Jacobian, pre-allocated workspace, optional feasibility extension
- **NEW v0.5.0**: Full JuliaPowerCase integration (shared types, keyword constructors)

Performance optimizations:
- Sparse Jacobian with pre-computed sparsity pattern
- Sparse LU factorization with symbolic reuse (UMFPACK)
- Pre-allocated SolverWorkspace with pre-computed entry map
- @inbounds and sin/cos caching for inner loops
- JuMP/Ipopt/NLsolve moved to optional extension

JuliaPowerCase integration (v0.5.0):
- All types (ACBus, ACBranch, DCBus, DCBranch, VSCConverter, Generator)
  imported directly from JuliaPowerCase
- Converter control modes use Symbol (:pq, :vdc_q, :vdc_vac)
- HybridSystem stores Generator objects separately from buses

Author: Tianyang Zhao
Version: 0.5.0 (JuliaPowerCase Integration)
"""
module HybridACDCPowerFlow

# Re-export everything from PowerSystem, MatpowerParser, and TestSystems
include("PowerSystem.jl")
include("MatpowerParser.jl")
include("TestSystems.jl")
include("PowerSystemEnhanced.jl")

using .PowerSystem
using .MatpowerParser
using .TestSystems
using .PowerSystemEnhanced

# Import detect_islands explicitly from PowerSystemEnhanced to resolve conflict
import .PowerSystemEnhanced: detect_islands, extract_island_subsystem

# Include JuliaPowerCase adapter for seamless data integration
include("JuliaPowerCaseAdapter.jl")

# Include feasibility checker stubs (implementations in extension)
include("FeasibilityChecker.jl")

# Re-export all public functions from base modules
export ACBus, ACBranch, DCBus, DCBranch, VSCConverter, HybridSystem, Generator,
       # JuliaPowerCase types (recommended for new code)
       HybridPowerSystem, IslandInfo,
       BusType, PQ, PV, SLACK, PQ_MODE, VDC_Q, VDC_VAC,
       build_admittance_matrix, solve_power_flow, power_flow_residual,
       get_bus_voltages, get_branch_flows, rebuild_matrices!,
       remove_ac_branch, extract_graph_data,
       build_ieee14_acdc, build_ieee24_3area_acdc, build_ieee118_acdc, build_ac_only_version,
       build_case33bw_acdc, build_case33mg_acdc, build_case69_acdc,
       build_case300_acdc, build_case2000_acdc,
       # Conversion functions
       to_hybrid_system, update_results!, to_solver_system,
       # Optimization exports
       SolverWorkspace, create_solver_workspace, build_jacobian_triplets!,
       compute_power_injections!

# Re-export enhanced functions
export detect_islands, solve_power_flow_islanded, solve_power_flow_adaptive,
       check_reactive_limits, pv_to_pq_conversion!, auto_select_swing_bus,
       auto_switch_converter_mode!, PowerFlowOptions, ReactiveLimit,
       create_default_Q_limits, print_island_summary, extract_island_subsystem,
       DistributedSlack, create_participation_factors, solve_power_flow_distributed_slack,
       solve_power_flow_distributed_slack_full,
       # Feasibility (requires extension)
       check_power_flow_feasibility, FeasibilityResult,
       check_power_flow_feasibility_nlsolve, check_power_flow_feasibility_jump,
       validate_against_nlsolve, validate_against_jump

end # module
