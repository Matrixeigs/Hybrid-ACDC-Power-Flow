"""
    HybridACDCPowerFlow

A Julia module for hybrid AC/DC power flow analysis.

# Architecture (v0.7.0 - Pure Application Layer)
This module is a pure **application layer** that operates on data structures 
from **JuliaPowerCase** (the information layer).

 ## Key Features
- AC power flow with Newton-Raphson solver
- DC network with linear analysis
- VSC converters with multiple control modes (PQ, VDC_Q, VDC_VAC)
- IEEE test systems with HVDC extensions  
- Island detection and adaptive solving
- Distributed slack bus model

## Performance Optimizations
- Sparse Jacobian with pre-computed sparsity pattern
- Sparse LU factorization with symbolic reuse (UMFPACK)
- Pre-allocated solver workspace
- @inbounds and sin/cos caching

## Usage Pattern
```julia
using JuliaPowerCase, HybridACDCPowerFlow

# Build or load system using JuliaPowerCase
hps = build_ieee14_acdc()  # Returns HybridPowerSystem

# Solve
result = solve_power_flow(hps)  # Automatically creates internal solver workspace
```

Author: Tianyang Zhao
Version: 0.7.0 (Pure Application Layer Architecture)
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
export ACBus, ACBranch, DCBus, DCBranch, VSCConverter, Generator,
       # JuliaPowerCase types (primary data structures)
       HybridPowerSystem, IslandInfo,
       BusType, PQ, PV, SLACK, PQ_MODE, VDC_Q, VDC_VAC,
       # DC bus types
       DC_V, DC_P,
       # Main solver functions
       solve_power_flow, solve_dc_power_flow,
       # Utility functions
       get_bus_voltages, get_branch_flows,
       # Test systems (return HybridPowerSystem)
       build_ieee14_acdc, build_ieee24_3area_acdc, build_ieee118_acdc, build_ac_only_version,
       build_case33bw_acdc, build_case33mg_acdc, build_case69_acdc,
       build_case300_acdc, build_case2000_acdc,
       build_from_matpower, build_from_matpower_with_limits,
       # Conversion functions (for legacy HybridPowerCaseData support)
       to_hybrid_power_system, to_hybrid_system, update_results!,
       to_solver_data

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
