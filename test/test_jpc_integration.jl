"""
Test suite for JuliaPowerCase (HybridPowerSystem) integration.

Tests the HybridPowerSystem overloads added in v0.5.0:
1. solve_power_flow(::HybridPowerSystem)
2. solve_power_flow_adaptive(::HybridPowerSystem)
3. solve_power_flow_islanded(::HybridPowerSystem)
4. solve_power_flow_distributed_slack(::HybridPowerSystem)
"""

using Test
using LinearAlgebra

# Load modules
if !isdefined(Main, :HybridACDCPowerFlow)
    include("../src/HybridACDCPowerFlow.jl")
    using .HybridACDCPowerFlow
end

using JuliaPowerCase
using JuliaPowerCase: AC, DC, Bus, Branch, Generator, VSCConverter, PowerSystem,
                      HybridPowerSystem, BusType, PQ_BUS, PV_BUS, REF_BUS

# ═══════════════════════════════════════════════════════════════════════════════
#  Test Fixtures - Build HybridPowerSystem from JuliaPowerCase types
# ═══════════════════════════════════════════════════════════════════════════════

"""
Build a simple 5-bus AC / 2-bus DC hybrid system using JuliaPowerCase types.
"""
function build_simple_jpc_system()
    # ── AC Subsystem ──────────────────────────────────────────────────────
    ac = PowerSystem{AC}()
    ac.base_mva = 100.0
    
    # AC Buses: 1=slack, 2=PV, 3-5=PQ
    push!(ac.buses, Bus{AC}(index=1, bus_type=REF_BUS, vm_pu=1.06, va_deg=0.0, 
                            pd_mw=0.0, qd_mvar=0.0, base_kv=230.0))
    push!(ac.buses, Bus{AC}(index=2, bus_type=PV_BUS, vm_pu=1.04, va_deg=0.0,
                            pd_mw=20.0, qd_mvar=10.0, base_kv=230.0))
    push!(ac.buses, Bus{AC}(index=3, bus_type=PQ_BUS, vm_pu=1.0, va_deg=0.0,
                            pd_mw=45.0, qd_mvar=15.0, base_kv=230.0))
    push!(ac.buses, Bus{AC}(index=4, bus_type=PQ_BUS, vm_pu=1.0, va_deg=0.0,
                            pd_mw=40.0, qd_mvar=5.0, base_kv=230.0))
    push!(ac.buses, Bus{AC}(index=5, bus_type=PQ_BUS, vm_pu=1.0, va_deg=0.0,
                            pd_mw=60.0, qd_mvar=10.0, base_kv=230.0))
    
    # AC Branches
    push!(ac.branches, Branch{AC}(index=1, from_bus=1, to_bus=2, r_pu=0.02, x_pu=0.06, b_pu=0.06, in_service=true))
    push!(ac.branches, Branch{AC}(index=2, from_bus=1, to_bus=3, r_pu=0.08, x_pu=0.24, b_pu=0.05, in_service=true))
    push!(ac.branches, Branch{AC}(index=3, from_bus=2, to_bus=3, r_pu=0.06, x_pu=0.18, b_pu=0.04, in_service=true))
    push!(ac.branches, Branch{AC}(index=4, from_bus=3, to_bus=4, r_pu=0.06, x_pu=0.18, b_pu=0.04, in_service=true))
    push!(ac.branches, Branch{AC}(index=5, from_bus=4, to_bus=5, r_pu=0.04, x_pu=0.12, b_pu=0.03, in_service=true))
    
    # Generators
    push!(ac.generators, Generator(index=1, bus=1, pg_mw=130.0, qg_mvar=0.0, vg_pu=1.06,
                                   pmax_mw=200.0, pmin_mw=0.0, qmax_mvar=100.0, qmin_mvar=-100.0, in_service=true))
    push!(ac.generators, Generator(index=2, bus=2, pg_mw=50.0, qg_mvar=0.0, vg_pu=1.04,
                                   pmax_mw=100.0, pmin_mw=0.0, qmax_mvar=50.0, qmin_mvar=-50.0, in_service=true))
    
    # ── DC Subsystem ──────────────────────────────────────────────────────
    dc = PowerSystem{DC}()
    dc.base_mva = 100.0
    
    # DC Buses (1-indexed for solver compatibility)
    push!(dc.buses, Bus{DC}(index=1, bus_type=PV_BUS, vm_pu=1.0, pd_mw=0.0))
    push!(dc.buses, Bus{DC}(index=2, bus_type=PQ_BUS, vm_pu=1.0, pd_mw=20.0))
    
    # DC Branch
    push!(dc.branches, Branch{DC}(index=1, from_bus=1, to_bus=2, r_pu=0.01, in_service=true))
    
    # ── VSC Converters ────────────────────────────────────────────────────
    vscs = VSCConverter[]
    push!(vscs, VSCConverter(
        index=1, bus_ac=3, bus_dc=1, p_mw=30.0, q_mvar=10.0,
        vm_ac_pu=1.0, vm_dc_pu=1.0, pmax_mw=100.0, pmin_mw=-100.0,
        qmax_mvar=50.0, qmin_mvar=-50.0, control_mode=:pq, in_service=true
    ))
    push!(vscs, VSCConverter(
        index=2, bus_ac=5, bus_dc=2, p_mw=15.0, q_mvar=5.0,
        vm_ac_pu=1.0, vm_dc_pu=1.0, pmax_mw=100.0, pmin_mw=-100.0,
        qmax_mvar=50.0, qmin_mvar=-50.0, control_mode=:pq, in_service=true
    ))
    
    # ── Assemble HybridPowerSystem ────────────────────────────────────────
    return HybridPowerSystem(
        ac=ac, dc=dc, vsc_converters=vscs,
        name="Simple 5AC/2DC Test System", base_mva=100.0
    )
end

"""
Build a system that can be islanded (for islanded solver tests).
"""
function build_islandable_jpc_system()
    hps = build_simple_jpc_system()
    
    # Add extra connectivity so we can island by removing one branch
    push!(hps.ac.branches, Branch{AC}(index=6, from_bus=2, to_bus=4, r_pu=0.05, x_pu=0.15, b_pu=0.03, in_service=true))
    
    return hps
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Tests
# ═══════════════════════════════════════════════════════════════════════════════

@testset "JuliaPowerCase HybridPowerSystem Integration" begin
    
    @testset "1. Basic solve_power_flow(HybridPowerSystem)" begin
        println("\n" * "="^70)
        println("【TEST 1】BASIC POWER FLOW ON HybridPowerSystem")
        println("="^70)
        
        hps = build_simple_jpc_system()
        
        # Solve using the HybridPowerSystem overload
        result = solve_power_flow(hps)
        
        println("  Converged: $(result.converged)")
        println("  Iterations: $(result.iterations)")
        println("  Final residual: $(result.residual)")
        
        @test result.converged == true
        @test result.iterations < 20
        @test result.residual < 1e-6
        
        # Check voltage vectors have correct dimension
        @test length(result.Vm) == 5  # 5 AC buses
        @test length(result.Vdc) == 2  # 2 DC buses
        
        println("  ✅ PASS\n")
    end
    
    @testset "2. solve_power_flow_adaptive(HybridPowerSystem)" begin
        println("="^70)
        println("【TEST 2】ADAPTIVE POWER FLOW ON HybridPowerSystem")
        println("="^70)
        
        hps = build_simple_jpc_system()
        
        # Create Q limits
        Q_limits = Dict{Int, ReactiveLimit}(
            2 => ReactiveLimit(-0.5, 0.5)  # Tight limits for PV bus 2
        )
        
        # Solve with adaptive features
        options = PowerFlowOptions(
            enable_pv_pq_conversion=true,
            verbose=false
        )
        result = solve_power_flow_adaptive(hps; options=options, Q_limits=Q_limits)
        
        println("  Converged: $(result.converged)")
        println("  Iterations: $(result.iterations)")
        
        @test result.converged == true
        @test result.iterations < 50
        
        println("  ✅ PASS\n")
    end
    
    @testset "3. solve_power_flow_islanded(HybridPowerSystem)" begin
        println("="^70)
        println("【TEST 3】ISLANDED POWER FLOW ON HybridPowerSystem")
        println("="^70)
        
        hps = build_islandable_jpc_system()
        
        # Disable one branch to potentially create an island
        hps.ac.branches[2] = Branch{AC}(
            index=2, from_bus=hps.ac.branches[2].from_bus, to_bus=hps.ac.branches[2].to_bus,
            r_pu=hps.ac.branches[2].r_pu, x_pu=hps.ac.branches[2].x_pu, b_pu=hps.ac.branches[2].b_pu,
            in_service=false  # Disable this branch
        )
        
        # Solve with islanded solver (returns array of island results)
        island_results = solve_power_flow_islanded(hps)
        
        println("  Number of islands: $(length(island_results))")
        all_converged = all(r.converged for r in island_results)
        println("  All converged: $all_converged")
        
        # Should converge (system still connected via other paths)
        @test all_converged == true
        @test length(island_results) >= 1
        
        println("  ✅ PASS\n")
    end
    
    @testset "4. solve_power_flow_distributed_slack(HybridPowerSystem)" begin
        println("="^70)
        println("【TEST 4】DISTRIBUTED SLACK ON HybridPowerSystem")  
        println("="^70)
        
        hps = build_simple_jpc_system()
        
        # Create DistributedSlack object with generator participation factors
        # Bus 1 and 2 have generators; we split slack between them
        dslack = DistributedSlack([1, 2], [0.6, 0.4])
        
        result = solve_power_flow_distributed_slack(hps, dslack)
        
        println("  Converged: $(result.converged)")
        println("  Iterations: $(result.iterations)")
        
        @test result.converged == true
        
        println("  ✅ PASS\n")
    end
    
    @testset "5. Result Parity: HybridPowerSystem vs HybridSystem" begin
        println("="^70)
        println("【TEST 5】RESULT PARITY CHECK")
        println("="^70)
        
        hps = build_simple_jpc_system()
        
        # Solve via HybridPowerSystem overload
        result_jpc = solve_power_flow(hps)
        
        # Solve via internal conversion (to_solver_system)
        sys_internal = to_solver_system(hps)
        result_internal = solve_power_flow(sys_internal)
        
        println("  JPC converged: $(result_jpc.converged), Internal converged: $(result_internal.converged)")
        println("  JPC iterations: $(result_jpc.iterations), Internal iterations: $(result_internal.iterations)")
        
        @test result_jpc.converged == result_internal.converged
        @test result_jpc.iterations == result_internal.iterations
        
        # Voltage results should be identical
        @test result_jpc.Vm ≈ result_internal.Vm atol=1e-10
        @test result_jpc.Va ≈ result_internal.Va atol=1e-10
        @test result_jpc.Vdc ≈ result_internal.Vdc atol=1e-10
        
        println("  Max Vm difference: $(maximum(abs.(result_jpc.Vm - result_internal.Vm)))")
        println("  Max Va difference: $(maximum(abs.(result_jpc.Va - result_internal.Va)))")
        
        println("  ✅ PASS\n")
    end
    
    @testset "6. to_solver_system Conversion" begin
        println("="^70)
        println("【TEST 6】to_solver_system CONVERSION")
        println("="^70)
        
        hps = build_simple_jpc_system()
        
        # Convert to solver system
        sys = to_solver_system(hps)
        
        # Check dimensions match
        @test length(sys.ac_buses) == length(hps.ac.buses)
        @test length(sys.dc_buses) == length(hps.dc.buses)
        @test length(sys.ac_branches) == length(hps.ac.branches)
        @test length(sys.dc_branches) == length(hps.dc.branches)
        @test length(sys.converters) == length(hps.vsc_converters)
        
        println("  AC buses: $(length(sys.ac_buses))")
        println("  DC buses: $(length(sys.dc_buses))")
        println("  AC branches: $(length(sys.ac_branches))")
        println("  DC branches: $(length(sys.dc_branches))")
        println("  Converters: $(length(sys.converters))")
        
        println("  ✅ PASS\n")
    end
    
end

println("\n" * "="^70)
println("ALL JuliaPowerCase INTEGRATION TESTS COMPLETED")
println("="^70)
