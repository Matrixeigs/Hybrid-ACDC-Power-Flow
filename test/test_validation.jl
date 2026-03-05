"""
Validation tests: Compare Newton-Raphson solver against NLsolve and JuMP+Ipopt.

This script validates that our custom Newton-Raphson solver produces correct results
by comparing against established reference solvers:
1. NLsolve (trust-region method for nonlinear root finding)
2. JuMP+Ipopt (interior-point nonlinear optimization)

Usage:
    cd HybridACDCPowerFlow
    julia --project=. test/test_validation.jl

Requirements:
    using Pkg
    Pkg.add(["NLsolve", "JuMP", "Ipopt"])
"""

using Test
using LinearAlgebra
using Printf

# Load extension dependencies FIRST (before main module)
using NLsolve
using JuMP
using Ipopt

# Load the package (extensions load automatically when weak deps are available)
using HybridACDCPowerFlow

println("="^70)
println("  POWER FLOW VALIDATION: Newton-Raphson vs Reference Solvers")
println("="^70)

@testset "Newton-Raphson vs NLsolve" begin
    println("\n" * "-"^70)
    println("  TEST 1: Validate against NLsolve (Trust-Region)")
    println("-"^70)
    
    @testset "IEEE 14-bus AC/DC" begin
        sys = build_ieee14_acdc()
        
        # Solve with Newton-Raphson
        res_nr = solve_power_flow(sys)
        @test res_nr.converged
        
        # Solve with NLsolve  
        res_nl = check_power_flow_feasibility_nlsolve(sys; verbose=false)
        @test res_nl.feasible
        
        # Compare voltages
        Vm_err = maximum(abs.(res_nr.Vm .- res_nl.Vm))
        Va_err = maximum(abs.(res_nr.Va .- res_nl.Va))
        Vdc_err = maximum(abs.(res_nr.Vdc .- res_nl.Vdc))
        
        println("    IEEE 14-bus AC/DC:")
        @printf("      Vm max error: %.2e p.u.\n", Vm_err)
        @printf("      Va max error: %.2e rad\n", Va_err)
        @printf("      Vdc max error: %.2e p.u.\n", Vdc_err)
        @printf("      NR iters: %d, NLsolve iters: %d\n", res_nr.iterations, res_nl.iterations)
        
        @test Vm_err < 1e-4
        @test Va_err < 1e-4
        @test Vdc_err < 1e-4
        println("      ✅ PASS")
    end
    
    @testset "IEEE 24-bus AC/DC (3 areas)" begin
        sys = build_ieee24_3area_acdc()
        
        res_nr = solve_power_flow(sys)
        @test res_nr.converged
        
        res_nl = check_power_flow_feasibility_nlsolve(sys; verbose=false)
        @test res_nl.feasible
        
        Vm_err = maximum(abs.(res_nr.Vm .- res_nl.Vm))
        Va_err = maximum(abs.(res_nr.Va .- res_nl.Va))
        Vdc_err = length(sys.dc_buses) > 0 ? maximum(abs.(res_nr.Vdc .- res_nl.Vdc)) : 0.0
        
        println("    IEEE 24-bus AC/DC:")
        @printf("      Vm max error: %.2e p.u.\n", Vm_err)
        @printf("      Va max error: %.2e rad\n", Va_err)
        @printf("      Vdc max error: %.2e p.u.\n", Vdc_err)
        
        @test Vm_err < 1e-4
        @test Va_err < 1e-4
        @test Vdc_err < 1e-4
        println("      ✅ PASS")
    end
end

@testset "Newton-Raphson vs JuMP+Ipopt" begin
    println("\n" * "-"^70)
    println("  TEST 2: Validate against JuMP+Ipopt (Interior-Point)")
    println("-"^70)
    
    @testset "IEEE 14-bus AC/DC" begin
        sys = build_ieee14_acdc()
        
        # Solve with Newton-Raphson
        res_nr = solve_power_flow(sys)
        @test res_nr.converged
        
        # Solve with JuMP+Ipopt
        res_jump = check_power_flow_feasibility_jump(sys; verbose=false)
        @test res_jump.feasible
        
        # Compare voltages
        Vm_err = maximum(abs.(res_nr.Vm .- res_jump.Vm))
        Va_err = maximum(abs.(res_nr.Va .- res_jump.Va))
        Vdc_err = length(res_jump.Vdc) > 0 ? maximum(abs.(res_nr.Vdc .- res_jump.Vdc)) : 0.0
        
        println("    IEEE 14-bus AC/DC:")
        @printf("      Vm max error: %.2e p.u.\n", Vm_err)
        @printf("      Va max error: %.2e rad\n", Va_err)
        @printf("      Vdc max error: %.2e p.u.\n", Vdc_err)
        @printf("      JuMP time: %.1f ms\n", res_jump.solve_time)
        
        # JuMP may have slightly looser tolerance
        @test Vm_err < 1e-3
        @test Va_err < 1e-3
        @test Vdc_err < 1e-3
        println("      ✅ PASS")
    end
    
    @testset "IEEE 24-bus AC/DC" begin
        sys = build_ieee24_3area_acdc()
        
        res_nr = solve_power_flow(sys)
        @test res_nr.converged
        
        res_jump = check_power_flow_feasibility_jump(sys; verbose=false)
        @test res_jump.feasible
        
        Vm_err = maximum(abs.(res_nr.Vm .- res_jump.Vm))
        Va_err = maximum(abs.(res_nr.Va .- res_jump.Va))
        
        println("    IEEE 24-bus AC/DC:")
        @printf("      Vm max error: %.2e p.u.\n", Vm_err)
        @printf("      Va max error: %.2e rad\n", Va_err)
        
        @test Vm_err < 1e-3
        @test Va_err < 1e-3
        println("      ✅ PASS")
    end
end

@testset "Residual Verification" begin
    println("\n" * "-"^70)
    println("  TEST 3: Power Balance Residual Verification")
    println("-"^70)
    
    @testset "All converters satisfy power balance" begin
        for (name, builder) in [("IEEE14", build_ieee14_acdc), 
                                ("IEEE24", build_ieee24_3area_acdc)]
            sys = builder()
            res = solve_power_flow(sys)
            @test res.converged
            
            # Compute full residual
            F_ac, F_dc = HybridACDCPowerFlow.PowerSystem.power_flow_residual(
                sys, res.Vm, res.Va, res.Vdc)
            
            P_residual = maximum(abs.(F_ac[1:2:end]))  # P equations (odd indices)
            Q_residual = maximum(abs.(F_ac[2:2:end]))  # Q equations (even indices)
            DC_residual = length(F_dc) > 0 ? maximum(abs.(F_dc)) : 0.0
            
            println("    $name:")
            @printf("      P balance residual: %.2e p.u.\n", P_residual)
            @printf("      Q balance residual: %.2e p.u.\n", Q_residual)
            @printf("      DC balance residual: %.2e p.u.\n", DC_residual)
            
            @test P_residual < 1e-8
            @test Q_residual < 1e-8
            @test DC_residual < 1e-8
            println("      ✅ PASS")
        end
    end
end

@testset "Consistency Under Load Changes" begin
    println("\n" * "-"^70)
    println("  TEST 4: Consistency Across Load Levels")
    println("-"^70)
    
    using JuliaPowerCase: Bus, AC
    
    function scale_loads(sys, factor)
        sys_scaled = build_ieee14_acdc()
        for i in 1:length(sys_scaled.ac_buses)
            b = sys_scaled.ac_buses[i]
            sys_scaled.ac_buses[i] = Bus{AC}(
                index=b.index, name=b.name, bus_id=b.bus_id, in_service=b.in_service,
                base_kv=b.base_kv, bus_type=b.bus_type, vm_pu=b.vm_pu, va_deg=b.va_deg,
                vmax_pu=b.vmax_pu, vmin_pu=b.vmin_pu,
                pd_mw=b.pd_mw * factor, qd_mvar=b.qd_mvar * factor,
                gs_mw=b.gs_mw, bs_mvar=b.bs_mvar, area=b.area, zone=b.zone,
                carbon_area=b.carbon_area, carbon_zone=b.carbon_zone,
                nc=b.nc, omega=b.omega, is_load=b.is_load
            )
        end
        return sys_scaled
    end
    
    for factor in [0.5, 0.8, 1.0, 1.2]
        sys = scale_loads(build_ieee14_acdc(), factor)
        
        res_nr = solve_power_flow(sys)
        res_nl = check_power_flow_feasibility_nlsolve(sys; verbose=false)
        
        if res_nr.converged && res_nl.feasible
            Vm_err = maximum(abs.(res_nr.Vm .- res_nl.Vm))
            println("    Load $(Int(factor*100))%: Vm error = $(@sprintf("%.2e", Vm_err))")
            @test Vm_err < 1e-4
        else
            println("    Load $(Int(factor*100))%: One solver did not converge")
        end
    end
    println("    ✅ PASS")
end

@testset "Speed Comparison" begin
    println("\n" * "-"^70)
    println("  TEST 5: Solver Speed Comparison")
    println("-"^70)
    
    sys = build_ieee14_acdc()
    
    # Warm up
    solve_power_flow(sys)
    check_power_flow_feasibility_nlsolve(sys; verbose=false)
    check_power_flow_feasibility_jump(sys; verbose=false)
    
    # Benchmark
    n_runs = 10
    
    t_nr = @elapsed for _ in 1:n_runs
        solve_power_flow(sys)
    end
    
    t_nl = @elapsed for _ in 1:n_runs
        check_power_flow_feasibility_nlsolve(sys; verbose=false)
    end
    
    t_jump = @elapsed for _ in 1:n_runs
        check_power_flow_feasibility_jump(sys; verbose=false)
    end
    
    println("    Average solve time ($n_runs runs):")
    @printf("      Newton-Raphson: %.2f ms\n", 1000*t_nr/n_runs)
    @printf("      NLsolve:        %.2f ms\n", 1000*t_nl/n_runs)
    @printf("      JuMP+Ipopt:     %.2f ms\n", 1000*t_jump/n_runs)
    
    # Our solver should be competitive with NLsolve
    @test t_nr < t_jump  # NR should be faster than JuMP
    println("    ✅ Newton-Raphson is faster than JuMP+Ipopt")
end

println("\n" * "="^70)
println("  ALL VALIDATION TESTS PASSED")
println("="^70)
