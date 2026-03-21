"""
Test all MATPOWER cases in the data/ folder.

This script tests that:
1. Each .m file can be parsed by MatpowerParser
2. The power flow converges (or at least runs without error) for each case
3. Gs/Bs shunt values are correctly included in the Ybus matrix

Usage:
    julia --project test/test_all_matpower_cases.jl
"""

using Test
using Printf

if !isdefined(Main, :HybridACDCPowerFlow)
    include(joinpath(@__DIR__, "..", "src", "HybridACDCPowerFlow.jl"))
    using .HybridACDCPowerFlow
end

const DATA_DIR = joinpath(@__DIR__, "..", "data")

# Discover all case*.m files (exclude contab_, scenarios_ files)
case_files = sort(filter(f -> startswith(basename(f), "case") && endswith(f, ".m"),
                          readdir(DATA_DIR, join=true)))

println("="^70)
println("  MATPOWER CASES POWER FLOW TEST")
println("  Testing $(length(case_files)) cases in $(DATA_DIR)")
println("="^70)

# ---------------------------------------------------------------------------
# Verify Gs/Bs shunt bug is fixed: Ybus diagonal should include bus shunts
# ---------------------------------------------------------------------------
@testset "Gs/Bs Shunt in Ybus – case14" begin
    # Build system with actual shunts (bus 9 has Bs=19 MVAR)
    sys_shunt = build_from_matpower(joinpath(DATA_DIR, "case14.m"))
    sd_shunt = HybridACDCPowerFlow.PowerSystem.to_solver_data(sys_shunt)

    # Build a reference system with explicit zero shunts to compare
    using JuliaPowerCase
    sys_noshunt = build_from_matpower(joinpath(DATA_DIR, "case14.m"))
    # Zero out bus 9 shunt by rebuilding without shunt
    for (i, bus) in enumerate(sys_noshunt.ac.buses)
        if bus.index == 9
            sys_noshunt.ac.buses[i] = JuliaPowerCase.Bus{JuliaPowerCase.AC}(
                index=bus.index, name=bus.name, bus_id=bus.bus_id,
                in_service=bus.in_service, base_kv=bus.base_kv, bus_type=bus.bus_type,
                vm_pu=bus.vm_pu, va_deg=bus.va_deg, vmax_pu=bus.vmax_pu, vmin_pu=bus.vmin_pu,
                pd_mw=bus.pd_mw, qd_mvar=bus.qd_mvar,
                gs_mw=0.0, bs_mvar=0.0,
                area=bus.area, zone=bus.zone, carbon_area=bus.carbon_area,
                carbon_zone=bus.carbon_zone, nc=bus.nc, omega=bus.omega, is_load=bus.is_load
            )
        end
    end
    sd_noshunt = HybridACDCPowerFlow.PowerSystem.to_solver_data(sys_noshunt)

    Y_with    = sd_shunt.Ybus
    Y_without = sd_noshunt.Ybus
    bs_pu = 19.0 / sd_shunt.baseMVA   # 0.19 pu

    # The imaginary part of Y[9,9] should be larger by bs_pu when shunt is included
    diff = imag(Y_with[9, 9]) - imag(Y_without[9, 9])
    @test diff ≈ bs_pu atol=1e-9
    println("  ✅ case14 Y[9,9] imag: with_shunt=$(round(imag(Y_with[9,9]),digits=5))  ",
            "without=$(round(imag(Y_without[9,9]),digits=5))  diff=$(round(diff,digits=6))  (expected=$bs_pu)")

    # Power flow with shunts should converge and give good voltages
    result = solve_power_flow(sys_shunt)
    @test result.converged
    # Bus 9 voltage reference: ≈1.056 (from MATPOWER)
    @test result.Vm[9] ≈ 1.056 atol=1e-3
    println("  ✅ case14 bus9 Vm=$(round(result.Vm[9],digits=4)) (expected≈1.056)")
end

# ---------------------------------------------------------------------------
# Run power flow on all cases
# ---------------------------------------------------------------------------
passed = 0
failed = 0
skipped = 0
failed_cases = String[]

@testset "All Matpower Cases" begin
    for filepath in case_files
        name = basename(filepath)
        # Skip very large cases to avoid excessive runtime (>25k buses)
        nbus_approx = 0
        try
            m = match(r"case.*?(\d+)", name)
            m !== nothing && (nbus_approx = parse(Int, m.captures[1]))
        catch; end
        if nbus_approx > 25000
            @testset "$name" begin
                @test_skip "Skipping large case ($nbus_approx buses)"
            end
            skipped += 1
            continue
        end

        @testset "$name" begin
            result_ok = false
            err_msg = ""
            local sys, result

            # Parse
            parse_ok = try
                sys = build_from_matpower(filepath)
                true
            catch e
                err_msg = "PARSE ERROR: $(sprint(showerror, e))"
                false
            end

            @test parse_ok
            if !parse_ok
                @printf("  ❌ %-40s  %s\n", name, err_msg)
                failed += 1
                push!(failed_cases, name * " [PARSE]")
                return
            end

            # Solve
            solve_ok = try
                result = solve_power_flow(sys; max_iter=200, tol=1e-6)
                true
            catch e
                err_msg = "SOLVE ERROR: $(sprint(showerror, e))"
                false
            end

            @test solve_ok
            if !solve_ok
                @printf("  ❌ %-40s  %s\n", name, err_msg)
                failed += 1
                push!(failed_cases, name * " [SOLVE]")
                return
            end

            # Check convergence
            @test result.converged

            if result.converged
                @printf("  ✅ %-40s  iters=%3d  Vm=[%.4f, %.4f]\n",
                        name, result.iterations,
                        minimum(result.Vm), maximum(result.Vm))
                passed += 1
            else
                @printf("  ⚠️  %-40s  DID NOT CONVERGE (iters=%d, res=%.2e)\n",
                        name, result.iterations, result.residual)
                failed += 1
                push!(failed_cases, name * " [NO CONV]")
            end
        end
    end
end

println()
println("="^70)
println("  SUMMARY")
println("="^70)
@printf("  Passed:  %d\n", passed)
@printf("  Failed:  %d\n", failed)
@printf("  Skipped: %d (>25k buses)\n", skipped)
@printf("  Total:   %d\n", passed + failed + skipped)

if !isempty(failed_cases)
    println("\n  Failed cases:")
    for c in failed_cases
        println("    - $c")
    end
end
println()
