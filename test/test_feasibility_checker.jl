"""
Test script for JuMP-based feasibility checker
"""

using HybridACDCPowerFlow

using Printf

println("="^80)
println("TESTING JUMP FEASIBILITY CHECKER")
println("="^80)

# Test 1: Base case (should be feasible)
println("\n" * "="^80)
println("Test 1: IEEE 24-bus base case")
println("="^80)

sys1 = build_ieee24_3area_acdc()
result1 = check_power_flow_feasibility(sys1; verbose=true, max_time=10.0)

println("\n📊 Result:")
println("  Feasible: ", result1.feasible)
println("  Status: ", result1.status)
@printf("  Solve time: %.3f s\n", result1.solve_time)
@printf("  Load margin: %.4f p.u.\n", result1.load_margin)

# Test 2: Heavy load (may be infeasible)
println("\n" * "="^80)
println("Test 2: Heavy load scenario (+50%)")
println("="^80)

sys2 = build_ieee24_3area_acdc()
for i in 1:length(sys2.ac_buses)
    bus = sys2.ac_buses[i]
    sys2.ac_buses[i] = ACBus(
        bus.id, bus.type,
        bus.Pd * 1.5, bus.Qd * 1.5,  # 50% load increase
        bus.Pg, bus.Qg, bus.Vm, bus.Va, bus.area
    )
end

result2 = check_power_flow_feasibility(sys2; verbose=true, max_time=10.0)

println("\n📊 Result:")
println("  Feasible: ", result2.feasible)
println("  Status: ", result2.status)
@printf("  Solve time: %.3f s\n", result2.solve_time)
@printf("  Load margin: %.4f p.u.\n", result2.load_margin)

# Test 3: Line outage
println("\n" * "="^80)
println("Test 3: N-2 contingency (2 line outages)")
println("="^80)

sys3 = build_ieee24_3area_acdc()
# Disable first two branches
sys3.ac_branches[1] = ACBranch(
    sys3.ac_branches[1].from, sys3.ac_branches[1].to,
    sys3.ac_branches[1].r, sys3.ac_branches[1].x, sys3.ac_branches[1].b,
    sys3.ac_branches[1].tap, false
)
sys3.ac_branches[2] = ACBranch(
    sys3.ac_branches[2].from, sys3.ac_branches[2].to,
    sys3.ac_branches[2].r, sys3.ac_branches[2].x, sys3.ac_branches[2].b,
    sys3.ac_branches[2].tap, false
)
sys3.Ybus = build_admittance_matrix(sys3)

result3 = check_power_flow_feasibility(sys3; verbose=true, max_time=10.0)

println("\n📊 Result:")
println("  Feasible: ", result3.feasible)
println("  Status: ", result3.status)
@printf("  Solve time: %.3f s\n", result3.solve_time)
@printf("  Load margin: %.4f p.u.\n", result3.load_margin)

println("\n" * "="^80)
println("FEASIBILITY CHECKER TESTS COMPLETE")
println("="^80)
