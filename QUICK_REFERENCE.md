# Quick Reference - HybridACDCPowerFlow.jl v0.2.0

⚡ **One-page cheat sheet for enhanced features**

---

## 🚀 Installation

```julia
using Pkg
Pkg.add(path="/path/to/HybridACDCPowerFlow")

using HybridACDCPowerFlow
```

---

## 📌 Basic Usage

### Standard Power Flow (v0.1.0 Compatible)
```julia
sys = build_ieee14_acdc()
result = solve_power_flow(sys)
println("Converged: $(result.converged)")
```

### Enhanced Adaptive Power Flow (v0.2.0)
```julia
sys = build_ieee24_3area_acdc()

Q_limits = create_default_Q_limits(sys; Qmin_default=-0.5, Qmax_default=1.0)

options = PowerFlowOptions(
    enable_pv_pq_conversion=true,
    enable_auto_swing_selection=true,
    enable_converter_mode_switching=true
)

result = solve_power_flow_adaptive(sys; options=options, Q_limits=Q_limits)
```

---

## 🏝️ Feature 1: Island Detection

```julia
# Detect islands
islands = detect_islands(sys)

# Print summary
print_island_summary(islands, sys)

# Access island data
for island in islands
    println("Island $(island.id):")
    println("  AC buses: $(island.ac_buses)")
    println("  DC buses: $(island.dc_buses)")
    println("  Slack: $(island.ac_slack_bus)")
end
```

**Output Example**:
```
Island 1: 14 AC buses, 2 DC buses
Island 2: 10 AC buses, 2 DC buses
```

---

## ⚡ Feature 2: PV → PQ Conversion

```julia
# Create reactive limits
Q_limits = Dict(
    2 => ReactiveLimit(-0.4, 0.8),  # Generator 2
    3 => ReactiveLimit(-0.3, 0.6)   # Generator 3
)

# Or use defaults
Q_limits = create_default_Q_limits(sys; Qmin_default=-0.5, Qmax_default=1.0)

# Solve with automatic conversion
options = PowerFlowOptions(enable_pv_pq_conversion=true)
result = solve_power_flow_adaptive(sys; options=options, Q_limits=Q_limits)

# Check violations
violations, Q_actual = check_reactive_limits(
    sys, result.Vm, result.Va, result.Vdc, Q_limits
)
println("Violations: $(length(violations))")
```

---

## 🎯 Feature 3: Auto Swing Bus Selection

```julia
# Convert slack to PV (simulate slack loss)
sys.ac_buses[1] = ACBus(1, PV, ...)

# Detect islands
islands = detect_islands(sys)

# Auto-select slack
if !islands[1].has_ac_slack
    slack_bus = auto_select_swing_bus(sys, islands[1])
    println("Auto-selected slack: bus $slack_bus")
end
```

**Selection Criteria** (priority order):
1. Largest Pg_max
2. Most connections
3. Lowest bus number

---

## 🔄 Feature 4: Converter Mode Switching

```julia
# Solve power flow
result = solve_power_flow(sys)

# Auto-switch converter modes
switched = auto_switch_converter_mode!(
    sys, result.Vm, result.Va, result.Vdc;
    V_low_threshold=0.95,   # Switch to VDC_VAC if V < 0.95
    V_high_threshold=1.05,  # Switch back if V > 1.05
    S_limit_frac=0.90,      # Switch to PQ if S > 90% Smax
    hysteresis=0.02         # Prevent oscillation
)

println("Switched $(length(switched)) converters")
for idx in switched
    println("  Converter $idx → $(sys.converters[idx].mode)")
end
```

**Mode Transitions**:
```
PQ_MODE ↔ VDC_VAC ↔ VDC_Q
```

---

## 🌐 Feature 5: Multi-Island Power Flow

```julia
# Solve each island independently
island_results = solve_power_flow_islanded(sys)

for (i, res) in enumerate(island_results)
    println("Island $i:")
    println("  Converged: $(res.converged)")
    println("  Voltage range: [$(minimum(res.Vm)), $(maximum(res.Vm))]")
end
```

---

## 📊 Data Structures

### PowerFlowOptions
```julia
PowerFlowOptions(
    max_iter=20,                           # Max Newton iterations
    tol=1e-8,                              # Convergence tolerance
    enable_pv_pq_conversion=false,         # Enable PV→PQ
    enable_auto_swing_selection=false,     # Enable auto slack
    enable_converter_mode_switching=false, # Enable mode switching
    verbose=false                          # Print details
)
```

### IslandInfo
```julia
island.id                 # Island number
island.ac_buses           # Vector{Int}
island.dc_buses           # Vector{Int}
island.converters         # Vector{Int}
island.has_ac_slack       # Bool
island.ac_slack_bus       # Int (or 0)
```

### ReactiveLimit
```julia
ReactiveLimit(Qmin, Qmax)  # Min/max reactive power (p.u.)
```

---

## 🧪 Test Systems

```julia
# IEEE 14-bus + 2 DC buses + 2 converters
sys = build_ieee14_acdc()

# IEEE 24-bus 3-area + 4 DC buses + 4 converters
sys = build_ieee24_3area_acdc()

# IEEE 118-bus + 6 DC buses + 5 converters
sys = build_ieee118_acdc()

# AC-only version (remove converters)
sys_ac = build_ac_only_version(sys)
```

---

## 🔍 Result Analysis

```julia
result = solve_power_flow_adaptive(sys; options, Q_limits)

# Convergence
result.converged    # Bool

# Voltages
result.Vm           # AC voltage magnitudes
result.Va           # AC voltage angles (rad)
result.Vdc          # DC voltages

# Islands (if adaptive)
result.islands      # Vector{IslandInfo}

# Iterations
result.iterations   # Total Newton iterations

# Extract voltages
Vm, Va, Vdc = get_bus_voltages(result)

# Branch flows
Pij, Qij = get_branch_flows(sys, result)
```

---

## ⚙️ Common Workflows

### 1. N-k Contingency Analysis
```julia
sys = build_ieee24_3area_acdc()

# Remove k branches
for i in 1:k
    sys.ac_branches[i] = ACBranch(..., false)  # Set status=false
end

# Solve adaptive
result = solve_power_flow_adaptive(sys; options, Q_limits)

# Check islanding
islands = detect_islands(sys)
println("System split into $(length(islands)) islands")
```

### 2. Generator Reactive Capability Study
```julia
# Tight limits
Q_limits = create_default_Q_limits(sys; Qmin_default=-0.3, Qmax_default=0.5)

options = PowerFlowOptions(enable_pv_pq_conversion=true)
result = solve_power_flow_adaptive(sys; options, Q_limits)

# Find violations
violations, Q_actual = check_reactive_limits(sys, result.Vm, result.Va, result.Vdc, Q_limits)

for (bus, Q) in violations
    lim = Q_limits[bus]
    println("Bus $bus: Q=$(Q) violates [$(lim.Qmin), $(lim.Qmax)]")
end
```

### 3. Extreme Event Simulation
```julia
# Progressive fault
for i in 1:10
    sys.ac_branches[i] = ACBranch(..., false)
    
    result = solve_power_flow_adaptive(sys; options, Q_limits)
    
    if !result.converged
        println("System unstable after $i faults")
        break
    end
    
    islands = detect_islands(sys)
    println("Fault $i: $(length(islands)) islands")
end
```

---

## 📖 Documentation

- **Full features**: `ENHANCED_FEATURES.md`
- **Changelog**: `CHANGELOG_v0.2.0.md`
- **Theory**: `docs/TheoreticalFoundations.pdf`
- **Tests**: `test/test_enhanced.jl`

---

## 🐛 Troubleshooting

### "UndefVarError: function not defined"
→ Restart Julia and reload module after updates

### "BoundsError in detect_islands"
→ Fixed in v0.2.0; update to latest version

### "Does not converge with tight Q limits"
→ Expected; relax limits or increase `max_iter`

### "Mode switching oscillates"
→ Increase `hysteresis` parameter

---

## 📊 Performance Tips

1. **Start simple**: Test with `build_ieee14_acdc()` first
2. **Enable selectively**: Only enable features you need
3. **Relax tolerance**: Use `tol=1e-6` for fault scenarios
4. **Check islands**: Run `detect_islands()` before adaptive solver

---

## 🎯 Validation

Run comprehensive test suite:
```bash
julia test/test_enhanced.jl
```

Expected:
```
Test Summary: | Pass  Total  Time
Enhanced      |   19     19  1.5s

ALL TESTS PASSED ✅
```

---

## 📜 License

MIT License - Free for academic and commercial use

---

**Version**: 0.2.0  
**Status**: Production Ready  
**Tests**: 19/19 Passing
