# Testing and Validation

## Main Test Entry

```bash
julia --project=HybridACDCPowerFlow HybridACDCPowerFlow/test/runtests.jl
```

Current scope includes:

- core AC/DC solve checks
- theorem validation workloads
- control-mode and loss-model sweeps
- islanding and adaptive-flow checks

## JuliaPowerCase Integration Suite

```bash
julia --project=HybridACDCPowerFlow HybridACDCPowerFlow/test/test_jpc_integration.jl
```

Covers:

- `HybridPowerSystem` overloads
- parity of converted `HybridSystem` solve path
- distributed slack and adaptive integration behavior

## Optional Extra Suites

`runtests.jl` conditionally includes extra tests when:

```bash
HYBRID_ACDC_RUN_EXTRA_TESTS=true
```

Extra files:

- `test/test_enhanced.jl`
- `test/test_distributed_slack.jl`

## Recommended CI Matrix

- Core tests always (`runtests.jl`)
- Integration tests always (`test_jpc_integration.jl`)
- Extension tests in dedicated job with JuMP/Ipopt/NLsolve
