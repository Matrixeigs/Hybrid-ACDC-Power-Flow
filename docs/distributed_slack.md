# Distributed Slack

## Data Structure

`DistributedSlack` defines:

- `participating_buses`
- `participation_factors` (must sum to `1.0`)
- `reference_bus`
- `max_participation_P`

## Factor Generation

```julia
dslack = create_participation_factors(sys; method=:capacity)
```

Methods:

- `:capacity`
- `:droop`
- `:equal`

## Solver Variants

### `solve_power_flow_distributed_slack`

- practical simplified flow
- uses standard solve path for robustness
- computes slack sharing summary post solve

### `solve_power_flow_distributed_slack_full`

- augmented Newton system with one slack variable `lambda_slack`
- includes participation terms directly in equations
- optional generator participation limits (`enforce_limits=true`)

## Returned Slack Metrics

Both variants return at least:

- `distributed_slack_P::Dict{Int,Float64}`
- `total_slack_P::Float64`

The full variant also returns `hit_limits`.

## Practical Recommendation

- use simplified variant for routine studies
- use full variant when participation-limit dynamics are required
