# 04 Flux check: which ABM flows deviate from the ODE rates, and by how much?

## Question

In experiment 03 (`beta` tuned to 500/100k, close to threshold), the ABM drifted below the ODE equilibrium. Prevalence was about 5% low at daily `dt`, 16% at weekly and 33% at monthly (mean prevalence, ABM ÷ ODE = 0.95, 0.84 and 0.67). That corresponds to only 0.5%, 1.8% and 3.8% less effective transmission, but prevalence there responds about 10× more strongly than `beta`. Which flows cause this: transmission (the `beta` conversion) or the natural-history transitions (step ordering and cascading)? How does each scale with `dt`?

## Design

- Same configuration as 03: tbsim default natural history, ODE `beta` = 9.917, ABM started at the ODE equilibrium, 50,000 agents.
- 10 years at `dt` = 1 day, 1 week and 1 month, with 4 seeds each. 10 years is short enough that the state mix stays near equilibrium.
- An analyzer records each living agent's state at the start and end of every step, and counts the transitions between them (`tbode.FluxCounts`).
- **Effective rates.** For each natural-history flow A→B: effective rate = (number of A→B transitions) / (stock of A at step start × dt). This is compared with the ODE rate. Cascades show up as direct "skip" flows such as INFECTION→ASYMPTOMATIC via NON_INFECTIOUS.
- **Transmission.** The realised per-step infection probability of SUSCEPTIBLE and CLEARED agents is compared with the ODE prediction, `1 − exp(−foi·rr·dt)`, using the ABM's own infectious stocks at the start of each step.

## Versions

tbsim main @ `b92d7e8` (v0.10.3), starsim @ `3521adb4` (v3.6.1).
