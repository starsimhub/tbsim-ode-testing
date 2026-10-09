# 02 Time step: is the baseline gap a discretisation effect?

## Question

In experiment 01, the ABM ran about 10–35% below the ODE for active prevalence and latent infection, at every `beta`, and CLEARED was about 20% high in the first year. Does the gap close as the ABM `dt` shrinks?

## Hypothesis

The gap comes from how a `tbsim.TB` step is ordered. Transmission comes first, then the transitions, which are applied immediately in a fixed order (INFECTION → NON_INFECTIOUS → ASYMPTOMATIC → SYMPTOMATIC). So newly infected agents can leave INFECTION in the step they enter it, and agents can pass through several states in one step, but only in the forward direction. With mean stays of about 4–6 months in INFECTION and ASYMPTOMATIC, a monthly step is coarse. If this is right, the gap should shrink roughly in proportion to `dt`.

## Design

Same setup as 01 (`tbode.make_abm`, ODE with `replace_tb_deaths=False`), ODE `beta` = 12 and 20, 1900–1960 (both are at equilibrium by about 1950).

- `dt` = 1 month and 1 week: 10 seeds × 20,000 agents
- `dt` = 1 day: 5 seeds × 20,000 agents

## Versions

tbsim main @ `b92d7e8` (v0.10.3), starsim @ `3521adb4` (v3.6.1).
