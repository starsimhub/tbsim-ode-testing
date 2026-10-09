# 01 Baseline: simplest matched ABM vs ODE

## Question

With the ABM set up as close to the ODE as possible, how well does single-strain `tbsim.TB` reproduce the ODE? Any remaining gap is the floor that later features add to.

## Design

- **ODE:** `tbode.TBODE` (a generalised `TB_ODE`) with LSHTM defaults. It is integrated with odeint and reported monthly.
  - Main comparison: `replace_tb_deaths=False`, so births = μN, matching the ABM's demography.
  - `replace_tb_deaths=True` (the original `TB_ODE`) is also shown, to measure how much that one choice matters.
- **ABM:** `tbode.make_abm`, with:
  - `tbsim.TB` using the same natural-history rates as the ODE;
  - `ss.RandomNet(n_contacts=10, dur=0)`, with `beta` converted by `beta_edge_from_ode`;
  - `ss.Births` and `ss.Deaths` at a constant μ = 1/70 per year;
  - a monthly `dt`;
  - 1% of agents starting in `SYMPTOMATIC`;
  - reinfection of latent and non-infectious agents off;
  - no treatment.
- **Scenarios:** ODE `beta` of 9 (the default, just above threshold), 12 and 20. 10 seeds × 20,000 agents each, plus 5 seeds × 100,000 agents at `beta` = 9.
- **Period:** 1900–2100.
- **Outputs:** proportion in each state (the ODE's CLEARED + RECOVERED + TREATED are compared with the ABM's CLEARED), active prevalence, and TB deaths per person-year. ABM results are shown as the median and the 10–90% range across seeds.

## Versions

tbsim main @ `b92d7e8` (v0.10.3), starsim @ `3521adb4` (v3.6.1).
