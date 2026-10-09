# 05 Exact per-step transitions: does fixing the step ordering remove the gap?

## Question

Experiment 04 traced the ABM's shortfall to tbsim's step ordering. Transitions are applied immediately in a fixed order, so agents can continue forward through several states within a step but not backward, and forward entrants to ASYMPTOMATIC and SYMPTOMATIC lose about half a step there. `tbode.TBExact` is a version of `tbsim.TB` that instead draws each agent's end-of-step state from the exact one-step transition probabilities of the continuous-time model, `expm(Q·dt)`, where Q holds the ODE's natural-history rates. With it, does the ABM match the ODE equilibrium at daily, weekly and monthly `dt`?

At daily `dt`, exp 03 found `tbsim.TB` at 0.95 ± 0.02 of the ODE prevalence (40 reps). Daily `TBExact` shows whether that remaining gap is step ordering (`TBExact` ≈ 1.00) or something both versions share (`TBExact` ≈ 0.95). `TBExact` still splits transmission and transitions into steps: the force of infection is fixed at its start-of-step value, and agents infected in a step also go through that step's transitions. Comparing the three `dt`s shows whether this leaves any `dt` dependence.

## Design

- **Setup:** same as 03 (tbsim default natural history, ODE `beta` = 9.917 for 500/100k, ABM started at the ODE equilibrium, 100 years, 50,000 agents). Uses `tb_class=tbode.TBExact`; transmission and everything else are unchanged.
- **Arms:** `dt` = 1 day, 1 week and 1 month, 40 reps each (seeds 0–39).
- **Comparison:** each rep's mean over the last 50 years. Report the mean ± SE across reps divided by the ODE equilibrium, a z-score, and the median and 10–90% range. Compare with `tbsim.TB` from exp 03: 40 daily reps, 10 weekly and 10 monthly.
- **Saved data:** `outputs.obj` holds the per-rep time series and the ODE run.
- **History:** the first version of this experiment ran 10 reps at weekly and monthly `dt` only. It was rerun with 40 reps and a daily arm (2026-10-09) after exp 03 showed the rep-to-rep spread needs about 40 reps for a ±2% SE. Seeds 0–9 are the same as in the first version.

**Flux check (`run_flux.py`):** experiment 04's per-step comparison against `expm(Q·dt)`, using the shared `tbode.run_flux` / `tbode.flux_tables`. Run for `TBExact` and `tbsim.TB` with the same seeds, at `dt` = 1 week and 1 month, over 10 years with 50,000 agents and 20 seeds each. Each flow gets a z-score (ratio − 1) / SE; rare flows also get exact Poisson tests.

**Unit tests:** `tests/test_exact.py`.

## Versions

tbsim main @ `b92d7e8` (v0.10.3), starsim @ `3521adb4` (v3.6.1).
