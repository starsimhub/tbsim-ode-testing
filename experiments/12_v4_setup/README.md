# 12 Setting up the v4 track: does tbsim v4's ODE reproduce ours, and how should agents be matched to it?

**Track: v4.** tbsim `origin/v4` @ `dcfb522` (unreleased; it still reports `__version__ = '0.10.3'`), starsim `origin/v4` @ `80c81f53` (4.0.0rc1 + 21 commits). The v0.10 reference values come from the v0.10 environment: tbsim v0.10.3 @ `b92d7e8`, starsim v3.6.1 @ `3521adb4`. Each script writes the exact versions to a `versions*.json` file next to its outputs.

## Background

In tbsim v4, `tbsim.TB` is declared as states and transitions, so the same model runs as an ODE (`solver='ode'`) or as agents (`solver='abm'`). Experiments 01–11 compared v0.10's agent-based `tbsim.TB` with our own ODE, `tbode.TBODE` (a generalisation of tbsim's reference `TB_ODE`), and with `tbode.TBExact`, our corrected ABM that samples each step exactly. In v4 the comparison is instead between the two solvers of one model.

Before anything in v4 can be compared with the earlier results, two things need checking:

- **The ODE is the same model.** Otherwise "v4 agents ÷ v4 ODE" can't be set beside "v0.10 agents ÷ `TBODE`".
- **The settings that matched agents to the ODE in v0.10 still hold in v4.** These are the mean degree of `ss.RandomNet`, the per-contact `beta`, and starting agents at the ODE equilibrium.

The baseline is exps 03–05's: tbsim's default natural history, no treatment, constant birth and death rates (`mu` = 1/70 per year, with TB deaths not replaced), and `beta` calibrated to 500/100k active (ASY + SYM) prevalence at equilibrium. That gives `beta` = 9.9173 in `TBODE`.

## Questions and design

### A. Does the v4 ODE reproduce our ODE?

- **A1. The reference.** v4 `TB(solver='ode')` with tbsim's defaults (births replacing every death, 1,000 of 100,000 symptomatic at the start), 1500–2020, against `tbsim.compartmental.TB_ODE` (which v4 keeps as its reference). Criterion: every state within 1e-4 relative.
- **A2. Our baseline, over time.** v4 ODE against `TBODE(replace_tb_deaths=False)` from a 1% symptomatic seed, 1900–2100, monthly output, for:
  - the baseline (tbsim defaults, `beta` = 9.9173);
  - set B's natural history (tb_LAI_TPT's best fit, with its reinfection risk ratios) at `k` = 0, with set B's `beta` (25.0). This exercises different rates and the reinfection terms. Set B proper has `k` = 5.846, which v4 has no exact counterpart for; that is a later experiment.

  `TBODE` is run in the v0.10 environment and its trajectories saved (`export_reference.py`). Criterion: state shares within 1e-4 relative from the second year on, and within 1e-5 absolute throughout. The first criterion was originally meant to hold from the start, but compartments filling from 0 in the first months carry the two integrators' first-step error (see SUMMARY).

- **A3. The calibration.** Started at `TBODE`'s equilibrium, does the v4 ODE stay there for 200 years? And does calibrating `beta` to 500/100k on the v4 ODE alone (equilibrium from a long run) give 9.9173?

### B. How should agents be set up in v4?

- **B1. Network degree.** At 50,000 agents, the mean degree, degree variance and edges per agent of `ss.RandomNet(n_contacts=10)`, and whether edges are replaced every step (the default `dur=0`). In v0.10 (starsim 3.6.1), `n_contacts` was the mean degree.
- **B2. Per-contact `beta`.**
  - Run the network sim with `solver='ode'`, which v4 approximates by homogeneous mixing, and check whether a per-contact `beta` of `beta_ode / K` reproduces the plain ODE.
  - As agents (monthly `dt`, 50k agents, started at the ODE equilibrium, 50 years, 4 seeds), compare the infections in each step with what the ODE's force of infection predicts from the state at the start of that step. There are three arms:
    - homogeneous mixing;
    - the network with v4's own mapping, `beta_ode / K` per contact;
    - the network with v0.10's matched mapping, `beta_edge = −ln(1 − beta_ode·dt/K)/dt` (`tbode4.beta_edge_from_ode`). It makes K·(1 − e^(−`beta_edge`·`dt`)) = `beta_ode`·`dt`, so each step's expected infections match the ODE's.

    The third arm was added after the second came out 4% low (see SUMMARY). Under homogeneous mixing, the step's transmission probability is expected to be `1 − e^(−x)`, x = `beta`·w·`dt`, which is a fraction of about x/2 below the ODE's flow (exp 10). Here that is ~0.2%, so this checks the mapping, not that bias.

- **B3. Starting at the equilibrium.** `init=` takes the share in each state. Check that agents start with those shares, so exp 03's `init_from_ode` isn't needed.
- **B4. Notes and run time.** The notes v4 prints about `dt` for agents at monthly, weekly and daily steps, and the time per step at 50k agents, under homogeneous mixing and on the network. `sim.check_dt()` at monthly `dt` (homogeneous mixing, 50 years), which reruns the sim at half the step and reports the largest change in each result.

## How to run

```bash
.venv/bin/python experiments/12_v4_setup/export_reference.py      # v0.10 environment: TBODE references
v4/.venv/bin/python experiments/12_v4_setup/run.py                 # v4 environment
```

`v4/tests/test_anchor.py` keeps A1 and A2 as regression tests, to rerun whenever the v4 pins are bumped.
