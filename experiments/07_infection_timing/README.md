# 07 Infection timing within the step: does it remove `TBExact`'s INFECTION artefact?

## Question

`tbode.TBExact` is a version of `tbsim.TB` without tbsim's step-ordering bias. Instead of applying transitions one state at a time, it draws each agent's end-of-step state from the exact one-step probabilities of the continuous-time model, `expm(Q·dt)`, where Q holds the ODE's natural-history rates. In [exp 05](../05_exact_transitions/SUMMARY.md) it matched the ODE at every `dt`, apart from one artefact: at monthly `dt`, the INFECTION stock was 0.89 ± 0.02 of the ODE (0.98 at weekly `dt`, 1.00 at daily).

Exp 05 explained this by how infections are timed. Transmission runs at the start of each step, and newly infected agents then go through that step's full transitions. So every infection is effectively dated to the start of its step, on average half a step earlier than in the ODE, where infections are spread through the step. Observed at the end of the step, each cohort is half a step further into its time in INFECTION. That predicts an INFECTION stock of `exp(−r·dt/2)` times the ODE value, with r ≈ 2.1/yr the exit rate from INFECTION: 0.92 at monthly `dt` and 0.98 at weekly.

If each new infection is instead given a time within the step, and only the remainder of the step is applied to it, does the INFECTION stock match the ODE, with everything else unchanged?

This matters for the next experiment, progression by time since infection. With `k` ≈ 5.85/yr, the progression hazard falls by a factor of e about every 2 months. Dating infections to the start of the step would then make the fastest progressors progress about half a month early. This experiment tests the change on its own, at the baseline, where the expected result is known.

## Change tested

`TBExact` gains an option, `infection_timing`:
- **`'start'` (default, as in exp 05):** newly infected agents go through the step's full transition matrix `expm(Q·dt)`.
- **`'within'`:** each agent infected in this step gets an infection time drawn uniformly within the step, and its end-of-step state is drawn from `expm(Q·(1 − u)·dt)`, where u is the fraction of the step that had already passed when it was infected. This is computed per agent.

Uniform timing is an approximation. With a constant force of infection λ during the step, the exact infection time is a truncated exponential. But λ·dt ≤ 0.4% here, so the mean infection time differs from uniform by less than 0.1% of a step.

Only agents entering INFECTION from an uninfected state count as new infections. Reinfection of latent agents is off, as in all experiments so far. `ti_infected` keeps tbsim's whole-step value; fractional infection times will be added with time-since-infection progression (exp 08).

## Predictions

At `k` = 0, every rate is constant, so infection timing only delays each infected agent's subsequent transitions by about half a step. It doesn't change any outcome probability. So:
- **INFECTION:** from 0.89 to about 1.00 of the ODE at monthly `dt`, and from 0.98 to about 1.00 at weekly.
- **Prevalence and mortality:** no change beyond noise (exp 05: 0.97 ± 0.02 at monthly `dt`, 1.00 at weekly). The force of infection is still fixed at its start-of-step value, which is the most likely cause of the monthly 0.97.
- **Incidence:** no change. Exp 05's monthly 0.94 is thought to come from missed brief visits to ASYMPTOMATIC between end-of-step snapshots, which this change doesn't touch.

## Design

- **Setup:** same as exps 03 and 05. tbsim default natural history, ODE `beta` = 9.917 (500/100k equilibrium), the ABM started at the ODE equilibrium (`init_from_ode`), 100 years (1900–2000), 50,000 agents, `ss.RandomNet(n_contacts=10, dur=0)` with the per-contact `beta` from `beta_edge_from_ode`, constant births and deaths.
- **Arms:**
  - `TBExact(infection_timing='within')`: monthly `dt`, 40 reps (seeds 0–39); weekly `dt`, 160 reps (seeds 0–159).
  - `TBExact` with the default `'start'` timing: monthly, 40 reps (seeds 0–39, from exp 05); weekly, 160 reps (seeds 0–39 from exp 05, seeds 40–159 run here).
  - No daily arm: the predicted INFECTION artefact there is 0.3%, too small to detect.
- **How the runs came about:** the first 40 weekly `'within'` reps came out about 5% low on every output, against the predicted no change. Weekly reps were added for both timings (seeds 40–79, then 80–159, `run_more_reps.py`) to tell noise from a real effect.
- **Reuse of exp 05:** `'start'` seeds 0–39 are reloaded from `../05_exact_transitions/outputs.obj` (same setup and versions). Each script first reruns seed 0 and checks it reproduces the saved time series exactly. `run.py` checks `'start'` monthly; `run_more_reps.py` checks both `'start'` monthly and `'within'` weekly, after `TBExact` gained latent classes for exp 08. This confirms the code changes leave `k` = 0 runs unaffected.
- **Comparison:** each rep's mean over the last 50 years. For each output, report the mean ± SE across reps ÷ the ODE equilibrium (1.00 = match), the z-score, and the median and 10–90% range, for `'within'` and `'start'`. Also report the paired difference, `'within'` − `'start'`, by seed. The seeds share random streams, so the pairs are correlated early in each run but not over 100 years. So this difference is roughly as noisy as an unpaired one.
- **Outputs:** states (INFECTION, CLEARED, NON_INFECTIOUS, ASYMPTOMATIC, SYMPTOMATIC), prevalence (ASY + SYM), incidence (first ASY or SYM entry per episode, tb_LAI_TPT definition), TB mortality, % asymptomatic.
- **Saved data:**
  - `outputs.obj`: per-rep time series for seeds 0–39 at both `dt`s, and the ODE run; table in `summary.csv`.
  - `outputs_more_reps_40-79.obj`, `outputs_more_reps_80-159.obj`: the extra weekly reps; pooled table in `summary_more_reps.csv`.

**Diagnostics for the weekly shortfall** (added after the first results):
- **`run_flux_diag.py`:** per-step flows inside weekly runs (12 seeds per timing, 15 years). Infections from SUSCEPTIBLE and CLEARED are compared with the ODE expectation, together with new infections' end-of-step states and every natural-history flow. Output: `run_flux_diag.log`.
- **`step_map.py`:** a deterministic version of the `TBExact` step, without sampling noise. It applies the start-of-step force of infection, then exact transitions with either timing, and is iterated to equilibrium at daily, weekly and monthly `dt`. Output: `step_map.log`.
- **`run_variants.py`:** two variants of `'within'`, 80 weekly reps each (seeds 0–79). `u0` draws infection times from the same random stream but ignores them, so its timing is the same as `'start'`. `np` applies the infection times but draws them from an independent numpy generator. Outputs: `outputs_variants.obj`, `run_variants.log`.

**Unit test (`tests/test_exact.py`):** in one monthly step, newly infected agents with `'within'` timing end in each state with the probability averaged over a uniform infection time, `∫₀¹ expm(Q·s·dt) ds`. This average is computed independently of `tbode.exact`, as a block of the matrix exponential of an augmented generator. A control test checks that `'start'` fails the same test.

## Versions

tbsim main @ `b92d7e8` (v0.10.3), starsim @ `3521adb4` (v3.6.1).
