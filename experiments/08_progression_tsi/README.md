# 08 Progression by time since infection: can the ODE match `k_asy`/`k_non` > 0?

## Question

In tbsim, the hazards of progressing from INFECTION to NON_INFECTIOUS and to ASYMPTOMATIC can decline with time since infection τ: `inf_non·e^(−k·τ)` and `inf_asy·e^(−k·τ)`. tb_LAI_TPT's best fit has `k_asy = k_non` = 5.846/yr, so the hazard falls by a factor of e about every 2 months, and progression is strongly front-loaded. This is the first of the two high-priority features. With it added, does the ODE still match the ABM?

As in every feature experiment, two ABM versions are run against the ODE:
- **`TBExact` ÷ ODE: can the ODE match the feature?** `tbode.TBExact` is a version of `tbsim.TB` that samples each agent's end-of-step state from the exact one-step probabilities of the continuous-time model, `expm(Q·dt)`. This removes tbsim's step-ordering bias: tbsim applies transitions one state at a time, so agents moving forward can carry on through several states in a step while agents moving backward stop. At baseline, `TBExact` matched the ODE (exp 05).
- **`tbsim.TB` ÷ `TBExact`: how much does step ordering add?** At baseline, `tbsim.TB` was 0.84 / 0.67 of the ODE prevalence at weekly / monthly `dt` (exps 03 and 05).

## How the ODE and `TBExact` represent the feature

**Poisson-mixture latent classes.** Let r₀ = `inf_non` + `inf_asy`. Without competing risks, a hazard r₀·e^(−kτ) gives a progression time with survival exp(−(r₀/k)(1 − e^(−kτ))). This expands exactly as Σₘ Poisson(m; r₀/k)·e^(−mkτ), which is a mixture of exponentials:
- each new infection joins latent class m with probability Poisson(m; r₀/k);
- class m progresses at the constant rate m·k (class 0 never progresses);
- the destination is NON_INFECTIOUS or ASYMPTOMATIC in the ratio `inf_non` : `inf_asy`.

Clearance and background death act identically on every class, so the mixture is exact with competing risks too. It is not an approximation. The ODE (`TBODE(k_asy=k, k_non=k)`) splits INFECTION into classes INFECTION_0, INFECTION_1, …, truncated where the remaining Poisson weight is below 1e-10: 7 classes for set A, 8 for set B. `TBExact` assigns each agent a class at infection and steps it with that class's constant-rate `expm(Q·dt)`.

[latent_classes.md](latent_classes.md) explains the method in more detail, with worked values for both parameter sets.

A chain of time-since-infection stages was considered and rejected. Agents leave each stage after an exponentially distributed time, which blurs the timing of progression. With 1-month stages, cohort progression is 4–7% too low; even 0.25-month stages are 1–2% off. That's too much near threshold, where prevalence amplifies errors in progression.

**Independent checks** (in `tests/test_progression.py` and `tests/test_exact.py`):
- `TBODE`'s cohort progression matches `ode_cohort_progression` to 1e-4. That function integrates the time-varying hazards directly, without classes.
- `TBExact`'s cohort progression matches it within sampling error (weekly `dt`, 200k agents).
- In one step, `TBExact`'s end-state frequencies match the class-weighted `expm`.
- Standard `tbsim.TB` uses the actual τ (from `ti_infected`), so its daily cohort run gives a further check that doesn't use classes.

**Equilibrium start.** The ABM starts at the ODE equilibrium (`init_from_ode`). Every INFECTION agent also needs a time since infection; otherwise all would look newly infected and progress in a burst. At equilibrium, the inflow to each class and its per-capita exit rate are constant, so time since infection within class m is exponential with rate inflow / stock (`TBODE.latent_age_rates`). Each agent gets a class from the ODE's class proportions and a τ from its class's distribution, stored as a fractional `ti_infected` before the start. `tbsim.TB` reads τ; `TBExact` reads the class.

**Infection timing within the step.** `TBExact` uses `infection_timing='within'` (exp 07): each new infection is given a uniform time within its step and goes through only the rest of the step. This gives correct end-of-step snapshots of INFECTION. It doesn't change how many progress, because the latent class shifts with the infection date. Exp 07's SUMMARY later found a trade-off. The force of infection is taken from the start-of-step state, so `'within'` dates each infection on average half a step after the infectiousness that caused it. That makes the generation interval half a step too long; `'start'` gets it right on average. This doesn't affect the equilibrium comparisons here. Whether to keep `'within'` is an open decision (top-level README, Plan).

## Expected bias in `tbsim.TB`

`tbsim.TB` computes τ at the start of each step and applies that hazard for the whole step. The hazard declines within the step, so this overstates progression; the integrated hazard is about kΔt / (1 − e^(−kΔt)) of the exact value, i.e. 1.26× at monthly `dt` and 1.06× at weekly. In a cohort (prototype, 100k agents), `tbsim.TB` gave about 10–19% more progression than the exact solution at monthly `dt`, 2–5% more at weekly, and about the same at daily. This works against the step-ordering bias, which lowers prevalence. Near threshold (set A), 10–19% more progression acts like a much larger `beta`, so `tbsim.TB` may overshoot the ODE. In a 10-year smoke test (set A, monthly `dt`, one seed), `tbsim.TB` climbed from about 480 to 680/100k while `TBExact` stayed near 500.

## Parameter sets (user decision, 2026-10-12)

In both sets, `k_asy = k_non` = 5.846 (tb_LAI_TPT row 0) and the ODE `beta` is calibrated so the equilibrium active prevalence (ASY + SYM, all ages) is 500/100k. Other settings are as before: constant births and deaths at `mu`, mass action, latent reinfection off, no treatment.

| | A: tbsim defaults + `k` | B: tb_LAI_TPT natural history + `k` |
|---|---|---|
| Natural-history rates | tbsim defaults (as exps 01–07) | tb_LAI_TPT row 0 (`collapsed_samples_eps1.2.csv`) |
| `trans_asymp`, `sym_dead` | 0.82, 0.34 | 0.941, 0.263 (pre-1995 value) |
| Reinfection risk ratios (CLEARED from latent / recovered / treated) | 1.0 / 0.21 / 3.15 | π = 0.242 for all three |
| ODE `beta` (500/100k) | 36.75 | 25.00 |
| Threshold `beta` | 31.7 | 7.59 |
| Prevalence elasticity (Δprev/prev ÷ Δβ/β) | 8.6 | 2.5 |
| ProgressionY1 / ProgressionRatio (targets 0.0523 / 0.663) | 0.011 / 0.586 | 0.051 / 0.674 |

- **A changes only `k` (and the recalibrated `beta`) relative to exps 01–07,** so it isolates the effect of `k` on ABM↔ODE consistency. Its dynamics are unrealistic: progression is a fifth of the target, so `beta` must be about 4× higher. But it stays near threshold, so it remains a sensitive test.
- **B is the reference set from here on:** tb_LAI_TPT's natural history and reinfection protection, matching both progression targets. It sits further above threshold, so prevalence responds less to any transmission bias.

## Design

- **ABM:** 50,000 agents, `ss.RandomNet(n_contacts=10, dur=0)` with the per-contact `beta` from `beta_edge_from_ode`, constant births and deaths. Started at the ODE equilibrium (with τ and classes, as above) and run for 100 years (1900–2000).
- **Arms:** parameter sets {A, B} × models {`TBExact` (`'within'`), `tbsim.TB`} × `dt` {1 week, 1 month}, with 40 reps each (seeds 0–39). The ODE is deterministic: its equilibrium is the reference.
- **Comparison:** each rep's mean over the last 50 years. For each output, report:
  - `TBExact` ÷ ODE and `tbsim.TB` ÷ ODE: mean ± SE across reps, z-score, median and 10–90% range;
  - `tbsim.TB` ÷ `TBExact`;
  - these against the baseline (`k` = 0) values from exps 03, 05 and 07.
- **Outputs:** states (INFECTION, CLEARED, NON_INFECTIOUS, ASYMPTOMATIC, SYMPTOMATIC), prevalence (ASY + SYM), incidence (first ASY or SYM entry per episode, tb_LAI_TPT definition), TB mortality, % asymptomatic.
- **Cohort progression (`run_cohort.py`):** ProgressionY1 and ProgressionRatio for each set × model × `dt`, from `abm_cohort_progression` (everyone latent at t = 0, no transmission, no deaths). 4 seeds × 100,000 agents each, against `ode_cohort_progression`. Adds `tbsim.TB` at daily `dt`, which is the check that doesn't use classes. Because cohort runs observe states only at the end of each step, entries to ASYMPTOMATIC that leave again within the step are missed, so monthly runs read slightly low.
- **Saved data:** `outputs.obj` holds the per-rep time series, the ODE equilibria and the calibration; `cohort.obj` holds the cohort curves; `summary.csv` and `cohort.csv` hold the tables.
- **Diagnostic (`step_map.py`, added after the first results):** a deterministic version of the `TBExact` step with latent classes, run for each set at weekly and monthly `dt`. It applies the expected start-of-step force of infection, then exact transitions with `'within'` timing, and is iterated to equilibrium. This gives what `TBExact` would produce without sampling noise or finite-population effects. Output: `step_map.log`.

## Versions

tbsim main @ `b92d7e8` (v0.10.3), starsim @ `3521adb4` (v3.6.1).
