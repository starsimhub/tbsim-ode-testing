# 27: Do v4's agents match its ODE with treatment, time-varying TB mortality and over-detection, across the calibrated posterior?

**Status:** done; see [SUMMARY.md](SUMMARY.md). Design confirmed by the user (2026-11-10). Added after the main run, to explain two monthly-step departures: `ramp_timing.py` (a time-varying rate in agents against the ODE) and `lag_effect.py` (the ODE with its time-varying inputs delayed by half a month).
**Track:** v4 (`v4/.venv`; tbsim `origin/v4` @ `619c296`, starsim `origin/v4` @ `d83f4f87`).

## Why this experiment

Exp 23, rerun at the current pins as exp 26, showed that the treatment design works in agents at one hand-picked parameter set, with constant TB mortality. This experiment asks whether that holds for the model as it is calibrated, the one that would carry into a new project. It is the treated counterpart of exp 21, which did the same check without treatment.

It reuses exp 23's machinery: `TBTreat` with a fixed course, agents ÷ v4's ODE at monthly and weekly steps, the same outputs, the pitfall checks, and a run from the untreated equilibrium through the scale-up. What it adds:

| | Exp 23 (and exp 26) | Exp 27 |
|---|---|---|
| Parameters | One hand-picked natural history (exp 21's set 1), not calibrated with treatment | 12 sets drawn from exp 24's calibrated posterior, covering both ends of its ridge |
| TB mortality (`sym_dead`) | Constant | Falls from 1995 to 2024, as calibrated |
| Over-detection | None | Included (affects notifications only) |
| Treatment over time | tb_LAI_TPT's best-fit coverage and efficacy | Each set's own calibrated coverage, efficacy and relapse |
| Scenarios | Scale-up, and constant treatment at equilibrium | Scale-up only, read in three windows (1995–2009, 2010–2024, 2025–2040) |
| Against the calibration targets | No | Yes: would agents fit them as the ODE does? |

## Background

- **Treatment in agents:** [exp 23](../23_treatment/SUMMARY.md) checked that agents match the ODE with exp 23's treatment design (CONVENTIONS.md 9), implemented as `tbode4.TBTreat`. [Exp 26](../26_v4_pin_bump/SUMMARY.md) repeated the check at the current pins.
  - People with symptomatic TB start treatment at θ = −ln(1 − coverage).
  - Whether a course will cure is drawn when it starts.
  - Failing courses face TB death, then return to symptomatic TB.
  - A share of cures relapse.

  Both experiments used one hand-picked natural history (exp 21's set 1), with constant TB mortality.
- **What exp 24 added:** [exp 24](../24_treatment_calibration/SUMMARY.md) calibrated the ODE with treatment to tb_LAI_TPT's 10 applicable targets. Two of its features haven't been checked in agents:
  - **TB mortality that falls over time:** `sym_dead` is constant until 1995, then falls linearly to 2024. It applies to symptomatic TB and to failing courses alike.
  - **Over-detection:** reported notifications are (1 + `over_detection`) × treatment starts. This affects only the output compared with data, not the dynamics.
- **How exp 24's ODE differs from agents:** exp 24 calibrated `tbode4.fastode`, a standalone copy of the ODE that runs the fixed 6-month course as 12 Erlang stages. v4's own ODE runs the course exactly, as a conveyor of whole steps. The two agree within 0.04% of prevalence over 1995–2040 at one set (exp 24).
- **The posterior:** its predictions of burden over time agree between exp 24's two IMIS seeds. Its parameters don't, along a ridge where `k_latent` (how fast recent infection becomes remote) and `inf_asy` vary together. That ridge isn't being resolved here (user, 2026-11-10). For this check, it only matters that the sets cover both stretches of it.
- **Exp 21** did the same check without treatment, at 12 sets from exp 20's posterior: agents matched within 0.3–1%.

## Question

At parameter sets drawn from exp 24's posterior, do v4's agents reproduce v4's ODE through tb_LAI_TPT's treatment scale-up and the fall in TB mortality?

1. **Over time:** prevalence, the share asymptomatic, and incidence, notifications and TB mortality per prevalent case, over 1995–2040.
2. **The targets:** each of exp 24's targets that comes from a run, with tb_LAI_TPT's definitions. Would agents fit them as the ODE does?

## Design

**Parameter sets:** 12 draws from exp 24's posterior, by weighted resampling without replacement, so each is a distinct point.

- **6 from each IMIS seed** (`real.obj`, `real_seed2.obj`), so that both stretches of the ridge are covered: fast front-loading with a high initial progression rate, and slow with a low one.
- **Each set counts equally.** The seed is fixed (27) and the sets are saved, for reuse.
- **Checked before running:** the table of drawn sets against each seed's 95% intervals.

**Model:**

- **Agents and the reference ODE:** `tbode4.TBTreat`, with a fixed 6-month course, exp 24's natural history (the two-state latent model, every reinfection risk at π), homogeneous mixing, and births replacing every death (μ = 1/70).
- **Time-varying inputs:** coverage, efficacy and `sym_dead`, as `ss.timeseries`, built with exp 24's knots and linear interpolation (`calib.schedules`).
- **Relapse:** a share `p_relapse` of cures, at mean 1.5 years.
- **The reference** is v4's own ODE of this model.
- **Before running, a mapping check:** for each set, v4's ODE against exp 24's `fastode` trajectory. Prevalence should agree within about 0.1%, the difference between Erlang stages and the exact course. Any larger gap is a mapping error.

**Runs for each set:**

- **Start:** 1980, at the set's untreated equilibrium (`fastode`). Everything is constant until the end of 1994, so the 15 years before treatment let agents settle from the shared starting state. Treatment starts on tb_LAI_TPT's first knot (31 December 1994). The run ends in 2041.
- **The ODE:** v4's ODE at monthly and weekly steps. First-entry incidence comes from `tbode4.ode_incidence`.
- **Agents:** 200,000 per run, at monthly and weekly steps, 20 reps each, with new seeds. Incidence is counted with `tbode4.Incidence`, and the pitfall checks with `tbode4.TxChecks`.
- **Why 200,000:** prevalence falls from about 530 per 100,000 in 1994 to about 200 by 2040, so this keeps 400–1,000 people with TB. 20 reps then give an SE of about 1–2% of prevalence for one set.
- **Over-detection:** applied to both sides' notifications in the same way. It cancels in agents ÷ ODE, so it is checked only through the targets.
- **Size:** 12 sets × 2 steps × 20 reps = 480 agent runs, about 30 minutes on 4 cores (exp 26's timings).
- **Flows:** dated by Starsim's current indexing, where index i is the step from t_i to t_(i+1) (exp 26).

**Measured:**

- **Agents ÷ ODE over time:** yearly means, then means over three windows: early scale-up (1995–2009), late scale-up with falling mortality (2010–2024), and constant treatment (2025–2040). For each set and step, mean ± SE across reps; then averaged over the 12 sets, with SE.
- **Outputs:**
  - active prevalence and the share asymptomatic;
  - per prevalent case: incidence (first entry to asymptomatic TB per episode, plus relapses), notifications, the number on treatment, and TB mortality (from symptomatic TB and on failing courses);
  - TB deaths on treatment, reported separately;
  - recent and remote latent infection, and the cleared share.
- **The targets:** each set's values for the 8 targets that come from a run, for agents (per rep) and the ODE, with tb_LAI_TPT's definitions (exp 24's README).
  - **The 8 targets:** prevalence and the share asymptomatic on 1 January 2020; prevalence ÷ incidence and prevalence ÷ notifications in 2020; mortality ÷ incidence in 2024; incidence and mortality 2024 ÷ 2015; prevalence 2040 ÷ 2020.
  - **Compared:** agents ÷ ODE, and each against the target's 95% interval.
  - **Left out:** the two progression targets, which come from the latent chain, not a run.
- **Expected:**
  - **Recent latent infection reads high** in monthly end-of-step snapshots (×1.2–1.7 in exp 21). It's reported, but not counted as a mismatch.
  - **TB deaths on treatment read about 0.4% low** in agents at monthly steps (exp 26, unexplained, not pursued). That moves TB mortality by under 0.1%.
- **The pitfall checks** from exp 23, pooled over sets, are run as a guard.

**Figures:**

- agents ÷ ODE by set and window, for each output, at both steps, with each set's position on the ridge (`k_latent`, `inf_asy`) shown alongside;
- prevalence, incidence, notifications and mortality over 1995–2040 for a few sets: agents' reps, their mean and the ODE;
- agents' target values against the target intervals, next to the ODE's.

The raw outputs behind each figure are saved with `sc.save`: the drawn sets, the ODE runs, and per-rep yearly values and target values.

## What would answer the question

- **Agents match the ODE across the posterior** if every set's agents ÷ ODE is within about 2 SE of 1 in each window, and the average over sets is within about 1%.
- **Ratios per prevalent case** should be within 0.5%, as in exps 23 and 26.
- **If some sets or windows don't match:**
  - check what sets them apart (position on the ridge, fast mortality decline, high coverage, large φ) before calling it a mismatch;
  - rerun those sets with more reps, and at weekly steps if only monthly differs, to separate step effects from noise.
- **A time-varying effect would show** as a ratio that drifts across the three windows while the constant window matches.

## Not covered

- **Model structure:** ages, demography, child transmission and the household network.
- **Treatment outcomes:** loss to follow-up and death on treatment beyond TB death on failing courses.
- **The posterior:** its ridge, and its tails beyond what 12 draws reach.

## Files

- `run.py`:
  - `sets` draws the 12 sets → `sets.obj`, `sets.log` (against each seed's 95% intervals);
  - `check` runs the mapping check → `mapping.csv`, `mapping.log`;
  - with no argument, it runs everything → `ode.obj`, `runs.obj`, `runs_by_group/`, `by_set.csv`, `pooled.csv`, `targets.csv`, `targets_pooled.csv`, `pitfall_checks.csv`, `agents_vs_ode.png`, `over_time_month.png`, `over_time_week.png`, `targets.png`, `run.log`, `versions.json`;
  - `summary` redraws the summaries and figures.
- Uses exp 24's `calib.py` (the parameter mapping, schedules and targets) and its `run.py` (target intervals), and `tbode4`'s `TBTreat`, `Incidence`, `TxChecks` and `ode_incidence`.
- `SUMMARY.md`: findings, written after running.
