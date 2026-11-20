# 32: Recalibrating the ODE with ages and India's demography

**Status:** run (findings in [SUMMARY.md](SUMMARY.md)).
**Track:** v4 (`v4/.venv`; tbsim `origin/v4` @ `619c296`, starsim `origin/v4` @ `d83f4f87`).
**Run on:** zebra (160 cores, non-spot).

## Why this experiment

[Exp 24](../24_treatment_calibration/SUMMARY.md) calibrated the ODE with treatment to tb_LAI_TPT's targets, with no ages and a constant population (births replacing every death at μ = 1/70). India's demography (WPP fertility and mortality by age and sex from 1950, a growing and ageing population) raises the ODE's burden by 1.2–1.6× at exp 24's posterior sets ([exp 28](../28_age_demography/SUMMARY.md)), so exp 24's posterior doesn't fit the targets once demography is in. At exp 24's best fit, 15+ prevalence in 2020 is 389 per 100,000 with demography, against the target of 310.

[Exp 31](../31_fastode_ages/SUMMARY.md) built `tbode4.fastode_age`, a standalone copy of v4's ODE with ages and demography that is fast enough for IMIS (about 1.1–1.4 s per draw on one core, against about 15 s through Starsim's ODE). This experiment reruns exp 24's calibration on it, with exp 24's checks.

## Question

With ages and India's demography, does the ODE still meet every applicable tb_LAI_TPT target, and does IMIS recover a known truth and converge (on parameters, not only predictions)? How do the posterior and its predictions change from exp 24's?

## Model

`tbode4.fastode_age.run_cohort()`, as [CONVENTIONS.md 8](../../CONVENTIONS.md#8-calibrating-the-ode-a-standalone-copy-at-equilibrium-with-imis) prescribes for calibration with ages:

- **TB and treatment:** exp 24's model (`fastode`'s): the two-state latent model, every reinfection risk at π, homogeneous mixing, every rate the same at every age, treatment as a fixed 6-month course (CONVENTIONS.md 9), tb_LAI_TPT's coverage, efficacy and `sym_dead` trajectories over 1995–2024, and reported notifications = (1 + `over_detection`) × treatment starts.
- **Demography** (CONVENTIONS.md 10): India's WPP 2024 fertility and mortality by age and sex, interpolated between decades; births as female as the 1950 population; TB deaths not replaced.
- **Aging and step:** exact aging by cohort (groups one step wide, moved up one group each step), monthly steps, the force of infection at mid-step (a predictor-corrector). Within 0.1% of exact aging and within 0.06% of weekly steps on every target (exp 31).
- **Each draw:**
  1. **The untreated equilibrium without ages** (`fastode.equilibrium`, births replacing deaths at 1/70), as shares the same at every age. A draw without an endemic equilibrium gets zero likelihood, as in exp 24.
  2. **A 100-year frozen start:** from those shares, 100 years at 1950's demography held fixed, without treatment, with `sym_dead` at its 1995 value (`settle(years=100)`). Its shares by single year of age and sex (`start_table()`) give the 1950 starting state on WPP's 1950 population. 100 years gives the same targets as 200 to 0.002% (exp 31). A draw whose frozen run loses TB (prevalence below 1e-9) gets zero likelihood.
  3. **1950–2041** with WPP's demography and the treatment scale-up from 1995.
- **Cost:** about 1.35 s per draw on one laptop core: 0.03 s equilibrium, 0.7 s frozen start, 0.6 s for 1950–2041.

## Targets and likelihoods

Exp 24's ten targets, means, σ and likelihoods (Binomial for the survey proportions, log-Normal for the rest), with one change: **every target built on the prevalence survey uses 15+ prevalence**, as tb_LAI_TPT computes them (`hm_common.extract_targets_one_sim`). Exp 24 had to use all-ages prevalence for these.

| Target | Model output | Ages | Change from exp 24 |
|---|---|---|---|
| Prevalence | Point prevalence, 1 January 2020 | 15+ | All ages → 15+ |
| Prop_Asymptomatic | Asymptomatic ÷ active TB, 1 January 2020 | 15+ | All ages → 15+ |
| Prev_to_Inc | 15+ prevalence 2020 ÷ incidence in 2020 | 15+ ÷ all ages | Numerator → 15+ |
| Prev_to_Notif | 15+ prevalence 2020 ÷ reported notifications in 2020 | 15+ ÷ all ages | Numerator → 15+ |
| Prev_2040_Ratio | 15+ prevalence, 1 January 2040 ÷ 2020 | 15+ | All ages → 15+ |
| Mort_to_Inc, Inc_Ratio, Mort_Ratio | As exp 24 (TB deaths, first-entry incidence with relapses; CONVENTIONS.md 5) | All ages | None |
| ProgressionY1, ProgressionRatio | As exp 24 (`tbode4.progression_targets`, exact) | — | None |

Incidence, notifications and TB deaths are per person of all ages, as tb_LAI_TPT's are. `fastode_age` gains one output for this, asymptomatic TB at 15+, which doesn't change the dynamics.

**Left out:** Prop_Transmission_HH (no household network). tb_LAI_TPT has no under-15 calibration targets (it dropped Notif_Ratio_u15 and only reports Prop_Incidence_u15).

**Accepted mismatches, for later rows:**

- **Child transmission:** tb_LAI_TPT sets `rel_trans` = 0 under age 10; here children transmit like adults. That's the next row of the consistency table.
- **Children's treatment coverage:** tb_LAI_TPT treats under-15s at 0.22 × adults' coverage; here one coverage applies at every age.

## Parameters and priors

Exp 24's 21 parameters with its priors (its README, "Parameters and priors"), and its parameterisations of the trends. The only change is that background mortality μ is gone: deaths come from WPP. `sym_dead_1995` applies before 1995 too, including in the frozen start.

## Method

As exp 24 (CONVENTIONS.md 8):

1. **Coverage check:** 5,000 prior draws. For each target, the share inside its 95% interval and the prior predictive's 5–95% range; the share of draws inside every interval at once; the share without an endemic state.
2. **Re-identification:** IMIS on synthetic targets, the truth's own outputs from this model, with the real targets' σ and survey sizes.
   - **The truth:** exp 24's (exp 21's set 1 natural history, with tb_LAI_TPT's best-fit treatment, `sym_dead` trend and `over_detection`).
   - **Checks:** every synthetic target fitted; the log-likelihood at the truth no higher than the best point found; which parameters cover the truth in their 95% interval, and for those that don't, whether points far from the truth fit as well (not identified) or not (a failure).
3. **Real targets,** IMIS seed 1.
4. **Real targets, IMIS seed 2:** convergence, compared on **parameter marginals** as well as predictions. For each parameter, both seeds' medians and 95% intervals. A parameter counts as agreeing if each seed's median lies within the other's central 50% interval. Exp 24's two seeds disagreed along the ridge between `k_latent` and `inf_asy` (the shape of early latent progression). The progression targets are unchanged, so the ridge is expected to remain. It is reported, not fixed, here (README Plan: deferred).

**IMIS settings:** exp 24's. 10,000 initial draws, 1,000 per added component, a stop at (1 − 1/e) × 3,000 expected unique points in a resample of 3,000, or after 300 iterations. The algorithm is exp 24's, with one change for speed that doesn't change its results: the mixture density at each point is updated as components are added rather than recomputed over every point and component each iteration. At the end it is recomputed in full and checked against the running value.

**Compute:** zebra, with a persistent pool of 160 single-threaded worker processes (Numba and BLAS at one thread each). 310,000 draws × 1.35 s ÷ 160 ≈ 45 minutes per IMIS run if it reaches the cap, so the whole protocol is about 3 hours. Zebra's checkout is partial (as in exp 30), so the files this experiment uses are copied there from the commit that adds them; `versions.json` records zebra's own commit, so the SUMMARY names the local one.

## Outputs

For each IMIS run:

- the posterior sample with weights, its effective sample size and expected unique points;
- parameter marginals against the priors (and the truth, for re-identification);
- the posterior predictive of each target against its interval;
- trajectories of 200 posterior draws, 1950–2040, by calendar year (stocks on 1 January, flows over the year; `fastode_age.calendar()`): prevalence at 15+ and all ages, incidence, reported notifications, TB mortality and population, against the data;
- a weighted draw of parameter sets, saved for a later check of agents against this ODE at posterior sets.

Across runs:

- **seed 1 against seed 2:** the parameter-marginal table above, and the targets' and trajectories' posterior predictives;
- **against exp 24:** marginals and predictions, to show what demography and the 15+ targets change. Exp 24 ran at earlier pins (tbsim `dcfb522`, starsim `80c81f53`), but `fastode` is standalone and `fastode_age` reproduces its TB model.

## Files

- `calib.py`: the model per draw (`simulate`), the targets, the likelihood (exp 24's), and IMIS with a persistent worker pool. Priors, parameterisations and target means are imported from exp 24's `calib.py`.
- `run.py`: `coverage`, `reid`, `real`, `real_seed2`, `summary` (tables and figures from the saved runs). Each run saves its raw results (`*.obj`) before summarising, and writes `versions.json`.
- `SUMMARY.md`: findings, written after running.
