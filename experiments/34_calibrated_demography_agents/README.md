# 34: Do v4's agents match its ODE with ages and demography, across the calibrated posterior?

**Status:** designed, not yet run.
**Track:** v4 (`v4/.venv`; tbsim `origin/v4` @ `a2ae583`, starsim `origin/v4` @ `0089ac0f`).
**Run on:** zebra if it's free when the run starts; otherwise the laptop's 4 cores, kept awake (`caffeinate`), about 3–4 hours.

## Why this experiment

[Exp 32](../32_demography_calibration/SUMMARY.md) recalibrated the ODE to tb_LAI_TPT's targets with ages and India's demography, using `tbode4.fastode_age`, a standalone copy of v4's ODE written for speed. Agents have been checked against v4's own ODE with ages and demography at only two hand-picked sets ([exp 30](../30_tb_corrected_demography/SUMMARY.md), and [exp 33](../33_starsim_rc5_demography/SUMMARY.md) at the current pins). Without demography, [exp 27](../27_calibrated_treatment/SUMMARY.md) checked agents across exp 24's posterior and found them within about 1%.

This experiment does exp 27's check for exp 32's posterior: the calibrated model with ages and demography, the one that would carry into a new project.

| | Exp 27 | Exp 33, stage 3 | Exp 34 |
|---|---|---|---|
| Parameters | 12 sets from exp 24's posterior | Exp 30's two sets, from exp 24's posterior, not recalibrated with demography | 12 sets from exp 32's posterior |
| Population | Constant, births replace deaths | India's WPP demography from 1950, by age and sex | The same |
| Start | Each set's untreated equilibrium, 1980 | Exp 30's frozen-start tables by age and sex, 1950 | Each set's frozen-start table from `fastode_age`, 1950, as the calibration starts |
| Reference ODE | v4's ODE | v4's ODE with exact aging (2 × eighth-year − quarter-year groups) | The same |
| Prevalence targets | All ages | — | 15+, as exp 32 |

## Background

- **The model:** `tbode4.TBTreat`.
  - **Natural history:** v4's `tbsim.TB` with the two-state latent model and every reinfection risk at π.
  - **Mixing and rates:** homogeneous mixing, every rate the same at every age.
  - **Treatment:** a fixed 6-month course (CONVENTIONS.md 9), with each set's coverage, efficacy, relapse and `sym_dead` trajectories over 1995–2024, and over-detection of notifications.
  - **Demography:** WPP 2024 fertility and mortality by age and sex, with Starsim's defaults (CONVENTIONS.md 10). Births are female at the 1950 population's fraction at ages 0–4.
- **The calibration's start:**
  - **How each draw starts:** the untreated equilibrium without ages, then 100 years at 1950's demography held fixed, without treatment (`fastode_age.settle()`). Its shares by single year of age and sex (`start_table()`) start the run in 1950.
  - **Agents need the same table**, but also its counts drawn rather than rounded (`tbode4.use_drawn_start()`; proposal 07, still open upstream).
  - **Exp 31's finding:** these tables differ from exp 30's eighth-year frozen ODE runs mostly in the first year of age, and the difference is gone by 1995.
- **Births' sex in the calibration:** exp 32 gave births the 1950 population's overall female share (48.63%). At the current pins, both v4's ODE and `fastode_age` use the fraction at ages 0–4 (49.06%). That moves exp 32's targets by under 0.5% (exp 33). The runs here use the current default.
- **The posterior:** exp 32's two IMIS seeds agree on the targets but on only 5 of 21 parameter marginals. They spread along the ridge between `k_latent` and `inf_asy` (the shape of early latent progression), and differ in starting treatment efficacy and in the burden before 1995 (15+ prevalence about 730 against 550 per 100,000 in 1995). For this check, it only matters that the sets cover that spread.
- **Known effects at monthly steps,** reported but not counted as mismatches:
  - **Treatment's first years:** agents apply time-varying rates as if half a step earlier than the ODE, so treatment starts about 9% low in 1995 (exp 27).
  - **Settling:** monthly agents drift by about 1% from the ODE's starting state over decades (exps 27, 33).
  - **Snapshots:** recent latent infection reads high in end-of-step snapshots.
- **Also known:** agents' 0–14 population share runs about 0.5% below the exact ODE's from 1970 (exp 33, unexplained, within 1%).

## Question

At parameter sets drawn from exp 32's posterior, do v4's agents reproduce v4's ODE with ages and demography, from 1950 through the treatment scale-up to 2040?

1. **Over time:**
   - all-ages prevalence and the share asymptomatic;
   - per prevalent case: incidence, notifications, the number on treatment and TB mortality;
   - prevalence in five age groups;
   - the population's size and age shares.
2. **The targets:** would agents fit exp 32's eight run-based targets as the ODE does? These are with 15+ prevalence and tb_LAI_TPT's definitions.

## Design

**Parameter sets:** 12 draws from exp 32's posterior, by weighted resampling without replacement.

- **Spread:** 6 from each IMIS seed (`real.obj`, `real_seed2.obj`), so the sets span both seeds' parts of the ridge and both levels of pre-1995 burden.
- **Fixed:** the draw's seed is fixed (34), and the sets are saved for reuse.
- **Checked before running:** a table of the drawn sets against each seed's 95% intervals.

**Each set's start:** `fastode_age.settle(years=100)` at monthly steps, then `start_table()`. This gives shares by single year of age and sex, and the share of non-infectious TB already asymptomatic this episode (for first-entry incidence). Both solvers start from it on WPP's 1950 population.

- **The ODE:** `init=` tables by single year of age and sex (`tbode4.init_table()`).
- **Agents:** stratified by single years of age and sex (needed to read the table), with drawn counts.

**Mapping check, before running agents:** for each set, v4's ODE (exact aging extrapolated, monthly) against `fastode_age.run_cohort()` from the same table.

- **Expected:** within about 0.2% in prevalence, incidence, TB mortality and notifications, 1950–2040, as in exp 31.
- **A larger gap** means the mapping from exp 32's parameters to `TBTreat` is wrong. The runs wait until it's found.

**Runs for each set:**

- **The ODE:** v4's ODE with quarter- and eighth-year age groups, at monthly and weekly steps. Exact aging is extrapolated as 2 × eighth − quarter. First-entry incidence comes from exp 30's `ode_incidence`.
- **Agents:** 100,000 at the start (about 480,000 by 2040), at monthly and weekly steps, 20 reps each, new seeds. Counting:
  - incidence with `tbode4.Incidence`, and the pitfall checks with `tbode4.TxChecks`;
  - prevalence by age group with exp 30's `AgeCounts`, extended to count asymptomatic TB, for the 15+ share asymptomatic;
  - population results under Starsim's current names (`death`, `birth`; exp 33).
- **Why 100,000 agents:** 15+ prevalence runs from about 400–1,200 per 100,000 before 1995 to about 140–270 by 2040 across the posterior, so a run has about 300–1,200 people with active TB throughout. 20 reps then give an SE of about 1–2% of prevalence per set and window, as in exps 27 and 30.
- **Years:** 1950–2041, as the calibration.
- **Size:**
  - **Agents:** 12 sets × 2 steps × 20 reps = 480 agent runs. About 40 s each monthly and roughly 2.5 minutes weekly, so about 13 core-hours.
  - **ODE:** 48 ODE runs of 6–45 s.
  - **Wall time:** about 5–10 minutes on zebra's 160 cores, or 3–4 hours on the laptop's 4.

**Measured:**

- **Agents ÷ ODE over time:**
  - **Windows:** yearly means, then means over four windows: before treatment (1950–1994), early scale-up (1995–2009), late scale-up with falling TB mortality (2010–2024), and constant treatment (2025–2040).
  - **Statistics:** for each set and step, mean ± SE across reps; then the mean over the 12 sets, with the SE between sets.
- **Outputs:**
  - all-ages and 15+ prevalence, and the share asymptomatic (all ages and 15+);
  - per prevalent case: incidence (first entry to asymptomatic TB per episode, plus relapses; CONVENTIONS.md 5), notifications, the number on treatment, and TB mortality (from symptomatic TB and on failing courses);
  - prevalence in five age groups (0–14, 15–29, 30–49, 50–64, 65+);
  - recent and remote latent infection, and the cleared share;
  - the population's size and its share in each age group.
- **The targets:**
  - **What:** each set's values for exp 32's eight run-based targets, for agents (per rep) and the ODE: 15+ prevalence and the 15+ share asymptomatic on 1 January 2020; 15+ prevalence ÷ incidence and ÷ reported notifications in 2020; mortality ÷ incidence in 2024; incidence and mortality 2024 ÷ 2015; 15+ prevalence 2040 ÷ 2020.
  - **Compared:** agents ÷ ODE, and each against the target's 95% interval.
  - **Left out:** the two progression targets, which come from the latent chain, not a run.
- **The pitfall checks** from exp 23, pooled over sets, as a guard.

**Figures** (raw outputs saved with `sc.save`: the drawn sets and their starting tables, the ODE runs, per-rep yearly values and target values):

- **By set:** agents ÷ ODE by set and window for each output, at both steps, with each set's position on the ridge (`k_latent`, `inf_asy`) and its pre-1995 burden alongside;
- **Over time:** prevalence (all ages and 15+), incidence, notifications and mortality, 1950–2040, for a few sets: agents' reps, their mean, and the ODE;
- **Targets:** agents' target values against the target intervals, next to the ODE's.

## What would answer the question

- **Agents match across the posterior** if, at weekly steps:
  - every set's all-ages prevalence is within about 2 SE of the ODE in each window;
  - the mean over sets is within about 1%;
  - ratios per prevalent case are within about 0.5% on average (as in exps 27 and 33);
  - prevalence in each broad age group is within about 2% per set (as in exp 30).
- **Monthly steps** are judged the same way, allowing for the known effects above.
- **The targets:** agents' values should be within about 1% of the ODE's, so that a set fitting the targets as the ODE fits them as agents too.
- **If some sets or windows don't match:**
  - **First,** check what sets them apart before calling it a mismatch: position on the ridge, a high pre-1995 burden, fast mortality decline, high coverage, or a large φ.
  - **Then,** rerun those sets with more reps, and at weekly steps if only monthly differs, to separate step effects from noise.

## Not covered

- **Model structure:** child transmission (the next row), age mixing, children's lower treatment coverage, and the household network.
- **Treatment outcomes:** loss to follow-up, and death on treatment beyond TB death on failing courses.
- **The posterior:** its convergence and ridge (exp 32), and its tails beyond 12 draws.
- **Births' sex as calibrated:** exp 32 used 48.63% girls; the runs here use the current default (exp 33: targets within 0.5%).

## Files

- `run.py`:
  - `sets`: the 12 sets and their starting tables → `sets.obj`, `sets.log`;
  - `check`: the mapping check → `mapping.csv`, `mapping.log`;
  - `ode`, `agents`: the runs → `ode.obj`, `runs.obj`;
  - `summary`: tables and figures → `by_set.csv`, `pooled.csv`, `targets.csv`, `targets_pooled.csv`, `pitfall_checks.csv`, `summary.log`, figures; also `versions.json`.
- **Reuses:**
  - exp 24's `calib.py` (the parameter mapping and schedules) and exp 32's (the posterior and targets);
  - exp 30's `TBTreat` setup, `AgeCounts` and `ode_incidence`;
  - exp 33's demography with Starsim's defaults;
  - `tbode4`'s `fastode_age`, `Incidence`, `TxChecks`, `init_table` and `use_drawn_start`.
- `SUMMARY.md`: findings, written after running.
