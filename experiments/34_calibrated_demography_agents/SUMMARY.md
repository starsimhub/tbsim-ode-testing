# 34: Summary — agents against the ODE with ages and demography, across exp 32's posterior

**Track:** v4 (tbsim `a2ae583`, starsim `0089ac0f`, 4.0.0rc5). Run on zebra (100 workers) from local commit `34934ec` plus this experiment's `run.py`. Zebra's checkout is partial, so `versions.json` names zebra's old commit.

## Answer

**Yes.** At 12 sets drawn from exp 32's posterior, v4's agents reproduce v4's ODE with ages and India's demography from 1950 to 2040, at weekly steps within noise:

- all-ages prevalence is within 0.4% of the ODE on average in every window;
- ratios per prevalent case are within 0.25%;
- prevalence by age group is within 0.6% on average.

At monthly steps agents run up to 0.8% high on average after 1995. That is the known slow settling of monthly agents. Agents' values of exp 32's targets are within about 0.5% of the ODE's on average. A set that fits the targets as an ODE fits them as agents too.

## What was run

- **Model:** `tbode4.TBTreat`, the model exp 32 calibrated.
  - **Natural history:** the two latent states, with every reinfection risk at π.
  - **Mixing and treatment:** homogeneous mixing, and a fixed 6-month treatment course.
  - **Time-varying inputs:** coverage, efficacy and `sym_dead`, as exp 24's schedules.
  - **Demography:** India's WPP fertility and mortality, with Starsim's defaults (exp 33).
- **Sets:** 12 by weighted resampling without replacement, 6 from each of exp 32's IMIS seeds (`sets.log`).
  - They span `k_latent` 3.5–8.5 and `inf_asy` 0.13–0.36 along the posterior's ridge, and 15+ prevalence of 391–889 per 100,000 in 1995.
  - **Set 0 doesn't fit:** an initial prior draw with log-likelihood −25, against a best of 3.1. It falls outside its seed's 95% intervals on 7 parameters and outside 3 targets' intervals.
  - **Why it was drawn:** exp 32's IMIS gave it a weight of 0.003, close to the largest of any point. A badly fitting point from the initial prior sample gets a large weight wherever the mixture's density is low. About 1% of each seed's weight sits on such points (log-likelihood more than 20 below the best).
  - **Kept:** it was kept as the README's draw. It is a valid check of agents against the ODE, but not a calibrated set.
- **Start:** each set's settled start by single year of age and sex.
  - **How it's made:** 100 years at 1950's demography from the untreated equilibrium (`fastode_age.settle()`, `start_table()`), as the calibration starts.
  - **Used by:** the ODE through `init=` tables, and agents with drawn counts (`tbode4.use_drawn_start()`).
  - **Settling:** the starting shares were still moving at the end. The largest change in any group's share of any state over the last 20 years was 0.03–2.8 percentage points, largest for set 10; both solvers start from the same table, so this doesn't affect the comparison.
- **Runs:**
  - **ODE:** v4's ODE with quarter- and eighth-year age groups, extrapolated to exact aging (2 × eighth − quarter).
  - **Agents:** 100,000 at the start (up to 471,000 by 2040), 20 reps per set and step, at monthly and weekly steps.
  - **Run time:** agent runs took 133 s monthly and 479 s weekly.
- **Counting agents by age:** this experiment's `AgeCounts` leaves out agents dying in the step. Analyzers run before Starsim removes the dead, so exps 30 and 33's version counted each step's deaths as alive at its end. With the fix, the counts equal Starsim's `alive` exactly.

## Mapping check

`mapping.log`, `mapping.csv`. The check compares v4's ODE (monthly, exact aging) with `fastode_age.run_cohort()` from the same table at each set, in yearly means 1951–2040.

- **Prevalence, mortality, notifications and the population:** within 0.09% (prevalence by age group within 0.17%).
- **Incidence:** within 0.25%. The largest gap comes in the mid-1950s; it is probably the ODE's first-entry count, which gives every non-infectious person the same background death rate (exp 30's `ode_incidence`).
- **The eight targets:** within 0.07%.
- **Conclusion:** the mapping from exp 32's parameters to `TBTreat` is right.

## Agents ÷ ODE over time

`summary.log`, `by_set.csv`, `pooled.csv`, `agents_vs_ode.png`, `agents_vs_ode_ages.png`, `over_time_{month,week}.png`.

Agents ÷ ODE (1.00 = match): the mean over 12 sets of each set's mean over 20 reps ± the SE between sets; windows of yearly means. Per case: ÷ all-ages prevalence.

| Output | Step | 1950–1994 | 1995–2009 | 2010–2024 | 2025–2040 |
|---|---|---|---|---|---|
| Prevalence, all ages | weekly | 0.996 ± 0.002 | 0.998 ± 0.002 | 1.003 ± 0.001 | 1.004 ± 0.002 |
| | monthly | 1.001 ± 0.002 | 1.005 ± 0.002 | 1.007 ± 0.003 | 1.008 ± 0.003 |
| Prevalence, 15+ | weekly | 0.996 ± 0.002 | 0.998 ± 0.002 | 1.003 ± 0.001 | 1.003 ± 0.002 |
| | monthly | 1.001 ± 0.002 | 1.005 ± 0.002 | 1.008 ± 0.003 | 1.009 ± 0.003 |
| Share asymptomatic | weekly | 1.000 | 1.000 | 1.000 | 1.000 |
| | monthly | 0.999 | 0.998 | 0.999 | 0.998 |
| Incidence per case | weekly | 0.999 | 1.000 | 1.000 | 1.000 |
| | monthly | 0.997 | 0.998 | 0.998 | 0.999 |
| Notifications per case | weekly | — | 0.998 ± 0.001 | 0.998 ± 0.001 | 1.000 ± 0.001 |
| | monthly | — | 0.992 ± 0.002 | 0.997 ± 0.001 | 0.999 ± 0.001 |
| TB mortality per case | weekly | 0.999 | 1.000 ± 0.001 | 0.999 ± 0.001 | 0.999 ± 0.002 |
| | monthly | 0.998 ± 0.001 | 0.997 ± 0.001 | 1.001 ± 0.002 | 0.996 ± 0.002 |
| Population | weekly | 0.999 | 0.997 ± 0.001 | 0.997 ± 0.001 | 0.996 ± 0.001 |
| | monthly | 1.000 | 1.001 ± 0.001 | 1.001 ± 0.001 | 1.001 ± 0.001 |

SEs below 0.001 are left out.

- **Weekly steps, set by set:** all-ages prevalence is beyond 2 SE of the ODE in 2 of 48 set-windows (z = 2.0 and 2.3), about what chance gives (2.4). Each set's SE is about 0.7% of prevalence.
- **Monthly steps, set by set:** 8 of 48 set-windows are beyond 2 SE, every one with agents high, all after 1995 and growing over time. The largest are sets 4 and 8 in 2025–2040, at 1.023 (z 2.9 and 3.2).
  - **The step causes it:** the same sets match at weekly steps (1.006 and 1.002). Monthly minus weekly averages +0.44% over the sets in 2025–2040.
  - **No parameter explains it:** none correlates with the monthly ratio by more than |r| = 0.41 across the 12 sets.
  - **What it is:** the slow settling of monthly agents seen in exps 27 and 33 (about 1% over decades). The README's follow-up, rerunning outlying sets at weekly steps, is already answered.
- **Monthly notifications per case in 1995–2009** (0.992) is the known half-step lead of time-varying rates (treatment starts low in its first years; exp 27).
- **Prevalence by age group** is within 0.6% on average at weekly steps and within 1.1% at monthly steps in every window.
  - **By set:** individual sets range 0.96–1.04. Each set's SE is 1–2% (20 reps), so 20 reps can't test the README's 2%-per-set criterion more finely.
  - **Beyond 2 SE:** at weekly steps, 14 of the 240 set-window-group ratios (6%), against about 11 by chance. At monthly steps, 28 (12%), mostly agents high after 1995, as for all-ages prevalence.
- **Population:** within 0.4% at weekly steps, as in exp 33.
  - **65+ share at monthly steps:** 0.993–0.997, but the ODE's quarter-year groups alone put it 0.4–0.6% higher than exact aging.
- **Recent latent infection** reads 35% high at monthly steps and 7% at weekly steps (end-of-step snapshots, known). Remote latent infection and the cleared share are within 0.5%.
- **Deaths on failing courses per case:** 0.99–1.01, within noise (SE about 0.6%).

## Targets

`targets.csv`, `targets_pooled.csv`, `targets.png`.

Agents ÷ ODE for exp 32's eight run-based targets: the mean over sets ± SE between sets.

| Target | Monthly | Weekly | Sets inside the 95% interval, ODE / agents (monthly) |
|---|---|---|---|
| 15+ prevalence 2020 | 1.005 ± 0.003 | 1.006 ± 0.003 | 11 / 10 |
| 15+ share asymptomatic 2020 | 1.000 ± 0.002 | 1.002 ± 0.001 | 12 / 11 |
| 15+ prevalence ÷ incidence 2020 | 1.005 ± 0.004 | 1.008 ± 0.002 | 12 / 12 |
| Mortality ÷ incidence 2024 | 1.003 ± 0.005 | 0.999 ± 0.010 | 11 / 11 |
| 15+ prevalence ÷ notifications 2020 | 1.015 ± 0.004 | 1.004 ± 0.004 | 12 / 12 |
| Incidence 2024 ÷ 2015 | 0.996 ± 0.004 | 1.004 ± 0.003 | 12 / 12 |
| Mortality 2024 ÷ 2015 | 1.013 ± 0.008 | 1.016 ± 0.010 | 10 / 10 |
| 15+ prevalence 2040 ÷ 2020 | 1.009 ± 0.006 | 1.001 ± 0.004 | 11 / 11 |

- **Averages:** within about 0.5%, except notifications (1.5%, monthly only) and the mortality ratio (1.3–1.6%).
- **Noise:** each set's values are noisy for these one-year quantities (SE 1–3% per set with 20 reps), so per-set ratios range 0.93–1.08.
  - **Beyond 2 SE:** 8 of 96 set-target pairs at monthly steps and 4 at weekly, against about 4 by chance.
  - **Monthly prevalence ÷ notifications** is high because monthly agents' prevalence is high and their 2020 notifications slightly low.
- **Inside or outside the target intervals:** agents and the ODE disagree only where the ODE sits at an interval's edge, in 5 of 192 cases.
  - Set 7's prevalence sits at the upper edge.
  - Set 0's monthly share asymptomatic and weekly mortality ratio sit at a lower edge.
  - Set 8's weekly mortality ÷ incidence sits at the upper edge.
- **Set 0:** the ODE misses 3 targets, as its log-likelihood says.

## Pitfall checks

`pitfall_checks.csv`, pooled over sets and reps:

- **Course length:** courses last 0.500 years at both steps.
- **Starts from the wrong state:** no starts from outside symptomatic TB, and no relapse except from cure.
- **Rare counts:** a start while on treatment in 104 of 3.2 million starts (monthly; 10 weekly), and a TB death from a successful course in 29 (5). These are the same rare counts as exps 30 and 33, near 3 per 100,000 starts.

## Caveats

- **Twelve sets, 20 reps:** these resolve about 0.7% in prevalence per set and window. Age groups and one-year target values are noisier (1–3%).
- **Set 0** is a poorly fitting prior draw that IMIS weighted heavily. Exp 32's posterior summaries include about 1% weight on such points.
- **The settled start** hadn't settled fully at some sets (up to 2.8 percentage points' change in a group's share over 20 years). This matters for the calibration's start, not for this comparison.

## Files

- **The sets:** `run.py`; `sets.obj` and `sets.log` (the sets, their starting tables and `fastode_age` runs).
- **The mapping check:** `mapping.csv`, `mapping.log`.
- **The runs:** `ode.obj` and `ode_parts/` (ODE runs, yearly values and targets); `runs.obj` and `runs_by_group/` (per-rep yearly values, target values and checks).
- **Summaries:** `by_set.csv`, `pooled.csv`, `targets.csv`, `targets_pooled.csv`, `ode_widths.csv`, `pitfall_checks.csv`, `summary.log`.
- **Logs and versions:** `check.out`, `runs.out`, `versions.json`.
