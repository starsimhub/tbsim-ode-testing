# 35: Summary — death and loss to follow-up on treatment

**Track:** v4 (tbsim `a2ae583`, starsim `0089ac0f`).

- **Where it ran:** stages 2–3 and the check on zebra (100 workers, alongside another user's 8–16 processes), from local commit `2550e63` plus this experiment's code; zebra's `versions.json` names its old checkout. Stage 0 ran on the laptop.
- **Run times:** agents about 18 minutes, stage 3 under a minute.
- **How the runs came about:** the design was revised once after review (user, 2026-11-23): all new cases in the outcome data, fitted schedules, treatment of latent infection, TB death by starting state, and a direct check of agents against `fastode_age`. Every stage was then rerun, and the numbers here are from that run.

## Answer

1. **The structure runs identically in both solvers.**
   - **Agents against v4's ODE, at weekly steps:**
     - all-ages prevalence within 0.5% on average in every window;
     - incidence, notifications, the number on treatment and TB mortality per prevalent case within 0.6%;
     - deaths on treatment, losses to follow-up and partial cures within noise;
     - every pitfall check zero, including the new paths for latent infection and for states without TB disease.
   - **Agents against `fastode_age` (the calibration's model) directly, at weekly steps:** prevalence 0.995–1.004 and every ratio per prevalent case within 0.4%.
   - **`fastode_age` against v4's ODE:** within 0.3% on every main output and 0.1% on the targets.
   - **Monthly agents** run up to 1.8% high in prevalence after 1995 against both ODEs: the known monthly settling, with fewer deaths and losses on treatment (1–3%).
2. **It changes the calibrated model little.** With overall cure held at each draw's calibrated efficacy, exp 32's targets move by a few percent (median +1.6% in 2020 15+ prevalence, +2.4% in mortality 2024 ÷ 2015). The log-likelihood is unchanged on average (+0.02), and every target is met by as many draws as before. The carve-out fails for 1.3% of draws in 1995, all with calibrated efficacy above about 0.86. So the posterior stands with the new structure until the reparameterisation.

## What was built

The design is in the README; the code is in `tbode4` (tested in `v4/tests/test_treat.py`; 26 tests pass).

- **`TBTreat`'s new options:**
  - **Course structure:** `tx_from` (the starting states and their rates), `tx_phase_dur` (phase width, or a list of lengths), and pathways as a loop over sets of treatment states (one used).
  - **Death and loss:** `tx_death_fail` and `tx_death_cure` (by phase, and by starting state as {state: rate}; symptomatic TB's otherwise, none for other states with TB disease), `tx_ltfu`, `tx_ltfu_cure` (by phase, or a function of the fraction completed).
  - **Latent infection:** `tx_clear_latent`.
  - **Bookkeeping:** generated state and transition names, listed in `tb.tx_info`.
- **The defaults reproduce the old model step for step:** v4's ODE to 4e-16, agents bit for bit, against reference runs saved before the change (`v4/tests/make_tbtreat_reference.py`).
- **What a course does depends on the starting state:**
  - **TB disease** (symptomatic, asymptomatic, non-infectious): cure and fail tracks, cure to `treated` or `relapsing`, everything else back.
  - **Latent infection** (recent or remote): a clearance track (probability `tx_clear_latent`) ending in `treated`, never `relapsing`, so people count as previously treated; the other track returns to the latent state; no TB death.
  - **Any other state without TB disease:** one track, no cure, no TB death, back to the starting state.
  - **In every case,** natural history pauses during the course, and nobody on treatment is infected.
- **Loss to follow-up on the cure (or clearance) track** is two competing exits, q·h to cure and (1 − q)·h back. That is the same as one exit at h with a split, and it counts partial cures directly.
- **Phases at weekly steps:** v4's ODE runs each fixed length as whole steps, so three 2-month phases became 3 × 9 weeks against agents' exact 8.7.
  - **Effect:** in a cohort, the ODE gave 2–4% more deaths and losses.
  - **Fix:** `tbode4.phase_lengths()` gives both solvers phases of whole steps adding up to the course (9 + 9 + 8 weeks; 2 + 2 + 2 months). With it, every transition matches within 1 SE.
- **`tbode4.outcome_shares()`:** a course's shares cured, dying (TB and background), lost and failed, exact in continuous time. v4's ODE reproduces them within 0.1 percentage points (a test).
- **`tbode4.carve_out()`:** solves for c, μ_f, μ_c,1 and h.
- **`tbode4.fastode_age`:** `tx_options()` for `run_cohort()`. Each slot of the conveyor covers time on treatment [j − ½, j + ½), so it takes each phase's hazards in proportion to its overlap. Assigning whole slots to phases had made phase 1 last 1.5 steps instead of 2.
- **`v4/tests/conftest.py`** sets one numba thread before any test imports Starsim, which otherwise made agent tests fail when the suite ran on its own.

## Stage 0: WHO's outcome data

`data.py`, `outcomes_all_new.csv`, `outcome_schedule.csv`, `outcome_schedule.png`; the options considered are in `outcome_options.png`.

**What's used: all new cases, as shares of evaluated cases.** "Not evaluated" is left out of the denominator.

| Years | Cohorts | Note |
|---|---|---|
| 2012 on | New and relapse | All forms of TB |
| 2004–2011 | New smear-positive + new smear-negative / extrapulmonary / smear unknown | Pooled |
| 1999–2003 | New smear-positive only | Scaled to all new cases by the 2004–2006 ratio: × 0.85 for died, × 1.05 for lost |

- **Why pool:** smear-negative and extrapulmonary cases die less often (3.0–3.3% against 4.1–4.6%), and the model doesn't distinguish the forms of TB.
- **Cohorts before 1999** are left out: under 30% of each was evaluated during the rollout of DOTS.

**The schedules** (user, 2026-11-23):

- **Died:** constant at **3.9%**, the mean of 1999–2023. The year-to-year changes (SD 0.55 points) look like noise.
- **Lost:** linear between knots fitted by least squares to 1999–2023 (residual SD 0.4 points): **20.6%** (1999), **6.9%** (2003), **5.6%** (2014), **1.9%** (2023). It is held at the 1999 value before 1999 and constant after 2023. This follows the fall after the rollout of DOTS, the plateau and the later fall.
- **These are notified, mostly public-sector cohorts** (README).

## The phase width

`phases.log`. This is arithmetic, with no simulation: for a single course on the cure track, at each knot's loss-to-follow-up rate h and a background death rate of 0.01 per year, it compares two ways of computing the share lost and cured with q(f) = f.

- **Continuous time:** someone lost after time t is cured with probability t/D.
- **Phases:** q is fixed at each phase's midpoint.

The step function overstates the cured-on-loss share as follows:

| Phases | 1995 (h = 0.47 per year) | 2023 (h = 0.04 per year) |
|---|---|---|
| One phase (6 months) | 4.2% | 0.4% |
| 2-month phases | 0.5% | 0.05% |
| 1-month phases | 0.1% | 0.01% |

- **Overall cure** is within 0.05 percentage points of continuous time for every width.
- **In stage 3,** the targets with 1-month phases equal those with 2-month phases to 0.1%.

## The check: `fastode_age` against v4's ODE

`check.log`, `check.csv`. Arm A (ω = 0.5, q(f) = f) at stage 2's three sets, monthly, with v4's ODE extrapolated to exact aging. The largest |ratio − 1| over yearly means, 1995–2040:

- **Main outputs** (prevalence, incidence, all-cause TB mortality, notifications, the number on treatment): 0.13%.
- **Deaths on treatment:** 0.29% (cure track 0.36%).
- **Losses:** 0.08%; **partial cures:** 0.12%.
- **Prevalence by age group:** 0.32% (65+, set 10).
- **The eight targets:** 0.10%.
- **Before 1995,** first-entry incidence is within 0.25%, as in exp 34.

## Stage 2: agents against the ODEs

`summary.log`, `by_set.csv`, `by_set_fastode.csv`, `pooled.csv`, `pitfall_checks.csv`, `agents_vs_ode.png`.

- **Sets:** 9, 10 and 2 of exp 34's, along the ridge (`k_latent` 4.2, 6.3, 8.5). Exp 34's set 0, a badly fitting prior draw, was left out.
- **Runs:** 20 reps of 100,000 agents per set, arm and step.
- **The arms:**
  - **Arm A:** treatment from symptomatic TB, with ω = 0.5 (deaths on the cure track, phase 1) and q(f) = f, so every exit is used.
  - **Arm B:** as A, plus placeholder starts from asymptomatic TB (a tenth of symptomatic TB's rate; no TB death on treatment), from recent latent infection (0.01 per year, clearance 0.9), and from susceptible (0.001 per year).

**Agents ÷ v4's ODE.** The mean over 3 sets of each set's mean over 20 reps. Prevalence is given by window (1950–94, 1995–2009, 2010–24, 2025–40); the other rows are ranges over windows.

| Output | Arm A, weekly | Arm A, monthly | Arm B, weekly |
|---|---|---|---|
| Prevalence, all ages | 1.004, 1.001, 0.995, 0.995 | 0.993, 1.008, 1.018, 1.014 | 0.998, 1.003, 1.004, 0.997 |
| Incidence per case | 0.998–1.000 | 0.998–0.999 | 0.999–1.000 |
| Notifications per case | 0.997–1.001 | 0.996–1.000 | 0.995–1.003 |
| TB mortality per case | 0.997–1.002 | 0.996–1.000 | 0.996–1.006 |
| Deaths on treatment per case | 0.992–1.010 | 0.993–1.002 | 0.983–1.008 |
| Losses per case | 1.002–1.019 | 0.986–0.990 | 0.997–1.008 |
| Partial cures per case | 0.99–1.06 | 0.97–1.01 | 0.98–1.00 |
| Prevalence by age group | 0.990–1.013 | 0.992–1.028 | 0.996–1.010 |
| Population | 0.995–0.997 | 1.000–1.001 | 0.995–0.997 |

**Agents ÷ `fastode_age` directly** (arm A, the calibration's model at monthly steps, the same 3 sets):

- **Weekly agents:** prevalence 1.004, 1.002, 0.995, 0.995 by window. Incidence, notifications and TB mortality per prevalent case are within 0.4%, deaths on treatment within 1%, and losses within 2% (± 1.2%).
- **Monthly agents:** the same as against v4's ODE, with prevalence up to 1.018 after 1995.

The comparison is direct, without going through v4's ODE.

**What the comparisons show:**

- **Weekly steps:** both arms match. Set-level misses are no more than chance: at most 1 set beyond 2 SE in any cell, except the population, which is low as in exp 34.
- **Monthly steps:** agents run up to 1.8% high in prevalence after 1995, against exp 34's 0.5–0.8% with the old treatment. They also record 1–3.5% fewer deaths and losses on treatment, as in the cohort test and the known 0.4% gap in exp 26.
  - **What it is:** a step effect from competing rate exits out of fixed-length states in monthly agents, together with slow settling.
  - **With 3 sets** it isn't resolved further (± 0.5–0.7% between sets).
  - **Not pursued.**
- **Pitfall checks,** pooled (`pitfall_checks.csv`):
  - **Zero:** nobody treated without TB disease was cured or died of TB; nobody treated for latent infection relapsed or died of TB; nobody started from symptomatic TB without being symptomatic; and nobody relapsed except after a cure.
  - **Latent clearances happened:** 31,000 in arm B, all to `treated`.
  - **Leaving treatment:** 97% (monthly) and 99.4% (weekly) of those leaving treatment other than by cure or death were in their starting state at the end of the step. The rest had moved on within the step.
  - **Rare counts:** a death from symptomatic TB or the fail track among people on the cure track at the step's start (93 of 760,000 starts monthly, 33 weekly), from leaving and dying in one step. A start while on treatment: 31 monthly, 2 weekly, as exps 30–34.
  - **Outcome shares among agents' courses** (arm A, 1995–2040 pooled): 3.6% died, of which 49–50% on the cure track (ω = 0.5); 4.1–4.2% lost, 32% of them cured on being lost; completed courses 0.500 years (monthly) and 0.499 (weekly).

## Stage 3: the calibrated model

`targets.log`, `targets.csv`, `carve_out.csv`, `targets.png`. Exp 32's 200 posterior draws from each IMIS seed, `fastode_age` at monthly steps. Each is run as the current model and four carve-out variants, each holding overall cure at the draw's calibrated efficacy:

- **base:** ω = 0, q(f) = f;
- **q0:** ω = 0, no partial benefit;
- **front:** ω = 0.5, cure-track deaths in phase 1;
- **fine:** as base, with 1-month phases.

**The targets.** New ÷ current, median (5–95%) over draws; the four variants agree to within 0.002 on every target:

| Target | base |
|---|---|
| 15+ prevalence 2020 | 1.016 (0.982–1.054) |
| 15+ share asymptomatic 2020 | 0.999 (0.994–1.004) |
| 15+ prevalence ÷ incidence 2020 | 1.003 (0.996–1.009) |
| Mortality ÷ incidence 2024 | 1.012 (0.943–1.072) |
| 15+ prevalence ÷ notifications 2020 | 1.000 (0.989–1.006) |
| Incidence 2024 ÷ 2015 | 0.994 (0.988–1.003) |
| Mortality 2024 ÷ 2015 | 1.024 (1.009–1.044) |
| 15+ prevalence 2040 ÷ 2020 | 0.984 (0.969–0.995) |
| Δ log-likelihood | +0.02 (−0.71 to 0.81) |

- **Inside the targets' 95% intervals:** the share of draws is the same as with the current model for every target, except mortality 2024 ÷ 2015, which rises from 0.96 to 0.98.
- **Why the targets move this way:**
  - **The data replace `sym_dead` on failing courses.** WHO's died share sets the fail track's death rate, μ_f: median 0.20 per year in 1995 and 0.21 in 2023, against `sym_dead` of about 0.4 and 0.2. So fewer people on failing courses die in the early years, more return to symptomatic TB, and prevalence is higher.
  - **Deaths on treatment no longer fall with `sym_dead`,** since the died share is constant. So mortality 2024 ÷ 2015 rises, in the direction of exp 32's tension (0.58 against a target of 0.66).
- **The carve-out's solution:**
  - **The cure-track share c** is above calibrated efficacy by about 13% in 1995 (median 0.57 against 0.51) and 1% in 2023, as losses fall.
  - **Capped or unsolved in the base variant:** c is capped at 1 in 0.5% of draws in 1995, and the solve fails in 1.3%, all at calibrated efficacy above about 0.86; 1999 the same; later knots at most 0.3%.
  - **Without partial benefit (q0),** 8.5% fail in 1995, at efficacy above about 0.76. There, deaths only on the fail track can't hold WHO's died share when losses are high and failing courses few. This is the failure of "deaths only among failures" that the README's rationale describes.
  - **LTFU hazard h:** 0.47 (to 1999), 0.15 (2003), 0.12 (2014) and 0.04 per year (2023 on).
  - **Cured on being lost:** 5.8% of courses in 1995 and 0.6% in 2023.

## What this means

- **The structure is ready to use** in both solvers and in `fastode_age`. `fastode_age` now has a direct check against agents too.
- **Exp 32's posterior can stand for now.** With overall cure fixed, deaths and losses move the targets by less than their uncertainty.
- **Partial benefit, cure-track deaths and the phase width** barely matter while overall cure is held fixed. They will matter once c is calibrated instead (the reparameterisation, next).
- **Calibrated efficacy near or above about 0.86** is where the carve-out breaks, so the reparameterisation is needed before efficacy can be near trial values.
- **Treatment from latent infection, asymptomatic TB or false positives** has been tested as plumbing only. Using it needs its own rates and, for latent infection, a trial-based clearance probability.

## Caveats

- **Stage 2 has 3 sets,** so between-set SEs are 0.3–1% in prevalence.
- **The monthly agents' excess** (1–2% after 1995) isn't resolved further.
- **WHO's outcomes are the notified, mostly public sector.** The years before 1999 are held at the 1999 values (20.6% lost), an assumption for the rollout of DOTS.
- **Stage 3 is `fastode_age` only,** at monthly steps, as the calibration runs.

## Files

- **Code:** `data.py` (`main`, `options`, `fitted`); `run.py` (`phases`, `schedules`, `check`, `ode`, `agents`, `summary`, `targets`).
- **Logs:** `outcomes.log`, `outcome_schedule.log`, `phases.log`, `check.log`, `summary.log`, `targets.log`; `runs.out`.
- **Tables:** `outcomes.csv`, `outcomes_all_new.csv`, `outcome_options.csv`, `outcome_schedule.csv`, `phases.csv`, `check.csv`, `by_set.csv`, `by_set_fastode.csv`, `pooled.csv`, `pitfall_checks.csv`, `targets.csv`, `carve_out.csv`.
- **Figures:** `outcomes.png`, `outcome_options.png`, `outcome_schedule.png`, `agents_vs_ode.png`, `targets.png`.
- **Raw outputs (`.obj`, local only):** `sets.obj`, `fastode.obj`, `ode.obj`, `ode_parts/`, `runs.obj`, `runs_by_group/`, `targets.obj`.
- **Data:** `data/` (WHO's CSVs).
- **Versions:** `versions.json`.
