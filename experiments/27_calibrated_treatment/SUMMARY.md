# 27: Summary

**Track:** v4 (tbsim `619c296`, starsim `d83f4f87`; `versions.json`).

**Setting:**

- **Sets:** 12 parameter sets drawn from exp 24's calibrated posterior, 6 from each IMIS seed (`sets.log`). Between them they span its ridge: `k_latent` (how fast recent infection becomes remote) runs from 2.4 to 11.1, and `inf_asy` from 0.13 to 0.46.
- **Model:** `tbode4.TBTreat` with a fixed 6-month course. Each set has its own coverage, efficacy and relapse over time, and its own fall in `sym_dead` from 1995 to 2024, all as `ss.timeseries`. Over-detection scales notifications.
- **Runs:**
  - start in 1980 at each set's untreated equilibrium; treatment begins at the end of 1994;
  - agents: 200,000 per run, 20 reps per set and step, at monthly and weekly steps;
  - reference: v4's own ODE.
- **Windows:** 1995–2009 (early scale-up), 2010–2024 (late scale-up, falling TB mortality) and 2025–2040 (constant treatment).

**Mapping check:** at every set, v4's ODE reproduces exp 24's calibrated `fastode` trajectories over 1995–2040 (`mapping.csv`):

- within 0.04% of prevalence with the same 12-stage course;
- within 0.08% with the exact fixed course that agents use.

## Finding 1: agents match the ODE across the posterior, within about 1%

| Agents ÷ ODE, mean over 12 sets ± SE between sets | Monthly 1995–2009 | 2010–2024 | 2025–2040 | Weekly 1995–2009 | 2010–2024 | 2025–2040 |
|---|---|---|---|---|---|---|
| Prevalence | 0.999 ± 0.002 | 1.007 ± 0.002 | 1.008 ± 0.004 | 0.999 ± 0.003 | 0.997 ± 0.004 | 1.003 ± 0.005 |
| Share asymptomatic | 0.999 | 0.998 | 0.999 | 0.999 | 1.000 | 1.000 |
| Incidence per prevalent case | 1.000 | 0.999 | 0.998 | 1.000 | 1.000 | 1.000 |
| Notifications per prevalent case | **0.992 ± 0.002** | 0.999 | 0.997 | 0.998 | 1.001 | 0.998 |
| On treatment per prevalent case | **0.993 ± 0.003** | 1.000 | 0.997 | 1.000 | 1.004 | 1.001 |
| TB mortality per prevalent case | 0.999 | 0.997 | 0.997 | 1.000 | 1.002 | 0.997 |
| TB deaths on treatment per prevalent case | 0.990 ± 0.008 | 0.987 ± 0.005 | 0.998 ± 0.007 | 1.003 ± 0.008 | 0.999 ± 0.003 | 0.985 ± 0.006 |
| Remote latent | 0.997 | 1.003 | 1.005 | 0.998 | 0.997 | 1.002 |
| Recent latent (end-of-step snapshot) | 1.37 | 1.38 | 1.38 | 1.08 | 1.07 | 1.08 |

- **The SEs:** 0.001–0.003 between sets for the ratios per prevalent case. Each set's ratio has an SE across its 20 reps of about 1–2% for prevalence.
- **No set stands out:** each set's prevalence is within 2 SE of 1 in every window except one, set 3 at monthly steps over 2025–2040 (1.041, 3.1 SE). There's no clear pattern by position on the ridge, i.e. by `k_latent` (`agents_vs_ode.png`). The source tables are `by_set.csv` and `pooled.csv`.
- **Mostly as expected:** recent latent reads high in snapshots, as in exp 21 (the snapshot effect is noted under Finding 4). TB deaths on treatment are within noise of exp 26's known 0.4% shortfall.
- **Two departures, both at monthly steps only,** are explained below:
  - notifications per case in the early scale-up (Finding 2);
  - prevalence about 0.7% high after 2010 (Finding 3).

## Finding 2: agents apply time-varying rates as if half a step earlier than the ODE

At monthly steps, agents' notifications per prevalent case are well below the ODE's in the first years of treatment, then match from about 2000:

| Year | 1995 | 1996 | 1997 |
|---|---|---|---|
| Agents ÷ ODE | 0.911 | 0.963 | 0.981 |

The number on treatment follows.

**Two checks isolate the cause:**

- **A rate rising linearly from 0 at the end of 1994,** in a one-step model with 4 × 1,000,000 agents (`ramp_timing.py`):
  - monthly: 0.917, 0.968, 0.983 in 1995–1997;
  - weekly: 0.979 in 1995.

  Taking the rate half a step earlier predicts 1 − (step/2) ÷ (mean time since the ramp began): 0.917, 0.972 and 0.983 monthly, and 0.981 weekly.
- **v4's ODE with every time-varying input delayed by half a month** (`lag_effect.py`): notifications per case 0.917, 0.972, 0.983, matching agents.

**So:** agents apply a time-varying rate as if half a step earlier than v4's ODE, whose value applies to the middle of the step.

**Size:** the relative error is (step/2) × (the rate's slope ÷ the rate). So it matters only while a rate rises from near zero, as treatment coverage does from 1995.

- **On mortality:** a falling `sym_dead` changes by about 0.1% per half month, which is negligible.
- **On prevalence:** the half-month delay moves the ODE's prevalence by under 0.15% in any year.
- **At weekly steps:** the effect is a quarter as large.

## Finding 3: monthly prevalence drifts up because agents settle slowly from the ODE's starting state, not because of treatment

At monthly steps, agents' prevalence over time against the ODE's (mean over sets; SE about 0.003):

| Year | 1980 | 1990 | 1994 | 2002 | 2010 | 2020 | 2040 |
|---|---|---|---|---|---|---|---|
| Prevalence, agents ÷ ODE | 0.990 | 0.989 | 0.993 | 0.997 | 1.003 | 1.013 | 1.012 |
| Remote latent, agents ÷ ODE | 0.973 | 0.988 | 0.992 | 0.997 | 1.001 | 1.009 | 1.003 |

**What the pattern shows:**

- **The early deficit predates treatment.** Prevalence is already 1% low in the first year, before any treatment, and remote latent infection starts 2.7% low.
- **The recovery is slow.** Both climb back over decades, the timescale on which remote latent infection turns over.
- **It ends near monthly agents' own level:** about 1.01 of the ODE, as exp 25 found without treatment (1.009 ± 0.005).
- **Weekly steps settle quickly:** they start at 0.996 and are at 1.000 by 1984.

**The cause:**

- **v4 dates infections at the end of their step,** so at monthly steps agents' latent infection in end-of-step snapshots doesn't match the ODE's; that is the recent latent ×1.37 in Finding 1.
- **Their equilibrium sits about 1% above the ODE's.** Starting them at the ODE's state shares puts them off their own equilibrium, mostly in remote latent infection, which settles over decades.
- **15 years before treatment isn't long enough** for that to settle at monthly steps.

**Effect on the trend targets:** the drift shows in the targets that span 2015–2024 (Finding 4).

## Finding 4: agents meet the targets as the ODE does

| Target, agents ÷ ODE (mean over 12 sets ± SE) | Monthly | Weekly | Sets inside the target interval: ODE / agents (monthly) |
|---|---|---|---|
| Prevalence 2020 | 1.012 ± 0.003 | 0.998 ± 0.006 | 7 / 8 |
| Share asymptomatic 2020 | 0.997 ± 0.002 | 1.000 ± 0.001 | 12 / 12 |
| Prevalence ÷ incidence 2020 | 1.006 ± 0.006 | 1.001 ± 0.004 | 12 / 12 |
| Mortality ÷ incidence 2024 | 1.007 ± 0.007 | 0.999 ± 0.013 | 10 / 9 |
| Prevalence ÷ notifications 2020 | 1.006 ± 0.005 | 0.999 ± 0.004 | 10 / 10 |
| Incidence 2024 ÷ 2015 | 1.014 ± 0.006 | 1.004 ± 0.005 | 12 / 12 |
| Mortality 2024 ÷ 2015 | 1.025 ± 0.011 | 1.002 ± 0.015 | 12 / 12 |
| Prevalence 2040 ÷ 2020 | 0.998 ± 0.006 | 1.015 ± 0.005 | 10 / 9 |

Sources: `targets.csv` and `targets_pooled.csv`.

- **The ODE and agents fit the same sets,** to within one set per target, at both steps. The ODE misses some targets at some sets because these are posterior draws, not best fits.
- **The two monthly trend targets read 1.4–2.5% high** in agents (2.2 and 2.4 SE). That is the drift of Finding 3: agents' prevalence rises against the ODE's between 2015 and 2024, and incidence and mortality follow. At weekly steps they're within 0.4%.

## Finding 5: no treatment pitfalls

These are exp 23's checks, pooled over sets and reps (`pitfall_checks.csv`; about 1.5 million courses at each step).

- **What passes:**
  - nobody starts treatment unless symptomatic, and nobody relapses without a cure;
  - every monthly course lasts exactly 6 steps, and weekly courses last 26 or 27 weeks (93% / 7%);
  - background death on treatment is 0.0143 per year, against μ = 0.0143;
  - the population stays within 1.1%.
- **Flagged agents:**
  - at monthly steps, 43 starts while on treatment and 13 TB deaths on courses that will cure (3.7 per 100,000 courses); at weekly steps, 2 and 0;
  - exp 23 traced every such case to a cure, relapse and restart, or a death, within one step. Continuous time allows that. These weren't traced individually here;
  - the rate is higher than exp 26's (2.2 per 100,000 courses), as these sets relapse more (`p_relapse` up to 0.096, against 0.034).

## What this means

- **v4's agents reproduce the calibrated ODE with treatment, falling TB mortality and over-detection**, across the posterior and along its ridge. Ratios per prevalent case are within 0.3%, and prevalence is within 1%.
- **Two monthly-step effects** stop agents matching the ODE more closely. Both are about 1% or less, and both shrink at weekly steps:
  - **Time-varying rates:** agents apply them as if half a step earlier than the ODE. This matters only while a rate rises from near zero, e.g. the first 2–3 years of a treatment scale-up from no coverage.
  - **Starting at the ODE's equilibrium:** that leaves monthly agents slightly off their own. They settle over decades, about +2% in prevalence from 1980 to 2020. A longer pre-treatment period or weekly steps would avoid it, and so would starting agents from their own equilibrium.
- **For calibration targets that are trends** (incidence and mortality 2024 ÷ 2015), monthly agents read 1.5–2.5% above the ODE. That is well inside the targets' intervals (σ about 0.11–0.12 in the log).

## Files

- `run.py`; `sets.obj`, `sets.log`; `mapping.csv`, `mapping.log`; `ode.obj`, `runs.obj`, `runs_by_group/`
- `by_set.csv`, `pooled.csv`, `targets.csv`, `targets_pooled.csv`, `pitfall_checks.csv`
- `agents_vs_ode.png`, `over_time_month.png`, `over_time_week.png`, `targets.png`
- `run.log`, `versions.json`
- `ramp_timing.py` → `ramp_timing_runs.csv`, `ramp_timing.csv`, `ramp_timing.log`: a time-varying rate in agents against the ODE
- `lag_effect.py` → `lag_effect.csv`, `lag_effect.log`: the ODE with time-varying inputs delayed by half a month
