# 23: Summary

**Track:** v4 (tbsim `dcfb522`, starsim `80c81f53`; `versions.json`).

**Setting:**

- **The treatment design:** CONVENTIONS.md 9, implemented as `tbode4.TBTreat`, a `tbsim.TB` subclass that runs as v4's ODE or as agents, and in `tbode4.fastode`.
  - People with symptomatic TB start treatment at θ = −ln(1 − coverage).
  - Whether a course will cure is drawn when it starts (adherence × efficacy).
  - Failing courses face TB death, then return to symptomatic TB.
  - A share `p_relapse` of cures relapse, at a rate with mean 1.5 years.
  - Incidence is the first entry to asymptomatic TB per episode, plus relapses.
- **Parameters:** exp 21's set 1 natural history, with tb_LAI_TPT's best-fit treatment values.
- **Scenarios:** 2024's treatment values held constant (`eq`, compared over 2005–2049), and tb_LAI_TPT's scale-up from 1995 (`ramp`, 1980–2035, compared over 1995–2035).
- **Agents:** 200,000, 40 reps per course, scenario and step (monthly, weekly).
- **ODE:** v4's own ODE of the same model, with first-entry incidence rebuilt from its transitions (`tbode4.ode_incidence`). `v4/tests/test_treat.py` checks it against `fastode`.
- **Caveat:** the treated equilibrium's prevalence is low, about 80 per 100,000, so agent results are noisier than in earlier experiments.

## Question 1: the design behaves as intended in agents

Every pitfall check passes (`pitfall_checks.csv`; about 100,000 courses per arm).

| Check | Result |
|---|---|
| Starts from anyone not symptomatic | 0 |
| Starts while on treatment; TB deaths from courses that will cure | 0 real cases. The checks flagged 12 agents over all arms. Each was cured, relapsed and then restarted or died, all within one step (`diagnose.py`), which continuous time allows (about one per 50,000 cures) |
| Relapse from anyone not cured | 0 |
| Restarts in the step of failure (fixed course; expected 1 − e^(−θ·dt/2)) | Monthly 4.44% (expected 4.60%); weekly 1.15% (1.08%) |
| The success share by year of starting, against efficacy then | Tracks it (e.g. 2005: 0.494 against 0.481; 2024: 0.537 against 0.530; about 2,000–3,700 courses a year) |
| Course length, fixed | Monthly: exactly 6 steps for every course. Weekly: 26 or 27 weeks (93% / 7%), mean 0.5000 years |
| Background death of people on courses that will cure, against μ = 0.0143 | 0.013–0.015 per year |
| Population | Constant within ±1% |

## Question 2: agents match v4's ODE with treatment

| Agents ÷ ODE (40 reps per arm, mean ± SE) | Prevalence | Incidence per prevalent case | Notifications per prevalent case | TB mortality per prevalent case |
|---|---|---|---|---|
| Scale-up, pooled over the 4 arms | 1.006 ± 0.004 | 0.998–1.002 | 0.995–1.007 | 0.993–0.998 |
| Equilibrium, pooled over the 4 arms | 0.983 ± 0.008 | 0.998–1.005 | 0.999–1.004 | 0.996–1.010 |

Each arm's SEs are about 0.002 for incidence, 0.004 for notifications and 0.004 for mortality per prevalent case. Per arm, prevalence is within 2.3 SE of the ODE (0.962–1.015; `agents_vs_ode.csv`, `per_prevalent_case.csv`).

- **Treatment's flows match:** per prevalent case, incidence, notifications, the number on treatment and mortality are all within 1% of the ODE, with either course, at either step.
- **Prevalence:** it matches in the scale-up. At equilibrium it reads 1.7% low (2 SE). The arms scatter on both sides (rate monthly 0.966, fixed monthly 1.015), and treatment's flows per prevalent case match there too. So this looks like noise, or the small-population effect of only about 160 people with TB at 200,000 agents, rather than a difference in how treatment works. More reps would settle it, if it matters.

## Question 3: a fixed 6-month course against a rate

| v4's ODE, per 100,000 at equilibrium (2024's treatment) | Rate | Erlang k = 2 | k = 4 | k = 12 | Fixed |
|---|---|---|---|---|---|
| Prevalence | 79.61 | 79.21 | 78.99 | 78.84 | 78.76 (−1.07%) |
| Incidence | 68.01 | 67.74 | 67.60 | 67.50 | 67.45 (−0.83%) |
| Notifications | 26.22 | 26.05 | 25.96 | 25.90 | 25.86 (−1.37%) |
| On treatment | 12.01 | 12.16 | 12.24 | 12.30 | 12.32 (+2.6%) |
| TB mortality | 11.20 | 11.22 | 11.22 | 11.23 | 11.23 (+0.22%) |
| Share of failing courses ending in TB death | 16.5% | 17.3% | 17.7% | 17.9% | 18.1% |

- **The difference comes from TB death on treatment.** With a fixed course, more people on failing courses die before the course ends (18.1% against 16.5%), so fewer return to symptomatic TB and transmit.
  - Prevalence and incidence are about 1% lower.
  - TB mortality is almost unchanged: deaths move from after a failed course to during it.
  - Without TB death on treatment, the expected difference was order 10⁻⁵.
- **During scale-up the gap grows gradually** (`fixed_vs_rate.png`):
  - prevalence is 0.2% lower in 1998, 0.8% in 2010, 1.2% in 2024 and 1.6% in 2035;
  - TB mortality stays within 0.3%;
  - the number on treatment is 22% higher in 1996, as courses accumulate, and 1.4% higher by 2024.
- **Erlang stages approach the fixed course slowly:** k = 12 gives 90% of the difference, and k = 4 gives 73%. For a fixed course, v4's ODE needs no stages, since it runs the course as a conveyor of whole steps.
- **Weekly against monthly:** the ODE's fixed course is 26 weeks at weekly steps (0.4986 years, against agents' mean of 0.5000), and its equilibrium differs from the monthly one by 0.03%.

## Other findings

- **This natural history isn't calibrated for treatment:** with exp 21's untreated calibration, tb_LAI_TPT's treatment takes prevalence from 312 per 100,000 in 1994 to 99 in 2020, far below the prevalence survey's 310.
  - **Mortality ÷ incidence** at the treated equilibrium is 0.165, against the target of 0.088–0.144.
  - **Next:** recalibrate with treatment.
- **Relapses are 0.7% of incidence** at the treated equilibrium.
- **Run cost:** at 200,000 agents with one thread per process, a run takes about 5 s (monthly) to 20 s (weekly). The whole experiment took 19 minutes on 4 cores. With Starsim's default threading it was twice as slow and used about 10 cores.

## What this means

- **The treatment design can be matched by an ODE exactly,** with a fixed course or a rate, and it has none of the intervention pitfalls checked.
- **A fixed course and a rate differ by about 1%** in prevalence and incidence at equilibrium, and by up to 1.6% during and after scale-up, because of TB death on treatment. That's small next to calibration uncertainty, so either can serve. A fixed course is what tb_LAI_TPT uses, and v4's ODE runs it exactly, so it costs nothing extra.
- **`fastode` approximates a fixed course** with Erlang stages (exactly a rate at k = 1). For calibration, k = 12 gets within 0.1% of prevalence.

## Files

- `v4/tbode4/treat.py` (`TBTreat`, `Incidence`, `TxChecks`, `ode_incidence`, `lai_tpt_treatment`); `v4/tbode4/fastode.py` (treatment added); `v4/tests/test_treat.py`
- `run.py` → `runs.obj`, `runs_by_group/`, `ode.obj`, `agents_vs_ode.csv`, `pitfall_checks.csv`, `success_by_year.csv`, `fixed_vs_rate.csv`, `ode_equilibria.csv`, `ode_fail_death_share.csv`, `over_time_*.png`, `fixed_vs_rate.png`, `agents_vs_ode.png`, `run.log`, `versions.json`
- `diagnose.py` → `diagnose_cases.csv`, `per_prevalent_case.csv`, `diagnose.log`
- `checks/`: the preliminary checks (README there)
