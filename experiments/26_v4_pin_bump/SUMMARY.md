# 26: Summary

**Track:** v4, at the v4 track's new pins: starsim `origin/v4` @ `d83f4f87` (4.0.0rc4) and tbsim `origin/v4` @ `619c296` (`versions.json`). Exps 12–25 ran at starsim `80c81f53` and tbsim `dcfb522`.

**Setting:** [exp 23](../23_treatment/SUMMARY.md)'s treatment checks, rerun at the new pins.

- **Treatment:** `tbode4.TBTreat` (CONVENTIONS.md 9), a `tbsim.TB` subclass that runs as v4's ODE or as agents.
- **Parameters:** exp 21's set 1 natural history, with tb_LAI_TPT's best-fit treatment values, and homogeneous mixing.
- **Agents:** 200,000 per run, 40 reps per arm. The arms are a fixed 6-month course or a rate, at monthly or weekly steps.
- **Scenarios:** 2024's treatment held constant (equilibrium, compared over 2005–2049), and tb_LAI_TPT's scale-up from 1995 (compared over 1995–2035).

## What moving the pins needed

- **All of `v4/tests/` passes:**
  - v4's ODE still reproduces tbsim's `TB_ODE` and our `TBODE`;
  - `fastode` still reproduces v4's ODE;
  - the latent-structure calculations pass;
  - so does `TBTreat`'s ODE against `fastode`.

  tbsim's default parameters are unchanged. Its two new ones, `immigration` and `emigration`, are 0 by default.
- **`TBTreat` needed four changes to how it declares the model, none to what it does,** all from Starsim's API:
  - `ss.State` takes the label first, so the state's name is passed as `name=`;
  - a transition with several destinations takes a dict of each destination and its probability, not a list with `p=`;
  - births that replace TB deaths come from `ss.unborn`, not `ss.null`;
  - `ss.timeseries` takes `interpolation=`, not `mode=`.
- **The starter script** (`starter/run_tbsim.py`) runs unchanged.

## Finding 1: Starsim now indexes a step's flows by the step's start; the dynamics are unchanged

**What changed:** at the new pins, a flow result at index i is the flow from t_i to t_(i+1). States at index i are still at t_i. At the old pins, the flow at index i was the step ending at t_i.

- **Same for both solvers:** this holds for the ODE and for agents. In a one-way test model, each flow equals its source state's change from i to i + 1 exactly.
- **Our own analyzer agrees:** `tbode4.Incidence` puts each agent's transitions at the same index as Starsim's flows.

**What this broke:** code that dated flows by the step's end then put each flow one step early. Two places did so:

- `tbode4.ode_incidence`, which rebuilds the ODE's first-entry incidence;
- exp 23's `yearly()`.

**The effect:** in the scale-up, the ODE's yearly incidence and TB mortality moved by up to 0.5% at monthly steps (0.1% weekly), and treatment appeared to start in 1994. Prevalence was identical.

**The fix:** both now follow the new indexing. With it, v4's ODE reproduces exp 23's yearly outputs to within 10⁻⁷ per 100,000 in every arm, over time and at equilibrium, incidence included.

Agents ÷ ODE barely changed, because agents and the ODE had shifted together. CONVENTIONS.md 7 now says to check this after any bump.

## Finding 2: agents still match the ODE with treatment

| Agents ÷ ODE (40 reps per arm, mean ± SE) | Prevalence | Incidence per prevalent case | Notifications per prevalent case | On treatment per prevalent case | TB mortality per prevalent case |
|---|---|---|---|---|---|
| Scale-up, pooled over the 4 arms | 1.003 ± 0.005 | 0.998–1.003 | 0.997–1.002 | 0.996–1.005 | 0.993–1.000 |
| Equilibrium, pooled over the 4 arms | 0.992 ± 0.008 | 0.998–1.005 | 0.994–1.002 | 0.993–1.007 | 0.994–1.003 |
| Exp 23 (old pins): scale-up / equilibrium | 1.006 ± 0.004 / 0.983 ± 0.008 | 0.998–1.005 | 0.995–1.007 | 0.995–1.011 | 0.993–1.010 |

The ranges are over the 4 arms. Per arm, the SEs are about 0.002 for incidence and 0.004–0.005 for the other ratios per prevalent case (`agents_vs_ode.csv`, `per_prevalent_case.csv`).

**Prevalence per arm** is 0.997–1.014, except at equilibrium with a fixed course at weekly steps: 0.952 ± 0.016. That arm was also low in exp 23 (0.962 ± 0.017).

- **Why more reps were run:** exp 26 reran exp 23's seeds, so 60 extra reps with new seeds were run for that arm.
- **The extra reps give 1.007 ± 0.013,** and all 100 reps 0.985 ± 0.011 (`extra_reps.csv`), so the low reading is noise.
- **Why the spread is large:** the treated equilibrium has only about 160 people with TB at 200,000 agents, so prevalence varies widely between reps (10–90%: 0.85–1.11).

**TB deaths on treatment read slightly low in agents.** Per prevalent case they are 0.968–0.999 by arm (exp 23: 0.996–1.011). The reps of the two experiments are uncorrelated (r = −0.01), so the two results are independent.

| TB deaths on treatment per prevalent case, agents ÷ ODE (4 arms × 40 reps, mean ± SE) | Exp 26 | Exp 23 |
|---|---|---|
| Scale-up | 0.983 ± 0.006 | 0.999 ± 0.006 |
| Equilibrium | 0.994 ± 0.006 | 1.005 ± 0.006 |

Two checks isolate treatment, each with 8 × 400,000 agents at monthly steps:

- **A cohort entering failing courses at a constant rate** (`fail_death.py`): agents give 0.998–0.999 (SE 0.001) of the exact share of courses that end in TB death, at both steps and with either course. The ODE gives 0.998–1.001.
- **Treatment varying in time** (`timeseries_check.py`). People start treatment from one state at tb_LAI_TPT's time-varying rate, with a time-varying share set to cure. Failing courses face TB death, and every course returns to that state.
  - **Arms:** the scale-up's `ss.timeseries` against 2024's values held constant, each with a fixed course or a rate.
  - **Deaths per person on a failing course:** agents give 0.996–0.997 (SE 0.001) of the ODE in all four arms (`timeseries_check.csv`).
  - **Starts and person-time on failing courses:** within 0.2%.

**What this shows:**

- **Time-varying treatment isn't the cause:** the scale-up and constant arms agree.
- **Agents show a small real shortfall:** about 0.4% fewer TB deaths on treatment per person-time at monthly steps, about 4 SE below 1. It isn't explained, and it is too small to matter.
- **The rest of the full model's scale-up shortfall is within noise:** 0.983 ± 0.006 against about 0.996 expected is 2 SE.
- **The effect on mortality is small:** deaths on treatment are about 19% of TB deaths at the treated equilibrium, so TB mortality moves by about 0.1%. That is well within the per-case agreement above.

## Finding 3: treatment has none of the pitfalls, and the ODE's fixed course is unchanged

Every pitfall check passes as in exp 23 (`pitfall_checks.csv`; 96,000–106,000 courses per arm):

- **Starts and relapses:** nobody starts treatment unless symptomatic, and nobody relapses without a cure.
- **Agents flagged:** 12 agents over all arms are flagged as starting while on treatment, or as dying of TB on a course that will cure.
  - Each was cured, relapsed, and then restarted or died, all within one step (`diagnose_cases.csv`), which continuous time allows.
  - Exp 23 also flagged 12.
- **Restarts in the step of failure** match their expectation: 4.53% against 4.60% at monthly steps, and 1.08% against 1.08% at weekly steps.
- **Course length:** a fixed course lasts exactly 6 steps at monthly steps. At weekly steps it lasts 26 or 27 weeks (93% / 7%; mean 0.5000 years).
- **Other checks:**
  - the share of courses set to cure tracks efficacy by year of starting (`success_by_year.csv`);
  - background death on treatment is 0.013–0.016 per year, against μ = 0.0143;
  - the population stays within ±1%.

**Fixed course against a rate:** in the ODE, the fixed course still gives about 1% less prevalence at equilibrium (78.76 against 79.61 per 100,000), and 1.25% less by 2024 in the scale-up (`fixed_vs_rate.csv`).

## What this means

- **The v4 track can move to the new pins.** Treatment, the anchor and `fastode` all carry over. Homogeneous agents were already checked at these pins (exp 25, finding 5).
- **Results of exps 12–25 stand for homogeneous mixing.** Rerunning those experiments themselves needs the old pins.
- **On a network, the new pins change the convention:** calibrate on Starsim's ODE of the network, which is the same model as agents on `RandomNet(dur=0)` (exp 25; CONVENTIONS.md 1).
- **Run cost is unchanged:** 18 minutes for the 320 agent runs on 4 cores, as in exp 23.

## Files

- **Exp 23's `run.py` and `diagnose.py`,** run with `TB_EXP23_FOLDER` pointing here. Outputs:
  - `runs.obj`, `runs_by_group/`, `ode.obj`
  - `agents_vs_ode.csv`, `per_prevalent_case.csv`, `pitfall_checks.csv`, `diagnose_cases.csv`, `success_by_year.csv`
  - `fixed_vs_rate.csv`, `ode_equilibria.csv`, `ode_fail_death_share.csv`
  - figures, `run.log`, `diagnose.log`, `versions.json`
- **`extra_reps.py`:** 60 more reps of the equilibrium, fixed-course, weekly arm → `extra_reps.obj`, `extra_reps.csv`, `extra_reps.log`.
- **`fail_death.py`:** the share of failing courses that end in TB death → `fail_death_runs.csv`, `fail_death.csv`, `fail_death.log`.
- **`timeseries_check.py`:** TB deaths on treatment with time-varying against constant treatment → `timeseries_check_runs.csv`, `timeseries_check.csv`, `timeseries_check.log`.
