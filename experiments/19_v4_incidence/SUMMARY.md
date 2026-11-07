# 19: Summary

**Track:** v4 (tbsim `origin/v4` @ `dcfb522`, starsim `origin/v4` @ `80c81f53`; `versions.json`).

**Setting:**

- **Definition:** incidence is the first entry to asymptomatic TB per episode (tb_LAI_TPT's definition, CONVENTIONS.md 5).
- **ODE:** a standalone copy of v4's ODE, with non-infectious TB split into "not yet asymptomatic this episode" (NON0) and "has been" (NON1), with identical rates.
- **Agents:** an analyzer counts the same from transitions, with a per-agent flag for "has been asymptomatic this episode".
- **Parameter sets:** set A (tbsim's defaults) and set B (exp 17's two-state model at φ = 0.01, latent reinfection on), both calibrated to 500/100k. The calibrated `beta`s, 9.9265 and 21.1434, match exps 13 and 17.
- **Runs:** 50,000 agents, homogeneous mixing, started at the ODE's equilibrium, 1900–2000, 40 reps per arm.

## Question 3: the standalone ODE reproduces v4's ODE

The copy writes tbsim's declared equations as one function and finds the equilibrium by root-finding. `v4/tests/test_fastode.py` checks it:

- **At equilibrium:** within 1e-14 for every state in both sets. v4's ODE, run from the copy's equilibrium, stays within 1e-9 over 50 years.
- **Over time:** over 100 years from a seed, within 2e-8 at monthly steps and 2e-11 at daily steps. The difference shrinks with v4's step, so it's v4's fixed-step integration, not the model.
- **Speed:** an equilibrium takes 6–14 ms, against 0.1–0.3 s through Starsim.

## Question 1: agents match the ODE

| Incidence per prevalent case, agents ÷ ODE, 1950–2000 (40 reps, mean ± SE) | Set A, monthly | Set A, weekly | Set B, monthly | Set B, weekly |
|---|---|---|---|---|
| First entry per episode | 1.000 ± 0.002 | 1.002 ± 0.002 | 0.996 ± 0.002 | 1.001 ± 0.002 |
| Every entry | 0.998 ± 0.002 | 1.002 ± 0.002 | 0.996 ± 0.002 | 1.000 ± 0.002 |
| Share of entries that are first entries: agents (ODE) | 0.539 ± 0.001 (0.538) | 0.538 ± 0.001 (0.538) | 0.500 ± 0.001 (0.500) | 0.500 ± 0.001 (0.500) |

Per prevalent case, because each rep's prevalence varies by about ±10% (10–90% range). Incidence itself, agents ÷ ODE, follows each arm's prevalence (0.97–1.01; `summary.csv`).

## Question 2: about half of all entries to asymptomatic TB are re-entries

| At 500/100k (ODE) | Set A | Set B |
|---|---|---|
| Incidence, first entry per episode | 274/100k per year | 275/100k per year |
| Every entry to asymptomatic TB | 510 | 551 |
| Re-entries (share of every entry) | 46% | 50% |
| Prevalence ÷ incidence, first entry | 1.82 years | 1.82 years |
| Prevalence ÷ incidence, every entry | 0.98 years | 0.91 years |

- **Why so many:** people move between asymptomatic and non-infectious TB quickly. `asy_non` is 1.7 per year in both sets, so an episode often passes through asymptomatic TB more than once.
- **For exp 20's targets:** tb_LAI_TPT's prevalence ÷ incidence target is 1.57 years (1.20–2.05). The first-entry count gives 1.82 at 500/100k, inside that interval. Every entry would give about 0.9, far outside it. So the definition matters for the calibration, and the first-entry count is the one that matches tb_LAI_TPT.
- **v4's existing output misses remote progression:** `v4/common.py`'s `outputs()` counts latent → asymptomatic and non-infectious → asymptomatic, but not remote latent → asymptomatic. With the two-state model that misses 2% of entries in set B (540 against 551/100k). Its numbers in exps 13–17 were for incidence per prevalent case in set A and in set B without remote progression, where nothing was missed.

## What this means

- **Exp 20 can use the standalone ODE** (`tbode4.fastode`): it's v4's model, about 30 times faster at equilibrium, and it gives incidence by tb_LAI_TPT's definition.
- **Agents can be compared on the same definition** with this experiment's analyzer (`Incidence` in `run.py`).
- **v4's every-entry count isn't incidence by our definition:** it's about twice as large. Comparisons and targets should use the first-entry count.

## Files

- `v4/tbode4/fastode.py`, `v4/tests/test_fastode.py`
- `run.py` → `results.obj`, `summary.csv`, `run.log`, `versions.json`
