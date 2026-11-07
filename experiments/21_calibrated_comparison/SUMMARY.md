# 21: Summary

**Track:** v4 (tbsim `origin/v4` @ `dcfb522`, starsim `origin/v4` @ `80c81f53`; `versions.json`).

**Setting:**

- **Parameter sets:** 12, drawn from exp 20's calibrated posterior by weighted resampling without replacement (seed 21; `sets.obj`, `sets.log`). They span the posterior: `beta` 10.9–29.0, `k_latent` 3.6–11.7, φ 0.012–0.093, π 0.20–0.66, each close to its posterior 95% interval.
- **Model:** v4's `tbsim.TB` with the two-state latent model, homogeneous mixing and no treatment.
- **Runs per set:** v4's ODE, and 50,000 agents × 20 reps at monthly and weekly `dt`, both started at the standalone ODE's equilibrium, 1900–2000.
- **Incidence:** the first entry to asymptomatic TB per episode (exp 19's analyzer for agents; the standalone ODE for the reference).
- **The check on the mapping:** v4's ODE stays at the standalone equilibrium to 4e-16, so the parameters map correctly.

## Finding: agents match the ODE across the posterior

| Agents ÷ ODE, 1950–2000, averaged over 12 sets (mean ± SE between sets; range of the 12 set means) | Monthly | Weekly |
|---|---|---|
| Active prevalence | 0.997 ± 0.007 (0.96–1.03) | 0.987 ± 0.007 (0.94–1.04) |
| Share asymptomatic | 0.998 ± 0.001 | 1.000 ± 0.001 |
| Incidence (first entry) per prevalent case | 0.998 ± 0.001 | 1.001 ± 0.001 |
| Prevalence ÷ incidence | 1.002 ± 0.001 | 0.999 ± 0.001 |
| TB mortality | 0.993 ± 0.007 | 0.986 ± 0.006 |
| Mortality ÷ incidence (not a fitted target) | 0.999 ± 0.003 | 0.998 ± 0.002 |
| Remote latent | 0.995 ± 0.006 | 0.989 ± 0.006 |
| Cleared (all three no-longer-infected states) | 0.996 ± 0.002 | 0.996 ± 0.001 |
| Recent latent (end-of-step snapshot; see below) | 1.42 ± 0.04 (1.20–1.67) | 1.08 ± 0.01 |

`pooled.csv`; by set, `by_set.csv` and `agents_vs_ode.png`.

- **By set:** every set's prevalence, mortality and ratios are within about 2 SE of 1. There are 12 sets × 2 steps × 8 outputs, so a few at about 2.5 SE are expected. The largest |z| for the ratios is 2.8; for prevalence it's 2.1.
- **Ratios per prevalent case match to 0.3%:** the share asymptomatic, incidence, prevalence ÷ incidence and mortality ÷ incidence. That holds across every part of the posterior these sets reach, including fast front-loading (`k_latent` up to 11.7), large φ (0.093) and high π (0.66).
- **Recent latent infection reads high in monthly snapshots, as predicted.** v4 dates infections at the end of their step, so a snapshot counts the last step's infections in full. The excess is x/(1 − e^(−x)), with x the rate of leaving the recent state × `dt`. That predicts each set's value within 2 SE (observed ÷ predicted 0.96–1.02). It's larger than in exp 17 (×1.18–1.33) because these sets include faster front-loading, so recent infections are shorter-lived. Total latent infection reads 1.04 monthly for the same reason. The mean over a step is right.
- **The cleared share is 0.4% low,** the small finite-population effect found in exp 17.

**The targets:** at every set, the agents' prevalence, share asymptomatic and prevalence ÷ incidence sit on the ODE's (`targets.png`). So agents fit the targets as the ODE does. A few sets lie outside a target's interval, e.g. set 5's prevalence (222/100k) and set 11's (404/100k). That's expected for posterior draws: the posterior predictive interval for prevalence is 233–391.

**Weekly prevalence here is 0.987 ± 0.007 (1.9 SE below 1), most likely noise.**

- **Across every v4 arm with homogeneous mixing,** pooled by inverse variance: weekly prevalence is 0.993 ± 0.003 (8 arms, exps 13, 17, 19 and 21), against monthly 1.001 ± 0.003 (9 arms).
- **Excluded from that pooling:** the reps exp 13 added *because* its first 40 looked low; only the other 160 are used.
- **Why it's probably chance:**
  - The difference is 2.5 SE at most.
  - It comes from many comparisons.
  - There's no plausible mechanism: the known step effects are larger at monthly steps, not weekly.
- **The test, if weekly runs ever matter:** one set, 200 reps at each step.

## What this means

- **v4's agents and ODE are consistent at the calibrated parameters,** across the posterior's range, with incidence counted as tb_LAI_TPT counts it. So the exp 20 calibration carries over from the ODE to the agents, for the targets, mortality and the latent shares. The exception is the known end-of-step snapshot of recent infection.
- **The 12 sets** (`sets.obj`) are a ready sample of the posterior, e.g. for the treatment work.
