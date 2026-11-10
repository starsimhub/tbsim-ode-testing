# 24: Summary

**Track:** v4 (tbsim `dcfb522`, starsim `80c81f53`; `versions.json`).

**Setting:**

- **The model:** `tbode4.fastode` with exp 23's treatment, as a fixed 6-month course (12 Erlang stages), with `sym_dead` falling linearly from 1995 to 2024, and with reported notifications = (1 + `over_detection`) × treatment starts.
- **Each draw:** starts at its untreated equilibrium at the end of 1994, then follows tb_LAI_TPT's coverage and efficacy trajectories to 2040. About 90 ms per draw.
- **Calibration:** 21 parameters, 10 targets, IMIS. Each IMIS run is 310,000 draws, stopped at its 300-iteration cap, about 1.6 hours on 4 cores.
- **Not in the model:** ages, demography and a household network. All-ages prevalence is compared with the 15+ survey.

## Coverage check: every target is within reach

Every target's interval lies inside the 5–95% range of 5,000 prior draws' outputs (`coverage.csv`). The narrowest is prevalence, which 3.7% of draws reach. No draw meets all 10 targets at once, against 1 of 5,000 for exp 20's 5. 46% of draws have no endemic equilibrium.

## Re-identification: the method finds the best fit; six parameters aren't identified

**The truth:** exp 21's set 1 natural history, with tb_LAI_TPT's best-fit treatment, `sym_dead` trend and `over_detection`.

**IMIS:** stopped at the cap, with an effective sample size of 1,027 and 1,646 expected unique points out of the 1,896 the stopping rule asks for.

- **Every synthetic target is fitted** (`reid_fit.csv`).
- **IMIS found the best-fitting region:** the log-likelihood at the truth, 4.26, is no higher than the best point found, 4.19.
- **15 of 21 parameters cover the truth in their 95% interval.**
- **The 6 that don't:** `inf_non`, φ, `non_rec`, `non_asy`, `p_relapse` and `over_detection`. Points far from the truth in these directions fit as well as the truth, so the targets don't identify them. Several truths sit near the edge of their prior (`non_asy` at 0.97 of its range, `non_rec` 0.90, `over_detection` 0.92), where the posterior's mass falls short of them.

## The real targets: all fitted, with three in some tension

Posterior predictive median (95% interval), first seed, against the target:

| Target | Posterior | Target (95% interval) |
|---|---|---|
| Prevalence 2020 (per 100k) | 306 (232–400) | 310 (247–385) |
| Share asymptomatic | 0.654 (0.611–0.696) | 0.654 (0.584–0.719) |
| Prevalence ÷ incidence (y) | 1.46 (1.26–1.74) | 1.57 (1.20–2.05) |
| **Mortality ÷ incidence** | **0.116 (0.091–0.148)** | 0.113 (0.088–0.144) |
| Prevalence ÷ notifications (y) | 3.08 (2.58–3.85) | 2.66 (2.08–3.40) |
| Incidence 2024 ÷ 2015 | 0.80 (0.69–0.89) | 0.79 (0.63–0.98) |
| Mortality 2024 ÷ 2015 | 0.58 (0.50–0.65) | 0.66 (0.52–0.83) |
| Prevalence 2040 ÷ 2020 | 0.65 (0.47–0.81) | 0.81 (0.56–1.17) |
| ProgressionY1 | 0.053 (0.045–0.065) | 0.052 (0.044–0.063) |
| ProgressionRatio | 0.645 (0.579–0.714) | 0.663 (0.577–0.762) |

- **Mortality ÷ incidence is now fitted.** It was out of reach without treatment (exp 20). Treatment brings deaths per case down to the WHO value.
- **Three targets are in some tension,** though their posteriors overlap the target intervals:
  - **Notifications are too few:** about 99 per 100,000 in 2020 against 117. The model's notifications stay roughly flat over 2015–2024, while the data rise to 173 in 2024. Only 2020 is a target.
  - **Mortality falls too fast** over 2015–2024: 0.58 against 0.66.
  - **Prevalence falls too fast** over 2020–2040: 0.65 against tb_LAI_TPT's judgement-based anchor of 0.81.
  - **Possible causes,** not tested here: the missing demography (India's population grew and aged), and the trajectories' fixed shapes (linear coverage and `sym_dead` between knots).
- **The burden before treatment:** prevalence about 530 per 100,000 in 1994 (330–950), falling to 306 by 2020 and about 200 by 2039 (`real_trajectories.png`). Exp 20 fitted the 2020 prevalence without treatment, so its untreated equilibrium sat at the 2020 level.
- **What's identified:**
  - `beta`: 12–32 per year.
  - π: 0.22–0.53.
  - Coverage in 2024: 0.55–0.69.
  - `over_detection`: 0.07–0.16.
  - `sym_dead` in 2024: about half its 1995 value.
  - φ: 0.02–0.09, now informed by the trend targets. Exp 20 couldn't identify it.

## Convergence: predictions agree between seeds, parameters along a ridge don't

The second seed also stopped at the cap, with an effective sample size of 961 and 1,574 expected unique points.

- **Predictions agree:**
  - Every target's posterior predictive is within a few percent between seeds, e.g. prevalence ÷ notifications 3.08 against 2.88, prevalence 2040 ÷ 2020 0.652 against 0.646.
  - The trajectories agree too, e.g. prevalence in 2024: 202–357 against 194–348 per 100,000; TB mortality in 2039: 9–23 against 10–25.
- **Several parameters don't agree:**

  | Parameter, median (95%) | Seed 1 | Seed 2 |
  |---|---|---|
  | `k_latent` | 8.3 (5.0–11.1) | 4.4 (1.7–8.0) |
  | `inf_asy` | 0.36 (0.24–0.49) | 0.22 (0.13–0.37) |
  | `inf_non` | 0.33 (0.23–0.43) | 0.21 (0.09–0.38) |
  | `p_relapse` | 0.049 (0.027–0.084) | 0.077 (0.045–0.096) |
  | `efficacy_frac` | 0.26 (0.02–0.64) | 0.15 (0.01–0.43) |

- **The ridge:** `k_latent` and `inf_asy` are strongly correlated in both seeds (0.76, 0.86). Faster front-loading with a higher initial progression rate meets the progression targets as well as slower front-loading with a lower rate. Each seed's IMIS run weighted a different stretch of this ridge, and neither reached its stopping rule.
- **So:**
  - The posterior's predictions of burden over time are reliable.
  - Its parameter marginals along the ridge are not. That includes the shape of early latent progression, which matters for any analysis of preventive therapy (tb_LAI_TPT's purpose).

## What this means

- **Treatment, time-varying TB mortality and over-detection let the ODE meet every tb_LAI_TPT target** that applies without ages, demography or a household network. Mortality ÷ incidence included.
- **For checking agents against the ODE, the posterior is good enough:** parameter sets drawn from it fit the targets, whichever seed they come from.
- **For parameter inference, it isn't converged yet.** The options:
  - run IMIS much longer;
  - reparameterise along the ridge, e.g. by the progression targets themselves;
  - or use a sampler that moves along ridges, such as MCMC or SMC.

## Files

- `calib.py` (the model, priors, targets, likelihood, IMIS); `run.py`
- `coverage.obj`, `coverage.csv`, `coverage.log`
- `reid.obj`, `reid_posterior.csv`, `reid_fit.csv`, `reid_marginals.png`, `reid_fit.png`, `reid_trajectories.png`, `reid_trajectories.csv`, `reid.log`
- `real.obj`, `real_posterior.csv`, `real_fit.csv`, `real_marginals.png`, `real_fit.png`, `real_trajectories.png`, `real_trajectories.csv`, `real.log`
- `real_seed2.obj`, `real_seed2_posterior.csv`, `real_seed2_fit.csv`, `real_seed2_trajectories.csv`, `real_seed2.log`
- `versions.json`
