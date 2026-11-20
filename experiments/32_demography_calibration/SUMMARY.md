# 32: Summary

**Track:** v4 (tbsim `619c296`, starsim `d83f4f87`). Run on zebra (160 workers) on 2026-11-20. Zebra's checkout is partial, so this experiment's files were copied there from the local working tree at commit `4ce4473` plus this experiment's uncommitted files; `versions.json` records zebra's own older commit (`e486c48`) as dirty.

**Setting:**

- **The model:** `tbode4.fastode_age.run_cohort()`, exp 24's TB model and treatment with India's WPP demography and exact aging by cohort, at monthly steps with the force of infection at mid-step (CONVENTIONS.md 8).
- **Each draw:** the untreated equilibrium without ages, a 100-year frozen start at 1950's demography, then 1950–2041.
- **Calibration:** exp 24's 21 parameters and priors and its 10 targets, with every target built on the prevalence survey using 15+ prevalence (prevalence, share asymptomatic, prevalence ÷ incidence, prevalence ÷ notifications, prevalence 2040 ÷ 2020), as tb_LAI_TPT computes them.
- **IMIS:** exp 24's settings. Each run was 310,000 draws and stopped at its 300-iteration cap, short of the stopping rule (1,621–1,655 expected unique points against 1,896), as exp 24's did. A run took 1.3–2.2 hours on 160 workers, against the 45 minutes the design estimated: a draw took about 2 s per worker on zebra, against 1.35 s on the laptop, and posterior draws were slower than the prior's average.
- **Not in the model:** child transmission (tb_LAI_TPT sets it to 0 under age 10), children's lower treatment coverage, and a household network.

## Coverage check: every target is within reach

Every target's 95% interval lies inside the 5–95% range of 5,000 prior draws' outputs (`coverage.csv`). The shares inside each interval are close to exp 24's, e.g. 4.2% for 15+ prevalence (3.7% for exp 24's all-ages prevalence). As in exp 24, no draw meets all 10 targets at once, and 46% of draws have no endemic state.

## Re-identification: as exp 24, the method finds the best fit and six parameters aren't identified

**The truth:** exp 24's (exp 21's set 1 natural history with tb_LAI_TPT's best-fit treatment, `sym_dead` trend and over-detection). With demography its 15+ prevalence in 2020 is 220 per 100,000 (the synthetic target).

- **Every synthetic target is fitted** (`reid_fit.csv`).
- **IMIS found the best-fitting region:** the log-likelihood at the truth is 4.05 and the best point found 3.91.
- **15 of 21 parameters cover the truth** in their 95% interval.
- **The 6 that don't:** `inf_non`, `inf_asy`, `k_latent`, φ, `non_asy` and `over_detection`. For each, the best point at least as far from the truth as the posterior median fits within 0.14–0.19 of the truth's log-likelihood (`reid_posterior.csv`, `max_ll_far_minus_truth`), so the targets don't identify them. `k_latent` and `inf_asy` lie along the ridge exp 24 found (below). `non_asy` (0.97 of its prior range) and `over_detection` (0.92) sit near the edge of their prior, as in exp 24.

## The real targets: all fitted, with the same three tensions as exp 24

Posterior predictive median (95% interval) for each seed, against the target and exp 24's first seed:

| Target | Seed 1 | Seed 2 | Exp 24, seed 1 (all-ages prevalence) | Target (95% interval) |
|---|---|---|---|---|
| Prevalence 15+, 2020 (per 100k) | 308 (231–394) | 306 (232–403) | 306 (232–400) | 310 (247–385) |
| Share asymptomatic, 15+ | 0.650 (0.603–0.697) | 0.641 (0.590–0.697) | 0.654 (0.611–0.696) | 0.654 (0.584–0.719) |
| Prevalence 15+ ÷ incidence (y) | 1.38 (1.13–1.73) | 1.48 (1.26–1.79) | 1.46 (1.26–1.74) | 1.57 (1.20–2.05) |
| Mortality ÷ incidence | 0.119 (0.094–0.151) | 0.116 (0.092–0.146) | 0.116 (0.091–0.148) | 0.113 (0.088–0.144) |
| Prevalence 15+ ÷ notifications (y) | 2.90 (2.40–3.50) | 2.80 (2.30–3.52) | 3.08 (2.58–3.86) | 2.66 (2.08–3.40) |
| Incidence 2024 ÷ 2015 | 0.77 (0.69–0.83) | 0.80 (0.71–0.90) | 0.80 (0.69–0.89) | 0.79 (0.63–0.98) |
| Mortality 2024 ÷ 2015 | 0.58 (0.49–0.65) | 0.60 (0.52–0.67) | 0.58 (0.50–0.65) | 0.66 (0.52–0.83) |
| Prevalence 15+, 2040 ÷ 2020 | 0.65 (0.50–0.79) | 0.67 (0.48–0.85) | 0.65 (0.47–0.81) | 0.81 (0.56–1.17) |
| ProgressionY1 | 0.052 (0.044–0.063) | 0.053 (0.044–0.065) | 0.053 (0.045–0.065) | 0.052 (0.044–0.063) |
| ProgressionRatio | 0.651 (0.589–0.729) | 0.656 (0.579–0.736) | 0.645 (0.579–0.714) | 0.663 (0.577–0.762) |

- **Demography and the 15+ targets leave the fit much as it was.** With ages and demography the ODE meets every target, and the posterior predictives are close to exp 24's.
- **The three tensions remain,** so the missing demography wasn't their cause:
  - **Notifications are too few:** 106 per 100,000 in 2020 (seed 1; 74–143), against 117. Somewhat closer than exp 24's 99. The model's notifications peak around 2009 and fall after, while the data rise to 173 in 2024 (`real_trajectories.png`).
  - **Mortality falls too fast** over 2015–2024: 0.58–0.60 against 0.66.
  - **Prevalence falls too fast** over 2020–2040: 0.65–0.67 against tb_LAI_TPT's judgement-based anchor of 0.81.
  - The other possible cause exp 24 named, the trajectories' fixed shapes (linear coverage and `sym_dead` between knots), is untested.
- **Burden over time** (`real_trajectories.png`, 200 posterior draws): 15+ prevalence about 730 per 100,000 from 1950 to 1995 (seed 1; 450–1,150), falling to 308 by 2020 and about 200 by 2040. All-ages prevalence in 2020 is 352 (249–446), above the 15+ value, because children transmit and progress like adults here. The population grows 3.7-fold by 2020, as WPP's does.

## Convergence: predictions in the targeted years agree between seeds; parameters agree less well than in exp 24

**Predictions:** every target's posterior predictive agrees between seeds within a few percent (table above). Outside the targeted years they agree less well: 15+ prevalence in 1995 has a median of 727 per 100,000 in seed 1 and 550 in seed 2, with wide, overlapping intervals (435–1,160 and 346–838). The best points found also differ: log-likelihood 3.10 in seed 1, 3.41 in seed 2.

**Parameters:** a parameter counts as agreeing if each seed's median lies within the other's central 50% interval (`seeds_parameters.csv`, `compare_marginals.png`).

- **5 of 21 agree:** `non_rec`, `sym_dead_1995`, `tx_covg_2010_frac`, `p_relapse`, `over_detection`. By the same criterion, exp 24's two seeds agreed on 8.
- **The largest disagreements,** median (95% interval):

  | Parameter | Seed 1 | Seed 2 | Exp 24, seed 1 | Exp 24, seed 2 |
  |---|---|---|---|---|
  | `k_latent` | 7.8 (4.7–10.9) | 5.7 (2.7–8.8) | 8.3 (5.0–11.1) | 4.4 (1.8–8.0) |
  | `inf_asy` | 0.33 (0.22–0.47) | 0.26 (0.16–0.38) | 0.36 (0.24–0.49) | 0.23 (0.13–0.37) |
  | `inf_non` | 0.30 (0.18–0.40) | 0.18 (0.08–0.35) | 0.33 (0.23–0.44) | 0.22 (0.09–0.39) |
  | `trans_asymp` | 0.85 (0.65–0.98) | 0.70 (0.55–0.85) | 0.86 (0.69–0.99) | 0.78 (0.64–0.93) |
  | `efficacy_1995` | 0.64 (0.40–0.86) | 0.42 (0.31–0.58) | 0.44 (0.31–0.76) | 0.42 (0.32–0.56) |
  | π | 0.25 (0.21–0.37) | 0.34 (0.22–0.58) | 0.33 (0.22–0.53) | 0.29 (0.21–0.48) |

- **The ridge remains:** the weighted correlation of `k_latent` with log `inf_asy` is 0.77 and 0.86 in the two seeds, as in exp 24 (0.77, 0.85). The progression targets are unchanged, so demography couldn't remove it.
- **New here:** the seeds also disagree on `efficacy_1995` and on `trans_asymp` and π, which exp 24's seeds agreed on or nearly did. Seed 1 puts more weight on a high initial treatment efficacy, which goes with a higher burden before 1995 (weighted correlation of `efficacy_1995` with log 15+ prevalence in 1995: 0.39 and 0.49 in the two seeds): a high early burden brought down faster by treatment fits the 2020 targets as well as a lower one brought down more slowly.
- **So:**
  - As in exp 24, neither IMIS run reached its stopping rule, and each weighted different parts of a broad, ridged posterior.
  - The posterior's predictions of the targeted quantities are reliable. Its parameter marginals, and its burden before 1995 (which no target constrains), are not.

## What this means

- **The recalibration with demography works as a fit:** with ages, India's demography and 15+ targets, the ODE meets every tb_LAI_TPT target that applies without child transmission or a household network, about as well as exp 24 did without demography.
- **For checking agents against this ODE** (the next step on this row), the posterior is good enough: parameter sets drawn from either seed fit the targets (`real_draws.csv`, `real_seed2_draws.csv`, 200 weighted draws each).
- **For parameter inference it isn't converged,** less so than exp 24. The options are exp 24's: run IMIS much longer, reparameterise along the ridge (e.g. by the progression targets themselves), or use a sampler that moves along ridges (MCMC or SMC). On zebra, a run several times longer is affordable (each run here took about 2 hours).
- **Before relying on this posterior,** exp 33 checks whether Starsim's new default for births' female share (the fraction at ages 0–4, not the starting population's 48.63%) moves these targets by more than about 1%.

## Files

- **Code:** `calib.py` (the model per draw, targets, IMIS with a persistent pool), `run.py` (`coverage`, `reid`, `real`, `real_seed2`, `summary`).
- **Raw runs (not committed):** `coverage.obj`, `reid.obj`, `real.obj`, `real_seed2.obj`.
- **Tables:** `coverage.csv`; `{reid,real,real_seed2}_posterior.csv`, `_fit.csv`, `_draws.csv` (200 weighted posterior draws), `_trajectories.csv` (those draws by calendar year, 1950–2040); `seeds_parameters.csv`, `seeds_targets.csv`.
- **Figures:** `{reid,real,real_seed2}_marginals.png`, `_fit.png`, `_trajectories.png`; `compare_marginals.png` (both seeds against exp 24's).
- **Logs:** `runs.out`, `coverage.log`, `{reid,real,real_seed2}_summary.log`, `compare.log`, `versions.json`.
