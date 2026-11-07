# 20: Calibrating the v4 ODE to tb_LAI_TPT's targets, without treatment

**Status:** design agreed with the user (2026-11-06); run (see SUMMARY.md). Built on exp 19 (incidence in v4).
**Track:** v4 (`v4/.venv`; tbsim `origin/v4` @ `dcfb522`, starsim `origin/v4` @ `80c81f53`). For homogeneous mixing, exp 18's newer commits leave the model and solver unchanged.

## Question

Calibrate the v4 TB ODE to the targets of tb_LAI_TPT's calibration that apply to a model without treatment, a household network or time trends. The model has the two-state latent structure (CONVENTIONS.md 2) and homogeneous mixing (CONVENTIONS.md 1).

**Purpose:** a parameter set, with its uncertainty, at which ODE and agents can be compared again, with reasonable burden levels. It doesn't describe untreated TB. The model has no treatment, so its natural history will absorb treatment's effect on 2020–24 burden (decided by the user, 2026-11-06). The model will be recalibrated with treatment before the treatment experiments.

## Targets and likelihoods

tb_LAI_TPT's calibration (`run_history_matching.py`, `run_trajectory_selection.py`, `hm_common.make_components_hm`, `cohort_progression.py`) has 11 targets. History matching scores each as Normal on the log scale, with the data σ. Trajectory selection uses BetaBinomial likelihoods for the two survey proportions and for four of the ratios (the "BB-on-split" form, moment-matched to the same log σ, which also allows for the agents' finite-cohort noise), and Normal-in-log for the progression targets.

**Six targets apply in principle; five are fitted.** Mortality ÷ incidence is dropped from the likelihood, but still reported (below). Means are tb_LAI_TPT's, pulled from its components (`utils.make_components_full`, `hm_common.calibration_sigmas`, `cohort_progression.TARGET_MEANS`/`TARGET_SIGMAS`).

| Target | Definition (tb_LAI_TPT) | Source, year | Mean | tb_LAI_TPT likelihood: history matching / trajectory selection | Likelihood for the ODE | Implied 95% interval |
|---|---|---|---|---|---|---|
| Prevalence | Asymptomatic + symptomatic TB, 15+, point prevalence | India national survey, 2020 | 310/100k | log-Normal, σ = 0.114 / BetaBinomial, x = 77, n = 24,815 (effective) | Binomial(x = 77, n = 24,815) | 247–385/100k (Binomial); 248–388 (log-Normal) |
| Prop_Asymptomatic | Asymptomatic ÷ (asymptomatic + symptomatic), 15+ | Same survey, 2020 | 0.654 | log-Normal, σ = 0.053 / BetaBinomial, x = 123, n = 188 | Binomial(x = 123, n = 188) | 0.584–0.720 (Binomial); 0.590–0.726 (log-Normal) |
| Mort_to_Inc (**not fitted**) | TB deaths ÷ incidence, all ages | WHO estimates, 2024 (21.1 ÷ 187.3 per 100k) | 0.113 | log-Normal, σ = 0.127 / BetaBinomial-on-split | Reported only | 0.088–0.144 |
| Prev_to_Inc | Prevalence (15+, survey) ÷ incidence (all ages, WHO), years | Survey 2020 ÷ WHO 2020 (310.3 ÷ 197.5) | 1.57 y | log-Normal, σ = 0.136 / BetaBinomial-on-split | log-Normal, σ = 0.136 | 1.20–2.05 y |
| ProgressionY1 | Share of a no-reinfection latent cohort that has reached asymptomatic TB by 1 year | Ferebee 1970, Sutherland 1968 | 0.0523 | log-Normal, σ = 0.094 / same | log-Normal, σ = 0.094 | 0.0435–0.0629 |
| ProgressionRatio | The same at 1 year ÷ at 10 years | Same | 0.663 | log-Normal, σ = 0.071 / same | log-Normal, σ = 0.071 | 0.577–0.762 |

- **Why Binomial and log-Normal, not tb_LAI_TPT's BetaBinomial forms:** the ODE is deterministic, so there is no finite-cohort simulation noise for the BetaBinomial to carry. With a model proportion p, a BetaBinomial against the survey reduces to Binomial(x; n, p). For the ratios, the data σ is all that's left of the BB-on-split form, i.e. log-Normal, which is also what history matching used.
- **Mortality ÷ incidence is out of reach without treatment,** so it's dropped from this calibration (user, 2026-11-06). It comes back when treatment is added. The coverage check found it at 0.18–0.46 across 90% of 5,000 prior draws (target 0.088–0.144), with a minimum of 0.123 in a corner of the box. Treatment lowers deaths per case; untreated natural history can't, within its priors. Widening the active-disease priors to reach it would make them implausible.
- **Ages:** the ODE has no ages, so its all-ages values are compared with the 15+ prevalence targets. The mismatch is accepted (user, 2026-11-06); all-ages prevalence is lower than 15+ prevalence where children have less TB.
- **Incidence:** the first entry to asymptomatic TB per episode (CONVENTIONS.md 5), as in tb_LAI_TPT, built and tested in exp 19.
- **Equilibrium:** the ODE runs to equilibrium with constant rates, and births replace deaths at `mu` = 1/70. Its equilibrium is compared with the 2020 and 2024 data.
- **The progression targets:** computed from a no-reinfection latent cohort with the exact chain (`tbode4.latent`, as exps 14–17), counted from transitions.

**Five targets don't apply:**

| Target | Why not |
|---|---|
| Prev_to_Notif | Notifications need treatment |
| Inc_Ratio (2024 ÷ 2015), Mort_Ratio (2024 ÷ 2015) | Trends over time; driven by treatment scale-up and time-varying mortality |
| Prev_2040_Ratio | A projection; driven by the same trends |
| Prop_Transmission_HH | No household network |

## Parameters and priors

tb_LAI_TPT's prior box (`run_history_matching.parameter_bounds`), with the treatment, household and trend parameters dropped, and its time-since-infection decay replaced by the two-state parameters. Uniform, or log-uniform where marked.

| Parameter | Prior | Note |
|---|---|---|
| `beta` (homogeneous, per year) | 2–100, log-uniform | tb_LAI_TPT's is per contact on its network, so it doesn't carry over; this range covers exps 13–18's calibrated values (9.9–22) |
| `trans_asymp` | 0.5–1.0 | |
| π (`rr_reinfection_cleared`, `_rec`, `_treat`, `_inf`) | 0.2–1.0 | One parameter for every reinfection risk, as in tb_LAI_TPT (CONVENTIONS.md 3) |
| `inf_cle` | 0.5–1.5, log-uniform | Clearance from both latent states |
| `inf_non`, `inf_asy` (recent latent) | 0.04–0.45; 0.01–0.60, log-uniform | The progression hazard at infection, as tb_LAI_TPT's rates at τ = 0 |
| `k_latent` | 0.1–12 | Replaces tb_LAI_TPT's `k_prog` (0–12): both set how fast early risk wanes. Above 0, because 0 switches off v4's remote state |
| φ (remote ÷ recent progression) | 0–0.1 | CONVENTIONS.md 2's upper range (user, 2026-11-06) |
| `non_rec`, `non_asy` | 0.14–0.23; 0.21–0.28 | |
| `asy_non`, `asy_sym`, `sym_asy` | 0.70–2.03; 0.56–0.99; 0.46–0.85 | |
| `sym_dead` | 0.20–0.50 | Constant; tb_LAI_TPT's 1995 range |

13 parameters for 6 targets, so several will be weakly identified; the posterior shows which.

## Method

**IMIS** (incremental mixture importance sampling): importance sampling from the prior, then repeatedly adding multivariate-normal components at the highest-weight points until the effective sample size is large enough. It suits a deterministic model with a posterior much narrower than its prior box, and the user knows it.

- **Settings:** 10,000 initial prior draws; 1,000 draws per added component, centred on the highest-weight point, with the weighted covariance of its 1,000 nearest points. Stop when a resample of 3,000 would contain (1 − 1/e)·3,000 ≈ 1,896 unique points in expectation, or after 150 iterations (re-identification) or 300 (real targets).
- **A fast copy of the ODE:** v4's solver takes about a second per equilibrium, too slow for IMIS's hundreds of thousands of evaluations. A standalone version of the same equations (our own code, outside Starsim's solver; exp 19), with its equilibrium found by root-finding, is checked against v4's ODE (equilibrium shares and calibrated targets within about 1e-6) before it is used.
- **Coverage check first:** 5,000 prior draws (`run.py coverage`). For each output, the share inside its target interval, and the share of draws inside every interval at once.
- **Re-identification first:** IMIS on synthetic targets from known parameters, to check it recovers them, before the real targets.
- **Outputs:** the posterior sample with weights, its effective sample size, marginals and pairwise plots against the priors, and the fit to each target. A small weighted set of parameter sets is then drawn for the ODE↔agents comparison.

## Files

- `calib.py`: priors, the model (`simulate`: `tbode4.fastode` at equilibrium plus the progression chain), the likelihood, and IMIS.
- `run.py`: `coverage` → `coverage.obj`, `coverage.csv`, `coverage.log`; `reid` and `real` → `reid.obj`/`real.obj` (every point, its outputs, log-likelihood and weight; the IMIS components and log), `*_posterior.csv`, `*_fit.csv`, `*_marginals.png`, `*_fit.png`; `real_seed2` → the real calibration again with another seed, to check convergence (`real_seed2.obj`, `real_seed2_posterior.csv`, `real_seed2_fit.csv`, `real_seed2.log`); `plots` redraws the summaries and figures. Writes `versions.json`.
- `SUMMARY.md`: findings, written after running.
