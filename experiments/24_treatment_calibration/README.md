# 24: Calibrating the v4 ODE with treatment, time-varying TB mortality and over-detection

**Status:** done; see [SUMMARY.md](SUMMARY.md).
**Track:** v4 (`v4/.venv`; tbsim `origin/v4` @ `dcfb522`, starsim `origin/v4` @ `80c81f53`).

## Question

Calibrate the v4 ODE, now with treatment, to tb_LAI_TPT's targets that apply to a model without age, demography or a household network. With treatment, the targets that exp 20 dropped or couldn't use come into reach: mortality ÷ incidence, prevalence ÷ notifications, and the trends over time.

The calibration includes two features that tb_LAI_TPT calibrates and that agents haven't yet been checked on (user, 2026-11-09):

- **TB mortality that changes over time:** `sym_dead` falls from 1995 to 2024, as tb_LAI_TPT lets it.
- **Over-detection:** reported notifications are (1 + `over_detection`) × treatment starts. It affects only the notifications compared with data, not the dynamics.

A later experiment checks that agents still match the ODE with both, at posterior draws, as exp 21 did without treatment.

## Model

- **The ODE:** `tbode4.fastode`, with exp 23's treatment.
  - People with symptomatic TB start treatment at θ = −ln(1 − coverage).
  - Whether a course will cure is drawn at the start, with probability efficacy (adherence 1, as in tb_LAI_TPT).
  - Failing courses face TB death, then return to symptomatic TB.
  - A share `p_relapse` of cures relapse, at rate 1/1.5 per year.
  - Treatment is a fixed 6-month course, as 12 Erlang stages: within 0.1% of the exact fixed course's prevalence at equilibrium (exp 23), and within 0.04% over 1995–2040 at one parameter set (checked against v4's ODE with its exact fixed course).
  - The rest is exp 20's model: the two-state latent model, homogeneous mixing, and births replacing every death at μ = 1/70.
- **Over time:** everything is constant before 1995, so each draw starts at its untreated endemic equilibrium at the end of 1994, and runs to 2040.
  - **Coverage:** 0 at the end of 1994, linear to `tx_covg_2010` in 2010 and `tx_covg_2024` in 2024, then constant.
  - **Efficacy:** `efficacy_1995` from 1995, linear to `efficacy_2024`, then constant.
  - **`sym_dead`:** `sym_dead_1995` until 1995, linear to `sym_dead_2024`, then constant. It applies to symptomatic TB and to failing courses alike.

  These are tb_LAI_TPT's trajectories (`utils.build_sim`, `TimeVaryingSymDead`).
- **Outputs, as tb_LAI_TPT defines them** (`hm_common.extract_targets_one_sim`):
  - **Point values on 1 January:** prevalence (asymptomatic + symptomatic TB) and the share asymptomatic.
  - **Totals over each calendar year:** incidence (first entry to asymptomatic TB per episode, plus relapses; CONVENTIONS.md 5), notifications (treatment starts × (1 + `over_detection`)) and TB deaths (from symptomatic TB and on failing courses).

  These yearly totals are integrated exactly along each trajectory.
- **Cost:** about 90 ms per draw (against 10 ms at equilibrium in exp 20), so an IMIS run of about 150,000 draws takes about an hour on 4 cores.
- **v4's own ODE agrees:** v4's `TBTreat` with the same inputs as `ss.timeseries`, including `sym_dead`, matches `fastode`'s prevalence over 1995–2040 exactly with 12 Erlang stages.

## Targets and likelihoods

The means are tb_LAI_TPT's. The σ are its data σ (`hm_common.calibration_sigmas`, run at its calibration's versions), with ratio σ combining numerator and denominator as it does. Likelihoods are as in exp 20: Binomial for the survey proportions, log-Normal for the rest, since the ODE is deterministic.

| Target | Definition | Data, year | Mean | Likelihood | 95% interval | In exp 20 |
|---|---|---|---|---|---|---|
| Prevalence | Point prevalence, 1 January 2020 | National survey, 2020 (15+) | 310/100k | Binomial(77, 24,815) | 247–385 | Yes |
| Prop_Asymptomatic | Asymptomatic ÷ prevalence, 2020 | Same survey | 0.654 | Binomial(123, 188) | 0.584–0.720 | Yes |
| Prev_to_Inc | Prevalence 2020 ÷ incidence in 2020 | Survey ÷ WHO, 2020 | 1.57 y | log-Normal, σ = 0.136 | 1.20–2.05 | Yes |
| Mort_to_Inc | TB deaths ÷ incidence, 2024 | WHO, 2024 | 0.113 | log-Normal, σ = 0.127 | 0.088–0.144 | Reported only |
| Prev_to_Notif | Prevalence 2020 ÷ reported notifications in 2020 | Survey ÷ notifications, 2020 (116.7/100k) | 2.66 y | log-Normal, σ = 0.125 | 2.08–3.40 | **New** |
| Inc_Ratio | Incidence 2024 ÷ 2015 | WHO (187.3 ÷ 237.2/100k) | 0.790 | log-Normal, σ = 0.112 | 0.634–0.983 | **New** |
| Mort_Ratio | TB deaths 2024 ÷ 2015 | WHO (21.1 ÷ 32.0/100k) | 0.658 | log-Normal, σ = 0.120 | 0.520–0.832 | **New** |
| Prev_2040_Ratio | Prevalence 2040 ÷ 2020 | tb_LAI_TPT's anchor (a judgement, not data) | 0.806 | log-Normal, σ = 0.190 | 0.556–1.170 | **New** |
| ProgressionY1 | Latent cohort reaching asymptomatic TB by 1 year | Ferebee, Sutherland | 0.0523 | log-Normal, σ = 0.094 | 0.0435–0.0629 | Yes |
| ProgressionRatio | The same at 1 ÷ 10 years | Same | 0.663 | log-Normal, σ = 0.071 | 0.577–0.762 | Yes |

**Left out:**

- **Prop_Transmission_HH:** there's no household network.
- **The under-15 targets:** there are no ages.

**Accepted mismatches:**

- **No ages:** the ODE's all-ages prevalence is compared with the 15+ survey, as in exp 20.
- **No demography:** India's population grew and aged over 1995–2040, which the constant-population ODE can't show (user, 2026-11-09). That affects the trend targets most.
- **Coverage for children:** tb_LAI_TPT treats children at 0.22 × adults' coverage. With no ages, one coverage applies to everyone.

## Parameters and priors

21 parameters for 10 targets. Exp 20's 13 natural-history parameters carry over, with `sym_dead` as its 1995 value. tb_LAI_TPT's treatment, mortality-trend and over-detection parameters are added with its bounds and its parameterisations, which keep each trend in its direction (`run_history_matching.parameter_bounds`). Uniform, or log-uniform where marked.

| Parameter | Prior | Note |
|---|---|---|
| `beta`, `trans_asymp`, π, `inf_cle`, `inf_non`, `inf_asy`, `k_latent`, φ, `non_rec`, `non_asy`, `asy_non`, `asy_sym`, `sym_asy` | As exp 20 | `beta`, `inf_cle` and `inf_asy` log-uniform |
| `sym_dead_1995` | 0.20–0.50 | Exp 20's `sym_dead` |
| `sym_dead_frac` | 0–1 | `sym_dead_2024` = `sym_dead_1995` − frac × (`sym_dead_1995` − 0.05) |
| `tx_covg_2024` | 0.50–0.95 | Annual coverage of symptomatic TB |
| `tx_covg_2010_frac` | 0–1 | `tx_covg_2010` = frac × `tx_covg_2024` |
| `efficacy_1995` | 0.30–0.95 | Chance a course cures |
| `efficacy_frac` | 0–1 | `efficacy_2024` = `efficacy_1995` + frac × (0.95 − `efficacy_1995`) |
| `p_relapse` | 0.02–0.10 | Share of cures that relapse |
| `over_detection` | 0–0.20 | Notifications only |

## Method

As exp 20 ([CONVENTIONS.md 8](../../CONVENTIONS.md#8-calibrating-the-ode-a-standalone-copy-at-equilibrium-with-imis)):

1. **Coverage check:** 5,000 prior draws. For each target, the share inside its interval, and the share of draws inside every interval at once.
2. **Re-identification:** IMIS on synthetic targets from a known point, with the real targets' σ and survey sizes.
   - **The truth:** exp 21's set 1 natural history (the starter script's), with tb_LAI_TPT's best-fit treatment, `sym_dead` trend and `over_detection`. Its `sym_dead_1995` is set 1's 0.396.
   - **What to check:** whether IMIS recovers it, and which parameters the targets identify.
3. **Real targets,** and again with a second seed, to check convergence.

**IMIS settings:** as exp 20. 10,000 initial draws, 1,000 per component, and a stop at (1 − 1/e) × 3,000 expected unique points in a resample of 3,000, or after 300 iterations.

**Outputs:**

- the posterior sample with weights, its effective sample size, and marginals against the priors;
- the fit to each target;
- the posterior's trajectories of prevalence, incidence, notifications and mortality over 1995–2040, against the data;
- a weighted draw of parameter sets for the agents check that follows.

## Files

- `calib.py`: priors, the model (`trajectories`, `simulate`), the likelihood and IMIS.
- `run.py`: `coverage`, `reid`, `real`, `real_seed2`, `plots`. Each saves its raw results (`*.obj`) before summarising, and writes `versions.json`.
- `SUMMARY.md`: findings, written after running.
