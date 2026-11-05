# 15: Summary

> **Superseded by [exp 16](../16_two_state_vs_lai_tpt/README.md)** (user decision, 2026-11-04).
>
> - **Why:** this experiment fits each posterior sample's own exact-model values and compares against the exact r₀·e^(−kτ) model. But tb_LAI_TPT computed its progression targets with v0.10's monthly ABM, which overstates early progression, so the exact model isn't what tb_LAI_TPT ran (see the note below the first table).
> - **Instead:** exp 16 fits the target means directly, sweeps φ, and compares the resulting dynamics with tb_LAI_TPT as run.

**Track:** v4 environment (tbsim `origin/v4` @ `dcfb522`, starsim `origin/v4` @ `80c81f53`). The calculation is numpy/scipy only and deterministic (`tbode4.latent`, checked by `v4/tests/test_latent.py`).

## Question

In tbsim v4's two latent states, recent infections progress and become remote at `k_latent`. Can a version that is genuinely two-state, with remote infections also progressing (at a fraction φ of the recent rates), do two things within tb_LAI_TPT's priors?

- **(a)** Fit the progression targets: ProgressionY1 = B(1) = 0.0523 and ProgressionRatio = B(1)/B(10) = 0.663, where B(t) is the share of a latent cohort that has ever reached asymptomatic TB.
- **(b)** Reproduce tb_LAI_TPT's progression by time since infection.

And how much does fitting only the progression targets shift cumulative progression, which would pull the model off the epidemiological targets?

**Design in brief:**

- **Samples:** for each of tb_LAI_TPT's 345 posterior samples, the sample's `inf_cle`, `non_rec`, `non_asy`, π and NON : ASY split are kept.
- **Fitted:** the recent progression rate r₀ and `k_latent`. They are fitted to reproduce that sample's own ProgressionY1 and ProgressionRatio, at φ = 0, 0.003, 0.01, 0.03, 0.1 and 0.3.
- **Bounds:** the agreed priors.
- **Reference:** each sample's exact r·e^(−kτ) model (see README).
- **Statistics:** weighted by `sample_weight`; weighted medians [10%, 90%].

## Findings

**1. (a) Remote progression up to about 3% of recent progression fits the targets as well as tb_LAI_TPT's posterior does. Beyond that, the fit worsens quickly.**

| φ (remote ÷ recent progression) | 0 | 0.003 | 0.01 | 0.03 | 0.1 | 0.3 |
|---|---|---|---|---|---|---|
| Posterior weight whose best fit is within 2σ of both target means (exact model: 68%) | 68% | 69% | 73% | 75% | 52% | 25% |
| Weighted median of the larger \|z\| to the target means (exact model: 1.36) | 1.36 | 1.35 | 1.28 | 1.40 | 1.71 | 3.37 |
| Posterior weight that reproduces the sample's own Y1 and ratio exactly | 100% | 55% | 38% | 22% | 10% | 1% |

z is log(model ÷ target mean) ÷ σ, with σ = 0.094 (Y1) and 0.071 (ratio) in log space, as in tb_LAI_TPT's calibration.

**What "exact model" means here, and a limitation.**

- **The term:** tb_LAI_TPT's own progression model, with hazard r₀·e^(−kτ) and each posterior sample's parameters, computed exactly with the Poisson-mixture classes (no time-step error).
- **The limitation:** that isn't what tb_LAI_TPT's calibration saw. It computed the progression observables by running v0.10's `tbsim.TB` on a cohort at monthly `dt` (`cohort_progression.py`). Exp 08 found that this overstates early progression by 10–19%, because the hazard at the start of each step is applied for the whole step.
- **So:** for the same parameters, the exact calculation gives a lower ProgressionY1 and a higher ratio than the calibration did. The 68% within 2σ is not the calibration's own fit. Nor is the exact model the dynamics tb_LAI_TPT actually ran, so comparing the two-state model against it compares against a reference tb_LAI_TPT never used.

- **The ratio is what limits the fit.** Remote progression adds progression after year 1, which lowers ProgressionRatio; Y1 is met throughout. For set B (row 0), the ratio's z goes from +0.24 at φ = 0 to −0.58 at 0.01, −1.65 at 0.03 and −3.8 at 0.1.
- **Small φ helps, at first.** tb_LAI_TPT's posterior mostly sits above the ratio target (weighted median 0.707 against 0.663). So up to φ ≈ 0.03, remote progression moves it toward the target, and slightly more of the posterior lands within 2σ than with the exact model.
- **The binding constraint is the model's structure, not the priors.** Only `inf_non` ≤ 0.45 binds, for 22% of the weight at φ = 0.003 and for set B there; no bound binds at φ ≥ 0.1. The structure has no room for later progression because the sample's route through non-infectious TB already supplies the progression after year 1 that the targets need (exp 14): in the exact model, 98% of it comes through non-infectious TB [62%, 99.95%]. Freeing `inf_cle`, `non_rec` and `non_asy` within their priors could make room, but it wasn't tried.
- **For set B, φ = 0.03 means remote rates of `inf_non_remote` = 0.0075 and `inf_asy_remote` = 0.0050 per year,** with the recent rates 0.25 / 0.17 and `k_latent` = 3.3.

**2. (b) Up to about a year after infection, the two-state model reproduces the dynamics. After that it has a step down to a constant floor, and the exact hazard declines to zero.**

- **The cohort curve B(t)** is matched at 1 and 10 years by construction, and is close in between for small φ. For samples that reproduce their own targets exactly, the largest gap over 3 months–50 years is 0.5% [0.08%, 1.4%] at φ = 0, and 1.0% [0.3%, 5.0%] at φ = 0.003.
- **The latent progression hazard** (figure, third panel):
  - **Early:** the two models are close in the first months.
  - **Later:** the exact hazard r₀·e^(−kτ) keeps falling, to about 10⁻¹³ by 5 years for set B (weighted median over the posterior about 10⁻¹⁰). The two-state hazard levels off at φ·r₀: 0.002 / 0.006 / 0.012 per year at φ = 0.003 / 0.01 / 0.03 for set B.
  - **So:** this is the lifelong reactivation that makes the model genuinely two-state, and it is the part of tb_LAI_TPT's dynamics the two-state model does not reproduce.
- **Where later progression comes from changes.**

  | Set B | Exact model | φ = 0.003 | φ = 0.01 | φ = 0.03 | φ = 0.1 |
  |---|---|---|---|---|---|
  | Share of progression between years 1 and 10 from people still latent at year 1 | 0.5% | 6% | 14% | 27% | 45% |

  The rest comes from people already in non-infectious TB.

**3. Side effects: beyond φ ≈ 0.003, fitting only the progression targets raises cumulative progression by more than ~1–2%.** The targets stop at 10 years, and remote progression continues beyond them.

| Best fits, all samples: two-state ÷ exact, weighted median [10%, 90%] | φ = 0.003 | 0.01 | 0.03 | 0.1 | 0.3 |
|---|---|---|---|---|---|
| Lifetime progression, B(50) | 1.000 [1.000, 1.018] | 1.015 [1.000, 1.06] | 1.055 [1.000, 1.15] | 1.16 [1.00, 1.35] | 1.36 [1.10, 1.71] |
| Ever reached asymptomatic TB by 30 years, with reinfection at λ = 0.12/y | 1.001 [1.000, 1.015] | 1.013 [1.000, 1.05] | 1.046 [1.000, 1.12] | 1.13 [1.00, 1.29] | 1.29 [1.09, 1.58] |

**How the reinfection row is defined:**

- **Cohort:** new infections followed for 30 years at a constant force of infection λ (the infection rate per year of a fully susceptible person), with background deaths at 1/70 per year.
- **Reinfection of people no longer infected:** at π·λ (π ≈ 0.24, the sample's reinfection risk ratio).
- **Reinfection of latent people:** also at π·λ. In the exact model this resets the clock; in the two-state model, remote goes back to recent and recent is unchanged.
- **λ = 0.12/y** is set B's force of infection at its 500/100k equilibrium. Results at 0.60/y (5×, a stress test) are in `fits.csv`.

- **Set B:** lifetime progression rises by 1.8% / 6.3% / 15% at φ = 0.003 / 0.01 / 0.03, and by 1.5% / 5.1% / 12% with reinfection.
- **Time spent latent** barely changes (within about 2% up to φ = 0.03; 0.983 for set B), because clearance dominates exit from latent infection.
- **Samples that reproduce their own targets exactly** have side effects of about 1.000. These are samples whose exact model already has late latent progression (small `k`).

## What this means

- **(a) Yes, for small φ.** A genuinely two-state model fits the progression targets as well as tb_LAI_TPT's posterior, within its priors, with remote progression up to about 3% of recent progression. That is remote rates of order 0.005–0.01 per year for set B.
- **(b) Not after the first year or so,** as anticipated. The two-state model replaces a hazard that declines to zero with a step to a constant floor. So it adds lifelong reactivation, and moves a substantial share of the progression after year 1 from non-infectious TB to latent infection.
- **The cost is cumulative progression beyond the 10-year horizon of the targets:** +1.5% at φ = 0.01 and +5% at φ = 0.03 (weighted medians), and more in set B. Near threshold, prevalence amplifies a change like this several-fold, so above φ ≈ 0.003 a two-state model fitted only to the progression targets would need recalibration against the epidemiological targets (`beta` at least).
- **Not settled here:**
  - whether freeing `inf_cle`, `non_rec` and `non_asy` within their priors allows larger φ. This would shift progression from non-infectious TB to remote latent infection, and would also affect the epidemiological targets.
  - what value of φ the long-term reactivation literature supports. The targets only constrain progression to 10 years.

## Files

- `run.py`; `run.log`
- `fits.csv` (per sample and φ); `summary.json`; `outputs.obj` (data and curves behind the figure); `two_state_fit.png`
- `versions.json`
- Shared code: `v4/tbode4/latent.py`; tests: `v4/tests/test_latent.py`

*How the runs came about: the first run printed ratios at 3 significant figures and hazard ratios at 5 years, which were uninformative because the exact hazard there is about 10⁻¹⁰. It was rerun with more digits, absolute hazards, and fit to the target means. The fits and the curves didn't change.*
