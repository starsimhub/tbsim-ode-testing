# 20: Summary

**Track:** v4 (tbsim `origin/v4` @ `dcfb522`, starsim `origin/v4` @ `80c81f53`; `versions.json`).

**Setting:**

- **Model:** v4's TB ODE with the two-state latent model and homogeneous mixing, no treatment, births replacing deaths. It runs as the standalone copy (`tbode4.fastode`, exp 19) at its endemic equilibrium.
- **Targets:** tb_LAI_TPT's targets that don't need treatment, trends or a household network (README table).
- **Method:** IMIS over 14 parameters with tb_LAI_TPT's priors (two-state version; φ up to 0.1).

## Coverage check: mortality ÷ incidence is out of reach without treatment

5,000 prior draws (`coverage.csv`):

- **Mortality ÷ incidence:** across 90% of draws it's 0.18–0.46, against a target of 0.088–0.144 (WHO, 2024). Only 0.6% of draws are inside, and none below.
  - **What it would take:** the lowest values need fast improvement from asymptomatic TB, slow progression to symptoms, fast remission and low `sym_dead`, all at the edges of their priors at once.
  - **The conflict:** the share asymptomatic needs high `sym_dead`, so the two targets pull against each other.
  - **Why:** treatment removes symptomatic cases without killing them, which no setting of untreated natural history can do.
  - **So it's dropped from this calibration** (user, 2026-11-06) and comes back with treatment. It is still reported.
- **The other five targets** are each reachable within the priors (4–51% of draws inside each interval).
- **No endemic equilibrium:** 47% of draws, mostly low `beta`. Those draws get zero likelihood.

## Re-identification: the pipeline recovers a known point, and φ isn't identified

**Setup:** synthetic targets from set B with φ = 0.01 (exp 17's model), with the real targets' sample sizes and σs. IMIS ran 160,000 model runs (150 iterations, effective sample size 697; `reid_posterior.csv`, `reid_marginals.png`).

- **The pipeline finds the best fit:** its best log-likelihood (−1.932) is essentially the truth's own (−1.918), and the posterior reproduces every synthetic target.
- **The truth is inside the 95% interval for 13 of 14 parameters.** `beta` is fairly well identified: its interval spans 22% of the prior's width. Most others are weakly identified, spanning 50–90%. Five targets can't pin down 14 parameters.
- **φ is the exception:** truth 0.01, posterior 0.015–0.094.
  - **The likelihood is flat in φ:** the best fit in each band of φ from 0.01 to 0.1 lies between −1.93 and −1.99; below 0.01 it's −2.39.
  - **Why the posterior leans higher:** more of the other parameters' prior volume fits at larger φ.
  - **So the targets can't identify φ,** as exp 16 found from the progression targets alone. Its posterior reflects the prior volume, not evidence for larger φ.

## Calibration to the real targets

**IMIS:** it stopped by its rule after 218 iterations: 228,000 model runs, effective sample size 1,529, and 1,902 expected unique points in a resample of 3,000 (`real.log`).

**Fit:** all five targets are reproduced within their intervals (`real_fit.csv`, `real_fit.png`).

| Target | Target mean (95% interval) | Posterior predictive, median (2.5–97.5%) |
|---|---|---|
| Prevalence (per 100k) | 310 (247–385) | 307 (233–391) |
| Share asymptomatic | 0.654 (0.584–0.719) | 0.631 (0.580–0.676) |
| Prevalence ÷ incidence (years) | 1.57 (1.20–2.05) | 1.60 (1.34–2.01) |
| ProgressionY1 | 0.0523 (0.0435–0.0629) | 0.0529 (0.0433–0.0635) |
| ProgressionRatio | 0.663 (0.577–0.762) | 0.666 (0.584–0.752) |
| Mortality ÷ incidence (not fitted) | 0.113 (0.088–0.144) | **0.234 (0.154–0.342)** |

- **The share asymptomatic sits on the low side,** as the coverage check suggested; without treatment, symptomatic cases aren't removed.
- **Mortality ÷ incidence is about twice the 2024 data,** which is what untreated natural history implies. It's roughly the WHO's 2000 ratio (0.29).

**Parameters** (`real_posterior.csv`, `real_marginals.png`): posterior median (2.5–97.5%), and the interval's share of the prior's width.

| Parameter | Prior | Posterior | Share of prior |
|---|---|---|---|
| `beta` | 2–100 | 19.4 (10.7–34.1) | 0.24 |
| π (reinfection risk) | 0.2–1.0 | 0.34 (0.21–0.66) | 0.56 |
| `inf_asy` | 0.01–0.60 | 0.31 (0.16–0.49) | 0.56 |
| `asy_non` | 0.70–2.03 | 1.64 (1.18–1.97) | 0.59 |
| `asy_sym` | 0.56–0.99 | 0.65 (0.56–0.82) | 0.60 |
| `k_latent` | 0.1–12 | 6.7 (2.3–11.6) | 0.78 |
| `sym_asy` | 0.46–0.85 | 0.72 (0.53–0.84) | 0.80 |
| `sym_dead` | 0.20–0.50 | 0.41 (0.25–0.49) | 0.81 |
| `inf_non` | 0.04–0.45 | 0.27 (0.09–0.43) | 0.84 |
| `inf_cle` | 0.5–1.5 | 0.82 (0.53–1.40) | 0.87 |
| φ | 0–0.1 | 0.054 (0.008–0.095) | 0.87 |
| `trans_asymp` | 0.5–1.0 | 0.73 (0.53–0.97) | 0.88 |
| `non_rec`, `non_asy` | 0.14–0.23; 0.21–0.28 | 0.19 (0.14–0.23); 0.24 (0.21–0.28) | 0.91 |

- **What the data constrain:** `beta`, π and the progression rates. The asymptomatic/symptomatic cascade is pushed towards keeping people asymptomatic: fast `asy_non`, slow `asy_sym`, fast `sym_asy`, and high `sym_dead`. That's how untreated natural history gets close to the share asymptomatic.
- **What they don't:** the non-infectious cascade, `trans_asymp` and φ are close to their priors. With five targets, that's expected.

## Convergence

- **Stopping rule:** Raftery and Bao's, the standard one for IMIS. A resample of 3,000 points should contain at least (1 − 1/e)·3,000 ≈ 1,896 unique points in expectation. The real run met it at iteration 218 (1,902), with an effective sample size of 1,529.
- **The weights aren't concentrated:** the largest single weight is 0.43%; the top 10 points carry 3% and the top 1,000 carry 55%.
- **No local optimum:** the best log-likelihood possible, every target matched exactly, is −1.638, and IMIS reached −1.650. With 14 parameters and 5 targets, many parameter sets fit exactly, and IMIS found that region.
- **A second run with another seed gives the same posterior** (`real_seed2.obj`, `real_seed2_posterior.csv`; 213,000 runs, stopped at iteration 203, effective sample size 1,553).
  - **The fit to every target is the same,** with medians within 0.6%.
  - **The constrained parameters agree closely.** The medians of `beta` (19.45, 19.44), π, `k_latent`, φ and `inf_asy` differ by at most 2% of the prior's width.
  - **The weakly identified parameters agree to within about 5–13% of their prior's width,** in medians and in interval ends: `inf_cle` 0.82 against 0.95, `sym_dead` 0.41 against 0.39. That's Monte Carlo error in broad, nearly flat marginals at an effective sample size of about 1,500, not separate modes: both runs give the same shapes.
  - **For sharper summaries of those parameters,** pool the two runs (an effective sample size of about 3,000), or run IMIS to a larger target.

## What this means

- **There is now a calibrated parameter set with realistic burden for the ODE↔agents comparison:** a weighted posterior sample (`real.obj`). The natural history partly absorbs the effects of treatment (the asymptomatic/symptomatic cascade above), and mortality per case is untreated-level.
- **φ needs data these targets don't have.** It stays open, as CONVENTIONS.md 2 anticipates, for a calibration with more targets.
- **The pipeline is ready for the recalibration with treatment:** add treatment to the model, bring mortality ÷ incidence back, and rerun the coverage check, re-identification and IMIS.
