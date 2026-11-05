# 14: Summary

> **Superseded by [exp 15](../15_two_state_fit/SUMMARY.md)** (user decision, 2026-11-04).
>
> - **Why:** this experiment fixes progression from the remote latent state at 0. That makes the two-state model effectively a single latent state with raised clearance, and we don't want to fix it at 0.
> - **Instead:** exp 15 fits a genuinely two-state model with remote progression over a range of values, and includes 0 as one point in that range.

**Track:** v4 environment (tbsim `origin/v4` @ `dcfb522`, starsim `origin/v4` @ `80c81f53`). The calculation is numpy/scipy only and deterministic.

## Question

Is tbsim v4's two-state latent model, with `k_latent` refitted and no remote progression (option C in an earlier draft of the [options note](../../proposals/02_latent_progression/options.md), now set aside there), close to tb_LAI_TPT's front-loaded progression r·e^(−kτ) across its whole posterior, not just at the best fit (set B)? Does its reinfection rule (remote → recent; a recent infection unchanged) still add the same progression as an exact clock reset?

The models compared:

- **A, exact:** Poisson-mixture latent classes, which reproduce r·e^(−kτ) exactly.
- **B:** two states with `k_latent` = k.
- **C:** two states with `k_latent` refitted so that ProgressionY1 matches A.

They are compared on a latent cohort for each of the 345 weighted posterior samples (see README for the design).

## Findings

**1. C is within about 1% of the exact model for most of the posterior, and never more than 2% off.**

| Largest gap from the exact model, 3 months–50 years (share ever reaching asymptomatic TB, \|structure ÷ exact − 1\|) | Weighted median [10%, 90%] | Max over 345 samples | Posterior weight with gap > 1% / > 2% / > 5% |
|---|---|---|---|
| B: `k_latent` = k | 3.8% [2.7%, 4.7%] | 6.3% | 1.00 / 0.99 / 0.06 |
| C: `k_latent` refitted | 0.85% [0.56%, 1.03%] | 2.0% | 0.14 / 0.004 / 0 |

Set B (row 0) is typical: 0.94% for C. The gap has the same shape in every sample (figure, middle panel):

- **Before 1 year:** C is slightly high, by up to 1.3% at 3 months.
- **At 1 year:** C is exact, by construction.
- **After 1 year:** C is slightly low, levelling off by about 5 years. Lifetime progression is low by 0.19% [0.03%, 0.96%], and by 2.0% at worst.

The largest gaps are in samples with r₀/k ≈ 0.2–0.5 and `k` ≈ 0.4–1.6 per year, where front-loading is moderate. They carry 7% of the posterior weight (r₀/k > 0.2). Where `k` is very small (down to 0.014), both models are nearly constant-rate, so the gap is small again (0.05%).

**2. The refit is close to a fixed rule:** `k_latent` = 0.951 k [0.937, 0.962] (max 0.973). Row 0's 5.533 = 0.947 × 5.846.

**3. The progression targets are unaffected.** ProgressionY1 matches A by construction. ProgressionRatio, B(1)/B(10), differs from A by at most 0.014 in any sample (0.021 for B). Its weighted median is 0.707 [0.631, 0.771]; the targets are 0.0523 / 0.663.

**4. Almost all progression after year 1 comes through non-infectious TB, not latent infection.**

| Exact model, no reinfection | Set B (row 0) | Weighted median [10%, 90%] |
|---|---|---|
| Progression between years 1 and 10, B(10) − B(1) | 0.025 | 0.021 [0.013, 0.031] |
| Share from people still latent at year 1 | 0.5% | 2% [0.05%, 38%] |
| Latent progression hazard at 1 year, r₀·e^(−k), per year | 0.002 | 0.007 [0.0003, 0.07] |

- **What the targets require:** later progression, about a third of 10-year progression. They measure ever reaching asymptomatic TB, from the Ferebee and Sutherland data.
- **Where the calibration put it:** in most of the posterior, in people who reached non-infectious TB early and move on to asymptomatic TB at `non_asy` over the following years.
- **The exception:** where `k` is small, late latent progression carries more. The share from latent infection is above 10% in 28% of the posterior weight, and above 25% in 16%.
- **Why the two-state model still meets the long-term targets with no progression from the remote state:** the route through non-infectious TB is unchanged.

**5. The reinfection rule adds the same progression as a clock reset, to within 0.0004 at set B's force of infection.**

| Progression added by latent reinfection by 30 years (share of the cohort) | λ = 0.12/y (set B at 500/100k) | λ = 0.60/y (5×) |
|---|---|---|
| A, exact clock reset: weighted median [10%, 90%] | 0.0032 [0.0014, 0.0055] | 0.029 [0.013, 0.048] |
| C − A: weighted median [10%, 90%] | 0.00001 [0.000005, 0.00005] | 0.00013 [0.00006, 0.0005] |
| C − A: max over samples | 0.0004 | 0.0024 |

- **Relative to what the reset adds,** C adds 0.4–0.5% more (weighted median), and 3% more at the 90th percentile.
- **Larger relative differences** (up to 21%) occur only where `k` is small. There front-loading is weak and a reset adds almost nothing: 0.00002 of the cohort in the 21% case.
- **By 10 years** the numbers are about three-quarters as large.

**6. Multiple strains.** In every candidate, as in v0.10, the clock belongs to the person, and superinfection carries no progression penalty in tb_LAI_TPT. So what C's rule loses with strains depends only on the total force of infection on latent people. Finding 5 is therefore also the answer for strains, at the total λ.

**7. Size of the exact model.** The classes need 6 [6, 7] states at a truncation of 1e-6 and 8 [8, 9] at 1e-10. Where `k` is very small they need up to 38 and 47.

## What this means for the decision

- **Option C holds across the posterior.** The ~1% found for set B is typical. Where calibrations put appreciable weight, C stays within ~1% of r·e^(−kτ), and the worst sample (2%) is in the 7% of posterior weight with moderate front-loading. Its reinfection rule matches a clock reset to 0.0004 or better at realistic forces of infection, with or without strains.
- **Whether 2% is "consistent enough" is the remaining judgment.** The gap is in lifetime progression, after fitting to ProgressionY1 only. A refit to a different criterion (e.g. the largest gap over time) would split it between early and late times. It wasn't tried.
- **With no progression from the remote state, a remote latent infection behaves almost like a cleared one.** Neither progresses, and in tb_LAI_TPT both are reinfected at π·λ. That follows from tb_LAI_TPT's fitted `k` (finding 4), not from the two-state approximation: the exact model's class 0 is the same. The remaining differences:
  - remote infections count as latent infection;
  - they carry strains;
  - TPT acts on them.

  How much late latent progression there should be is part of the open question on late reactivation.

- **Option B is 3–5% low for most of the posterior,** so a proposal for C should give the refit rule (≈ 0.95 k, or a root-find on ProgressionY1). Without it, tb_LAI_TPT's `k` can't be used directly.

*How the runs came about: the first run computed ProgressionRatio as B(1)/B(5). It was rerun with tb_LAI_TPT's definition, B(1)/B(10), and with finding 4's split added. No other numbers changed.*

## Files

- `run.py`; `run.log` (console output)
- `samples.csv` (per sample); `summary.json`; `outputs.obj` (data and cohort curves behind the figure); `posterior_progression.png`
- `versions.json`
