# 16: Summary

**Track:** both environments.

- **v4** (tbsim `origin/v4` @ `dcfb522`, starsim `origin/v4` @ `80c81f53`): the two-state fits and the exact reference, from `tbode4.latent`'s chain. The chain matches v4's own `tbsim.TB(solver='ode')` cohort to within 2.3e-6 in every state's share (`param00`, arm 1, φ = 0.03).
- **v0.10** (tbsim v0.10.3 @ `b92d7e8`, starsim v3.6.1 @ `3521adb4`): tb_LAI_TPT as run, an explicit exception to the no-v0.10-arms rule (user, 2026-11-04).
  - **It reproduces tb_LAI_TPT exactly:** at its own settings (20k agents, seeds 0 and 1), our analyzer's curve B equals `cohort_progression.FirstPassage`'s to 0. ProgressionY1 = 0.05390 and ProgressionRatio = 0.63825, as the calibration saved.

## Question

tbsim v4's two-state latent model has a recent state and a remote state, which progresses at a fraction φ of the recent rates. Within tb_LAI_TPT's priors:

1. Can it meet tb_LAI_TPT's progression targets?
2. What dynamics result, and how do they compare with tb_LAI_TPT as run (v0.10 `tbsim.TB`, monthly `dt`)?

**Design in brief:**

- **φ:** swept over 0, 0.003, 0.01, 0.03, 0.1 and 0.3, not calibrated.
- **Arm 1:** the active-disease cascade fixed at each sample's values.
- **Arm 2:** `inf_cle`, `non_rec` and `non_asy` also freed within their priors.
- **Fitting:** to the target means (ProgressionY1 = B(1) = 0.0523, ProgressionRatio = B(1)/B(10) = 0.663; log-space σ 0.094 / 0.071).
- **Transmission:** off. Outputs come from a no-reinfection cohort, for the 12 parameter sets tb_LAI_TPT's face-validity checks use.
- **The 12 sets** are the 12 most likely posterior samples (the file is sorted by log-likelihood, best first), not a representative draw. Their `sample_weight` is renormalised over the 12, as tb_LAI_TPT does. Together they hold 34% of the posterior's weight, and the top two hold 45% of the renormalised weight (25% and 20%). So finding 2's weighted means lean heavily on those two.
- **Which sets the results cover:** the tables are across all 345 samples (finding 1) or all 12 sets (finding 2). Where a line gives one set's numbers, it is **set B**: `param00`, row 0 of tb_LAI_TPT's posterior and its best fit. Earlier experiments (08–13) used set B alone as their reference parameter set. Here it is one of the 12, used for single-set examples and for the single-set figure `latent_dynamics.png`. With transmission off, only its natural-history rates matter.
- **Reference line:** the exact r₀·e^(−kτ) model, to separate the latent structure from v0.10's monthly step.

## Findings

**1. The two-state model can meet the progression targets, up to φ ≈ 0.01 with the cascade fixed, and up to φ ≈ 0.1 if clearance is allowed to rise within its prior.**

| Posterior weight meeting both targets (345 samples), within 1σ / 2σ | φ = 0 | 0.003 | 0.01 | 0.03 | 0.1 | 0.3 |
|---|---|---|---|---|---|---|
| Arm 1, cascade fixed | 97% / 100% | 95% / 100% | 92% / 99% | 74% / 97% | 46% / 65% | 18% / 30% |
| Arm 2, cascade freed within priors | 99% / 100% | 98% / 100% | 97% / 100% | 98% / 100% | 96% / 100% | 86% / 96% |

References:

- **tb_LAI_TPT as run** (the calibration's own saved values): 47% / 87%.
- **The exact r₀·e^(−kτ) model** at the same parameters: 26% / 68%.

The two-state model fits the targets better than tb_LAI_TPT itself because it is fitted to them alone. tb_LAI_TPT's calibration balanced them against the epidemiological targets.

- **What limits arm 1** is ProgressionRatio: remote progression adds progression after year 1. At φ = 0.1 the ratio's weighted median z is −1.3, and at 0.3 it is −3.2.
- **How arm 2 makes room:**
  - **The means:** it raises clearance of latent infection and recovery from non-infectious TB, so remote infections clear before they progress.
  - **Up to φ = 0.03,** the cascade barely moves: the weighted median of each freed parameter ÷ its sample value is 1.00, with at most 15% of the weight at a prior bound.
  - **At φ = 0.1,** the median `inf_cle` is ×1.34 the sample's value, with 30% of the weight at the bound of 1.5.
  - **At φ = 0.3,** it is ×1.65, with 66% at the bound.
  - **Set B** reaches φ = 0.1 with `inf_cle` 0.54 → 1.50 and `non_rec` 0.17 → 0.23, both at their upper bounds.
- **The calibration's own values carry Monte Carlo noise.** With 400k agents, `param00`'s ProgressionY1 as run is 0.0563 ± 0.0002, against 0.0539 from the calibration's 40k agents.

**2. For φ ≤ 0.03 in arm 1, cumulative progression is close to tb_LAI_TPT as run. The main differences are in how latent infection behaves.**

All values in this table are from the no-reinfection cohort, with weighted means across the 12 sets. Rows marked ÷ are the model ÷ tb_LAI_TPT as run, where 1.00 is a match.

| Output | tb_LAI_TPT as run | Exact (reference) | Arm 1, φ = 0.003 | Arm 1, φ = 0.01 | Arm 1, φ = 0.03 | Arm 2, φ = 0.03 | Arm 2, φ = 0.1 |
|---|---|---|---|---|---|---|---|
| B(50), progression to ASY (÷) | 0.081 | 0.92 | 1.00 | 1.01 | 1.03 | 1.00 | 1.00 |
| A(50), any progression (÷) | 0.109 | 0.87 | 0.94 | 0.95 | 0.98 | 0.96 | 0.96 |
| C(50), symptomatic (÷) | 0.044 | 0.89 | 0.97 | 0.98 | 1.00 | 0.95 | 0.94 |
| First exit from latent to NON (÷) | 6.1% | 0.78 | 0.84 | 0.85 | 0.87 | 0.85 | 0.86 |
| First exit from latent to ASY (÷) | 4.8% | 0.98 | 1.05 | 1.06 | 1.09 | 1.07 | 1.08 |
| First exit from latent to cleared (÷) | 89.2% | 1.02 | 1.01 | 1.01 | 1.00 | 1.01 | 1.01 |
| Expected time spent latent (÷) | 1.24 y | 1.02 | 1.01 | 1.01 | 1.01 | 0.87 | 0.70 |
| Still latent at 10 years (÷) | 0.16% | 1.24 | 1.17 | 1.15 | 1.10 | 0.82 | 0.44 |
| Share of years 1–10 progression from people still latent at year 1 | 12.5% | 14.7% | 27% | 30% | 35% | 31% | 32% |
| Latent progression hazard at 5 years (per year) | ≈ 0 | 0.0003 | 0.003 | 0.006 | 0.012 | 0.014 | 0.036 |

- **Curves A, B and C** (`progression_by_time_arm*.png`):
  - **The two-state curves** at φ ≤ 0.03 lie close to tb_LAI_TPT's (B within 3% by 50 years, C within 5%). A is 2–6% lower.
  - **Larger φ in arm 1** pushes B and C above tb_LAI_TPT's: B(50) ×1.12 at φ = 0.1 and ×1.33 at 0.3. Arm 2 keeps them close at every φ.
  - **The largest gap in B** between the two-state model and tb_LAI_TPT as run (3 months–20 years, by set) is 14–19% at φ ≤ 0.03. It falls at 3 months, where B is small.
- **Fewer first exits to non-infectious TB** (−13 to −16%) is mostly v0.10's monthly step, not the latent structure. The exact model, with tb_LAI_TPT's own parameters, shows −22%. v0.10 applies each month's start-of-step progression hazard for the whole month, which overstates early progression (exp 08), most visibly to NON. For set B by year 1, as run 8.0% of infections have left to NON, against 6.1% in the exact model.
- **Latent infection itself** (`latent_dynamics.png`):
  - **Hazard:** tb_LAI_TPT's latent hazard declines to about 0 within 2 years, while the two-state hazard levels off at φ·r₀.
  - **Where later progression comes from:** about 2–3 times as much of the progression in years 1–10 comes from people still latent at year 1 (27–35%, against 12.5%).
  - **How much latent infection there is:** in arm 1 it is almost unchanged (time spent latent within 1%).
  - **Arm 2's price** for larger φ is less latent infection: at φ = 0.1, time spent latent is ×0.70 and people still latent at 10 years ×0.44. That's because it gets there by raising clearance.
- **Recent vs remote:** in the two-state model, most latent infections are remote within a few months of infection (set B, arm 1: 4–18% still recent at 6 months for φ ≤ 0.03). v0.10 has no counterpart.

**3. Comparisons for individual parameter sets include each set's own miss of the targets.** The two-state model is fitted to the target means, but tb_LAI_TPT as run isn't at them. Set B as run has ProgressionY1 = 0.0563 (z ≈ +0.8), so its two-state B curve is 7–11% below tb_LAI_TPT's from year 1 on, even at φ = 0. Across the 12 weighted sets these misses largely cancel, which is why the table above uses the weighted means.

## What this means

- **The targets can be met by a genuinely two-state model within tb_LAI_TPT's priors.**
  - **With the cascade unchanged:** φ up to about 0.01, with remote progression about 1% of recent; set B's remote rates are 0.0035 / 0.0023 per year.
  - **Allowing clearance to rise within its prior:** φ up to about 0.1, but at the cost of 30% less time spent latent.
- **Cumulative progression and the A/B/C curves stay close to tb_LAI_TPT as run for φ ≤ 0.03.** This differs from exp 15, which compared against the exact model, because v0.10's monthly step had already pushed tb_LAI_TPT's progression up.
- **The underlying latent dynamics don't match, as expected.**
  - The two-state model has a lifelong reactivation floor, φ·r₀, where tb_LAI_TPT's hazard declines to about zero.
  - Its later progression comes 2–3 times as often from latent infection rather than from non-infectious TB.
  - That matters for interventions on latent infection, such as TPT. tb_LAI_TPT's TPT also acts on non-infectious TB, which limits the difference for this project.
- **Not addressed here:** burden with transmission on. φ ≤ 0.03 keeps cumulative progression within about 3% of tb_LAI_TPT as run, but the epidemiological targets still need checking with a recalibration.

## Next step

**Calibrate φ in the full recalibration,** if the two-state approach is kept (user, 2026-11-04).

- **Why the sweep isn't enough:** the progression targets alone don't identify φ, which is why it was swept here.
- **What the recalibration adds:** it fits the progression targets jointly with the epidemiological targets, which respond to cumulative progression beyond the targets' 10 years. So φ should be a calibrated parameter there.
- **A prior from this experiment:** about 0–0.03 with the cascade fixed, or up to about 0.1 if clearance may rise within its prior.
- **Until then:** φ may only be identifiable with all the targets, burden data included. So later experiments can sweep a couple of values in that range (user, 2026-11-05).

## Files

- `run_lai_tpt.py` (v0.10) → `lai_tpt.obj`, `versions_v010.json`, `run_lai_tpt.log`
- `run_fit.py` (v4) → `fits.csv`, `curves.obj`, `versions_v4.json`, `run_fit.log`
- `compare.py` (v4) → `summary.json`, `summary.obj`, `dynamics_by_set.csv`, `compare.log`, and the figures: `targets_met.png`, `progression_by_time_arm1.png`, `progression_by_time_arm2.png`, `latent_outcomes.png`, `latent_dynamics.png`, `scalars.png`
- Shared code: `v4/tbode4/latent.py` (now with the full cascade, `first_passage` and `latent_exits`); tests: `v4/tests/test_latent.py`
