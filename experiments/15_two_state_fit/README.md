# 15: Can a genuine two-state latent model fit the progression targets and reproduce tb_LAI_TPT's progression dynamics?

> **Superseded by [exp 16](../16_two_state_vs_lai_tpt/README.md)** (user decision, 2026-11-04).
>
> - **Why:** this experiment fits each posterior sample's own exact-model values and compares against the exact r₀·e^(−kτ) model, which isn't what tb_LAI_TPT ran: its progression targets were computed with v0.10's monthly ABM.
> - **Instead:** exp 16 fits the target means directly, sweeps φ, and compares with tb_LAI_TPT as run.

**Track:** v4 environment (`v4/.venv`; tbsim `origin/v4` @ `dcfb522`, starsim `origin/v4` @ `80c81f53`). The calculation uses only numpy and scipy (`tbode4.latent`, checked by `v4/tests/test_latent.py`), so it doesn't depend on either package's code; versions are recorded anyway.

## Question

tbsim v4 models latent infection as two states. Recent infections progress at `inf_non` / `inf_asy` and become remote at rate `k_latent`. Remote infections progress at `inf_non_remote` / `inf_asy_remote` and clear at `inf_cle`, as recent ones do.

[Exp 14](../14_posterior_progression/SUMMARY.md) matched tb_LAI_TPT's front-loaded progression r·e^(−kτ) to within ~1% with this model, but only with **no progression from the remote state**. Then remote infections never progress, and for progression they are the same as cleared ones. So that version amounts to a single latent state with clearance raised to about `inf_cle` + 0.95k, which is well outside tb_LAI_TPT's prior for clearance (0.5–1.5 per year). The user doesn't want an effectively single-state model, and wants to stay within the calibration's priors where possible.

With remote progression > 0, within priors:

- **(a)** Can the two-state model fit the progression targets (ProgressionY1 = B(1) = 0.0523, ProgressionRatio = B(1)/B(10) = 0.663)? Here B(t) is the share of a latent cohort that has ever reached asymptomatic TB. How much remote progression do the targets allow?
- **(b)** Does it reproduce tb_LAI_TPT's progression by time since infection: the curve B(t), the progression hazard among people still latent, and where late progression comes from? If (a) holds and (b) doesn't, that is an acceptable outcome (user, 2026-11-04).
- **Side effects:** fitting only the progression targets could change cumulative progression and so pull the model off tb_LAI_TPT's epidemiological targets. How large is that shift? A full recalibration is a possible later step.

## Design

**The two-state model and φ.** Rates per year. The two latent states are recent (R) and remote (M):

| From | To NON (non-infectious TB) | To ASY (asymptomatic TB) | Total progression | To remote | To cleared |
|---|---|---|---|---|---|
| Recent | `inf_non` | `inf_asy` | r₀ = `inf_non` + `inf_asy` | `k_latent` | `inf_cle` |
| Remote | `inf_non_remote` = φ · `inf_non` | `inf_asy_remote` = φ · `inf_asy` | φ · r₀ | — | `inf_cle` |

- **What φ is:** the ratio of remote to recent progression, φ = `inf_non_remote` / `inf_non` = `inf_asy_remote` / `inf_asy`.
- **One φ for both routes.** Remote infections progress at a fixed fraction φ of the recent rates, with the same NON : ASY split.
- **The ends of the range:**
  - φ = 0: remote infections never progress, so for progression they are the same as cleared (exp 14's model).
  - φ = 1: remote infections progress like recent ones, so there is no front-loading and `k_latent` has no effect.
- **The resulting hazard:** among people still latent τ years after infection, the progression hazard is h(τ) = r₀ · [q(τ) + φ · (1 − q(τ))], where q(τ) is the share of them still recent.
  - In a cohort without reinfection, q falls from 1 toward 0 at rate about `k_latent` + r₀. So h steps down from r₀ to a floor of φ · r₀, and stays there for as long as people stay latent.
  - In tb_LAI_TPT's exact model, the hazard is r₀ · e^(−kτ), which declines to 0.

**Posterior samples.** All 345 samples of `../tb_LAI_TPT/results/hm_real_20260728/trajectory_selection/collapsed_samples_eps1.2.csv`, weighted by `sample_weight`. Each sample's exact model is the reference (structure A of exp 14): Poisson-mixture classes, which reproduce r·e^(−kτ) exactly.

**Fitting.** Changing one thing, the latent structure. For each sample, the two-state model:

- **Keeps the sample's other rates:** `inf_cle` (for both latent states), `non_rec`, `non_asy`, π, and the split of progression between NON and ASY (`inf_non` : `inf_asy`).
- **Fits two parameters,** the total progression rate of recent infections r₀ = `inf_non` + `inf_asy`, and `k_latent`. They are fitted so that ProgressionY1 and ProgressionRatio equal the sample's own exact-model values.
  - So the two-state model fits the targets exactly as well as that posterior sample does.
  - Each sample's fit to the targets is also reported as log-space z-scores, z = log(model ÷ target mean) ÷ σ, with σ = 0.094 / 0.071.
- **Sets φ** (defined above) to 0, 0.003, 0.01, 0.03, 0.1 and 0.3. For each value, r₀ and `k_latent` are refitted, so the recent rates change with φ, and the remote rates follow as φ times the refitted recent rates. φ = 0 is exp 14's model, refitted on both parameters.
- **Uses the agreed priors** (user, 2026-11-04):
  - **Recent rates:** within tb_LAI_TPT's widened ranges, `inf_non` 0.04–0.45 and `inf_asy` 0.01–0.60.
  - **`k_latent`:** 0–12, like `k`.
  - **Remote rates:** above 0, at most the original constant-rate upper bounds (`inf_non` ≤ 0.23, `inf_asy` ≤ 0.20), and at most the recent rates. The caps can't bind on this grid.
  - **`inf_cle`:** stays at the sample's value, within its prior of 0.5–1.5. Freeing it could be a later variant.
- **Solving:** bounded least squares on the two log-ratios. A fit is **feasible** if both are within 1e-4, far below the target σs. Otherwise the remaining misfit is reported in σ units.

**Outputs, per sample and φ, against the sample's exact model:**

- **(a) Fit:**
  - feasibility, the misfit, and the fitted r₀, `k_latent`, and recent and remote rates;
  - for each sample, the largest φ on the grid that is feasible.
- **(b) Progression dynamics:**
  - B(t) on a grid every 1/48 year to 50 years: the largest |two-state ÷ exact − 1| over 3 months–50 years, and the gaps at 2, 5, 10 and 50 years;
  - the progression hazard among people still latent at τ = 0, 0.5, 1, 2, 5 and 10 years, against r₀·e^(−kτ);
  - the share of progression between years 1 and 10 that comes from people still latent at year 1, rather than already in non-infectious TB.
- **Side effects (cumulative progression):**
  - lifetime progression, B(50);
  - expected time spent latent (latent stock per infection);
  - **progression with reinfection:** a cohort of new infections followed for 30 years at a constant force of infection λ (the infection rate per year of a fully susceptible person), with background deaths at 1/70 per year.
    - People no longer infected (cleared, or recovered from non-infectious TB) are reinfected at π·λ.
    - Latent people are also reinfected at π·λ: in the exact model that resets the clock to τ = 0; in the two-state model, remote goes back to recent and recent is unchanged.
    - λ = 0.12/y is set B's force of infection at its 500/100k equilibrium (`beta` × (`trans_asymp` × asymptomatic prevalence + symptomatic prevalence), exp 13's calibration); 0.60/y (5×) is a stress test. λ isn't each sample's own equilibrium, which would need tb_LAI_TPT's full model.
    - Output: the share ever reaching asymptomatic TB by 30 years.
  - all as ratios to the exact model.

**Summaries, for each φ:**

- the posterior weight with a feasible fit;
- over feasible samples, the weighted median and 10–90% interval of each output, and the maximum.

**Figure:**

- feasibility and misfit against φ;
- row 0's (set B's) B(t) at each φ, against the exact model and against the Ferebee 1970 and Sutherland 1968 data behind the targets;
- row 0's latent hazard by τ;
- the side effects against φ.

**What would answer the question:**

- **(a):** the posterior weight with feasible fits, and how large φ can be. If the targets force φ to about 0, the data prefer an effectively single-state model within these priors.
- **(b):** the gaps in B(t) and in the latent hazard. The hazard is where a two-state model is expected to differ most: a step down from recent to remote, rather than a smooth decline.
- **Side effects:** lifetime progression and progression with reinfection. A shift of more than ~1–2% would be worth a recalibration, since prevalence amplifies changes in progression.

## Files

- `run.py`: the calculation. Run from the repo root: `v4/.venv/bin/python experiments/15_two_state_fit/run.py`; add `plot` to redraw the figure from `outputs.obj`
- `fits.csv`: per sample and φ; `summary.json`: the summaries; `outputs.obj`: everything behind the figure (`sc.save`)
- `two_state_fit.png`; `versions.json`; `run.log`
- `SUMMARY.md`: findings, written after running
