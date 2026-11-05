# 16: Can a two-state latent model meet tb_LAI_TPT's progression targets, and how do its dynamics compare with tb_LAI_TPT as run?

**Track:** both environments, with results passed between them as saved files.

- **v4 environment** (`v4/.venv`; tbsim `origin/v4` @ `dcfb522`, starsim `origin/v4` @ `80c81f53`): the two-state fits and the exact reference curves. These use `tbode4.latent`, a deterministic Markov-chain calculation of the same model v4 declares.
- **v0.10 environment** (repo root `.venv`; tbsim v0.10.3 @ `b92d7e8`, starsim v3.6.1 @ `3521adb4`): tb_LAI_TPT as run. This is an explicit exception to the rule of not rerunning v0.10 arms (user, 2026-11-04). The pinned environment reproduces tb_LAI_TPT's calibration values exactly: for `param00`, ProgressionY1 = 0.0539 and ProgressionRatio = 0.63825, from `cohort_progression.py` at its default settings.

## Question

tbsim v4 models latent infection as two states, recent and remote. tb_LAI_TPT (on v0.10) instead has a progression hazard that declines with time since infection, r₀·e^(−kτ). We want a **genuinely two-state** model: remote infections also progress, at a fraction φ of the recent rates. Without that, the remote state is the same as cleared for progression, and the model is effectively single-state (exp 14). Within tb_LAI_TPT's priors:

1. **Can it meet the progression targets?** These are ProgressionY1 = B(1) = 0.0523 and ProgressionRatio = B(1)/B(10) = 0.663, where B(t) is the share of a latent cohort that has ever reached asymptomatic TB. For which values of φ?
2. **What underlying dynamics result, and do they match tb_LAI_TPT as run?** The reference is v0.10's `tbsim.TB` at monthly `dt`, the model tb_LAI_TPT calibrated and runs.

**Why this replaces exp 15.** Exp 15 compared against the exact r₀·e^(−kτ) model. But tb_LAI_TPT computed its progression targets with v0.10's monthly ABM, which overstates early progression by 10–19% (exp 08), so the exact model isn't what tb_LAI_TPT ran. Here:

- the two-state model is fitted to the target means;
- its dynamics are compared with tb_LAI_TPT as run;
- the exact model is shown only as a reference line, to separate the effect of the latent structure from the effect of v0.10's monthly step.

**Transmission is off throughout,** in a cohort with no reinfection. So the comparison is of natural history only, not of overall burden, which will need a recalibration later (see exp 15).

## Design

### The two-state model and φ

Rates per year. The two latent states are recent (R) and remote (M):

| From | To NON (non-infectious TB) | To ASY (asymptomatic TB) | Total progression | To remote | To cleared |
|---|---|---|---|---|---|
| Recent | `inf_non` | `inf_asy` | r₀ = `inf_non` + `inf_asy` | `k_latent` | `inf_cle` |
| Remote | φ · `inf_non` | φ · `inf_asy` | φ · r₀ | — | `inf_cle` |

- **What φ is:** φ = `inf_non_remote` / `inf_non` = `inf_asy_remote` / `inf_asy`. Remote infections progress at a fixed fraction of the recent rates, with the same NON : ASY split.
- **The ends of the range:** φ = 0 means remote infections never progress (equivalent to cleared, for progression); φ = 1 means no front-loading.
- **The resulting hazard:** among people still latent τ years after infection, h(τ) = r₀ · [q(τ) + φ · (1 − q(τ))], where q(τ) is the share still recent. It steps down from r₀ to a floor of φ · r₀.
- **Active TB:** the cascade after progression is the same in all models: NON → recovered (`non_rec`) or → ASY (`non_asy`); ASY → NON (`asy_non`) or → SYM (`asy_sym`); SYM → ASY (`sym_asy`). TB deaths are off (`sym_dead` = 0), as in tb_LAI_TPT's cohort.

**φ is swept, not calibrated** (user, 2026-11-04). With the cascade fixed, the two-state model has three free parameters (r₀, `k_latent`, φ) for two targets, so φ isn't identifiable from them. The sweep shows which values meet the targets, and what dynamics each implies. Grid: φ = 0, 0.003, 0.01, 0.03, 0.1, 0.3.

### Parameter sets

- **Question 1 (targets):** all 345 samples of `../tb_LAI_TPT/results/hm_real_20260728/trajectory_selection/collapsed_samples_eps1.2.csv`, weighted by `sample_weight`.
- **Question 2 (dynamics):** the 12 parameter sets tb_LAI_TPT's face-validity checks use (`param00`–`param11` = the first 12 rows of that file).
  - **What they are:** the file is sorted by log-likelihood, best first, so these are the 12 most likely samples, not a representative draw from the posterior.
  - **Their weights:** together they hold 34% of the posterior's weight. `sample_weight` (set by how often trajectory selection drew each sample, not by rank) is renormalised over the 12, and the top two then hold 45% of it (25% and 20%).
- **Set B** is `param00`, row 0 and the best fit. Earlier experiments (08–13) used it alone as their reference set. Here it is one of the 12, used for single-set examples and figures. With transmission off, only its natural-history rates matter, not its `beta` or `sym_dead`. They are weighted by `sample_weight`, renormalised over the 12, as tb_LAI_TPT does.

### Models

1. **tb_LAI_TPT as run** (v0.10 environment).
   - **Settings:** `cohort_progression.py`'s cohort, run in a script of ours that builds the sim the same way: v0.10 `tbsim.TB`, monthly `dt`, `beta` = 0, everyone in INFECTION at the start, no demographics, `sym_dead` = 0, 50 years.
   - **Size:** 100k agents × 4 seeds per parameter set, against tb_LAI_TPT's 20k × 2. The SE of B(1) is then about 0.7%.
   - **Check:** at tb_LAI_TPT's settings (20k agents, seeds 0 and 1), our analyzer's B(1) and B(10) must equal `cohort_progression.FirstPassage`'s.
2. **Two-state, fitted to the targets** (v4 environment). For each sample and φ, fitted by bounded least squares on the two log-space z-scores, z = log(model ÷ target mean) ÷ σ (σ = 0.094 for Y1 and 0.071 for the ratio, as in tb_LAI_TPT's calibration). Two arms:
   - **Arm 1 (primary), cascade fixed.** It keeps the sample's `inf_cle`, `non_rec`, `non_asy`, `asy_non`, `asy_sym`, `sym_asy` and NON : ASY split, and fits r₀ and `k_latent`.
   - **Arm 2, cascade freed within the priors.**
     - **Also fitted:** `inf_cle` (0.5–1.5), `non_rec` (0.14–0.23) and `non_asy` (0.21–0.28).
     - **Tie-break:** with more parameters than targets, a small penalty on each one's log-distance from the sample's value, scaled by its prior width, picks the solution that moves the cascade least. Its weight is 0.01 of a target σ.
     - **Why:** exp 15 found that the route through non-infectious TB already supplies the progression after year 1, which leaves no room for remote progression. This arm tests whether staying within the priors allows more.
   - **Priors, in both arms** (user, 2026-11-04):
     - recent `inf_non` 0.04–0.45 and `inf_asy` 0.01–0.60;
     - `k_latent` 0–12;
     - remote rates at most 0.23 / 0.20 and at most the recent rates (these caps can't bind on this grid).
   - **Check:** for `param00` at one φ, v4's own `tbsim.TB(solver='ode')` cohort must match the chain's state occupancies.
3. **Exact r₀·e^(−kτ), reference line only** (v4 environment). tb_LAI_TPT's hazard computed exactly, with the Poisson-mixture classes, at each parameter set's calibrated values. The gap between it and model 1 is v0.10's monthly-step bias; the gap between it and model 2 is the latent structure.

### Outputs

All from the no-reinfection cohort, for each of the 12 parameter sets, each model, and each φ and arm:

- **Progression by time since infection,** tb_LAI_TPT's `progression_by_time.pdf` for a cohort:
  - three cumulative curves to 20 years: A (any progression: first reaching NON, ASY or SYM), B (at least asymptomatic) and C (symptomatic);
  - the Ferebee 1970 and Sutherland 1968 data overlaid on B.
  - tb_LAI_TPT's figure is in-population, with reinfection and treatment; here there is neither.
- **Latent-infection outcomes,** the bottom-right panel of `tb_state_outcomes.pdf`:
  - how each infection first leaves latent infection (cleared, to NON, to ASY), over 50 years;
  - the share still latent at 10 and 25 years.
  - For the two-state model, "latent" is recent plus remote: moving from recent to remote isn't an exit.
- **Latent prevalence:** the share of the cohort still latent by time since infection. For the two-state model, it is split into recent and remote.
  - In a cohort without reinfection, everyone still latent at time t was infected t years ago. So recent vs remote has no time-based analogue in v0.10, and is shown for the two-state model only.
- **Cumulative progression by year 50:** A(50), B(50) and C(50).
- **Expected time spent latent** per infection, over 50 years.
- **From question 1, all 345 samples:**
  - each target's best-fit z, and the fitted parameters (and arm 2's cascade values);
  - the posterior weight meeting the targets (both |z| ≤ 1, and ≤ 2) at each φ;
  - the largest φ that meets them.
  - **Reference:** the same shares for tb_LAI_TPT as run, from the calibration's own saved values (`posterior_samples.csv`, `mod__Progression*`).
- **Progression dynamics:**
  - the progression hazard among people still latent by τ (in the ABM, transitions out of INFECTION to NON or ASY in each month ÷ people in INFECTION);
  - the share of progression in years 1–10 that comes from people still latent at year 1, rather than already in NON;
  - the largest gap in curve B between the two-state model and tb_LAI_TPT as run, over 3 months–20 years.

**Summaries:**

- **12 parameter sets:** weighted mean and weighted 25–75% range of each output across the sets (as tb_LAI_TPT summarises across parameter sets), plus each set's values. For the ABM, the SE across seeds.
- **Question 1:** weighted shares over the 345 samples.

**Figures:**

- the posterior weight meeting the targets against φ, for both arms, with tb_LAI_TPT as run for reference;
- progression by time since infection (A, B, C, with the data), one panel per φ, each showing tb_LAI_TPT as run, the two-state model and the exact reference; one figure per arm;
- latent-infection outcomes as bars, by model and φ;
- latent prevalence (with the recent/remote split) and the latent hazard by τ, for `param00`;
- the scalar outputs against φ: B(50), time spent latent, and the share of later progression from latent infection.

### What would answer the question

- **Question 1:** the posterior weight meeting the targets, for each φ and arm, against tb_LAI_TPT as run.
- **Question 2:** for each φ that meets the targets, how far the two-state model's curves A, B and C, latent outcomes, latent prevalence and the other outputs are from tb_LAI_TPT as run. The exact line separates the part due to the structure from the part due to v0.10's monthly step. A mismatch in the dynamics is acceptable if the targets are met (user, 2026-11-04); the point is to describe it.

## Files

- `run_lai_tpt.py`: tb_LAI_TPT as run, in the v0.10 environment: `.venv/bin/python experiments/16_two_state_vs_lai_tpt/run_lai_tpt.py`. Writes `lai_tpt.obj`, `versions_v010.json`
- `run_fit.py`: the fits, the two-state and exact curves, and v4's ODE check, in the v4 environment: `v4/.venv/bin/python experiments/16_two_state_vs_lai_tpt/run_fit.py`. Writes `fits.csv`, `curves.obj`, `versions_v4.json`
- `compare.py`: summaries and figures from the saved files (v4 environment); add `plot` to redraw the figures only. Writes `summary.json`, `*.png`
- `run*.log`; `SUMMARY.md`: findings, written after running
