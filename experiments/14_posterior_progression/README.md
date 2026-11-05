# 14: Is the two-state model close to front-loaded progression across tb_LAI_TPT's posterior?

> **Superseded by [exp 15](../15_two_state_fit/README.md)** (user decision, 2026-11-04).
>
> - **Why:** this experiment fixes progression from the remote latent state at 0. That makes the two-state model effectively a single latent state with raised clearance, and we don't want to fix it at 0.
> - **Instead:** exp 15 fits a genuinely two-state model with remote progression over a range of values, and includes 0 as one point in that range.

**Track:** v4 environment (`v4/.venv`; tbsim `origin/v4` @ `dcfb522`, starsim `origin/v4` @ `80c81f53`). The calculation uses only numpy and scipy, so it doesn't depend on either package's code; versions are recorded anyway.

## Question

tb_LAI_TPT (on tbsim v0.10) makes progression from latent infection front-loaded. The rates to non-infectious and asymptomatic TB decline as r·e^(−kτ) with time since infection τ, and reinfection resets τ. tbsim v4 replaces this with two latent states: recent infections progress and become remote at rate `k_latent`, and remote ones don't progress.

For set B (row 0 of tb_LAI_TPT's posterior), the [options note](../../proposals/02_latent_progression/options.md) found:

- **Progression:** with `k_latent` refitted so that ProgressionY1 matches the exact model (5.533 vs `k` = 5.846), the two-state model's cohort progression is within about 1% of the exact model's at every time from 3 months to 50 years.
- **Reinfection:** its reinfection rule (remote → recent; a recent infection is unchanged) adds the same progression as a clock reset, to within 0.0001.

That closeness depends on r₀/k being small, where r₀ = `inf_non` + `inf_asy` (0.12 for set B). **Does it hold across the posterior?** This is the open question the note leaves for the decision on which latent structure to propose for v4.

The same calculation also answers the multiple-strains question the note left open. Two facts hold in every candidate and in v0.10:

- **One clock per person.** The progression hazard depends on the person's time since infection, not on which strains they carry.
- **No progression penalty for superinfection** (`rr_prog_super` = 1). So it isn't penalised in tb_LAI_TPT's resistance runs.

So what the two-state rule loses (no reset for recent infections, whether the new strain is the same or different) depends only on the total force of infection on latent people. That is exactly what the single-strain reinfection scenario below varies.

## Design

**Parameters.** All 345 samples in `../tb_LAI_TPT/results/hm_real_20260728/trajectory_selection/collapsed_samples_eps1.2.csv`, weighted by `sample_weight`:

- natural history: `inf_cle`, `inf_non`, `inf_asy`, `non_rec`, `non_asy`;
- `k` = `k_asy` (= `k_non` in every sample);
- π = `rr_reinfection_rec`, which tb_LAI_TPT applies to every reinfection risk ratio.

Across samples, `k` runs from 0.01 to 11.8 per year and r₀/k from 0.06 to 14.5 (weighted 10–90%: 0.08–0.21).

**Structures.** Each is a continuous-time Markov chain for one person, solved exactly with `expm_multiply`, as in [`compare_structures.py`](../../proposals/02_latent_progression/compare_structures.py):

- **A. Exact:** Poisson-mixture classes. A new infection joins class m with probability Poisson(m; r₀/k) and progresses at m·k. The number of classes is set per sample so the truncated weight is below 1e-12. A reinfected latent person is reassigned across the classes with the same weights, which is an exact clock reset.
- **B. Two states, `k_latent` = k.** Recent infections progress at r₀ and become remote at `k_latent`; remote infections don't progress. Reinfection moves remote → recent and leaves recent unchanged.
- **C. Two states, refitted:** as B, with `k_latent` chosen (root-finding) so that ProgressionY1 matches A's.

**Scenarios.**

1. **Cohort.** Everyone latent at t = 0, with no reinfection or deaths, and asymptomatic TB absorbing. This is how tb_LAI_TPT's targets are defined. Output: the share that has ever reached asymptomatic TB, B(t), on a grid every 1/48 year to 50 years.
2. **Reinfection.** The same cohort with background deaths (1/70 per year) and a constant force of infection λ.
   - People no longer infected are reinfected at π·λ.
   - Latent people are reinfected at π·λ by each structure's own rule, or not at all.
   - The difference is the progression the rule adds, by 10 and 30 years.
   - λ = 0.12 per year (set B at 500/100k) and 0.60 per year (5×), as in the options note. Each sample's own equilibrium λ isn't used: it would need tb_LAI_TPT's network and treatment.

**Outputs, per sample:**

- **Progression gaps:** the largest relative gap |B_C(t)/B_A(t) − 1| over 3 months to 50 years, and the signed gaps at 3 months and 1, 5, 10 and 50 years; the same for B.
- **Fitted values:** ProgressionY1 = B(1) and ProgressionRatio = B(1)/B(10) for each structure, as tb_LAI_TPT defines them (`cohort_progression.py`), and `k_latent`/k for C.
- **Where later progression comes from:** in the exact model, the progression between years 1 and 10, split into people still latent at year 1 and people already in non-infectious TB.
- **Size:** the number of classes A needs at truncations of 1e-6 and 1e-10.
- **Reinfection:** the progression each structure's rule adds, and C − A.

**Summaries:**

- weighted median and 10–90% interval, and the maximum over samples;
- the share of posterior weight with the largest gap above 1%, 2% and 5%;
- the gaps plotted against r₀/k.

**What would change the decision.** If C stays within about 1% wherever the posterior has appreciable weight, the result for set B generalises. If not, the samples where it fails show which parameter values need the exact classes.

## Files

- `run.py`: the calculation. Run from the repo root: `v4/.venv/bin/python experiments/14_posterior_progression/run.py`
- `samples.csv`: per-sample results; `summary.json`: the summaries; `outputs.obj`: the cohort curves behind the figure (`sc.save`)
- `posterior_progression.png`; `versions.json`
- `SUMMARY.md`: findings, written after running
