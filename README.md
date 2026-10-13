# tbsim-ode-testing

Building an ODE version of TBsim that is reasonably consistent with the agent-based `tbsim.TB` as it is used in practice (e.g. `tb_LAI_TPT`), and documenting which ABM features can and cannot be matched by an ODE.

## Approach

Start from the simplest matched case: single-strain `tbsim.TB` vs an ODE, with constant demographic rates, mass-action mixing and no treatment. Then add ABM features one at a time, extending the ODE (`tbode.TBODE`) to match each one where possible. For each feature, record whether the ODE can match it, how, and what gap remains.

## Consistency by feature

### Step ordering in `tbsim.TB` (affects every row)

`tbsim.TB` itself differs from the ODE in how it advances time, separately from any feature. Within each step, it applies transitions in a fixed order (INFECTION → NON_INFECTIOUS → ASYMPTOMATIC → SYMPTOMATIC). Agents moving forward can carry on through several states in one step, while agents moving backward stop for the step. That cuts time infectious, which acts like a lower `beta`:
- 0.5% / 1.8% / 3.8% less effective `beta` at daily / weekly / monthly `dt`;
- that is, 5% / 16% / 33% less prevalence at the near-threshold baseline calibration.

This is nearly the whole ABM↔ODE gap at baseline. The rest is transmission: each agent's per-step infection probability is `1 − e^(−x)` with a susceptible pool fixed at the start of the step, a fraction of about x/2 below the ODE's flow (exp 10). That matters near threshold at monthly `dt` (up to 6% of prevalence in set A) and is fixed by mass action with a linear probability. `tbode.TBExact`, a version of `tbsim.TB` that samples each agent's end-of-step state from the exact one-step probabilities `expm(Q·dt)`, removes it (exp 03–05), at no extra run time (exp 06). From exp 07 on, `TBExact` also dates each new infection within its step. Feature experiments run both `tbsim.TB` and `TBExact`, to show how much the step ordering adds as complexity grows (see Plan).

### Features

Each row is an ABM feature of the tb_LAI_TPT setup, added one at a time on top of the rows above it. Where the two ABM versions differ, the status is given for each.

Status: ⬜ not started · 🟡 in progress · ✅ matched · ⚠️ approximate only · ❌ not achievable

| Feature (ABM) | Priority | ODE approach | Status | Notes |
|---|---|---|---|---|
| Baseline: natural history, single strain, constant birth and death rates, no treatment, `ss.RandomNet(dur=0)` | — | `TBODE` (generalised `TB_ODE`), `beta` tuned to 500/100k. The random network, which gives fresh contacts each step, becomes mass action, with the per-contact ABM `beta` set by `beta_edge = -ln(1 - beta·dt/K)/dt`, where K = `n_contacts` is the mean number of contacts | ✅ `TBExact` · ⚠️ `tbsim.TB` | ABM ÷ ODE prevalence at daily / weekly / monthly `dt`: `TBExact` 1.00 / 1.00 / 0.97 ± 0.02 (40 reps); `tbsim.TB` 0.95 / 0.84 / 0.67 (step ordering, above). Per-step transmission is `1 − e^(−x)` (x = `beta`·w·`dt`) in both, a fraction of about x/2 below the ODE's flow: about 1% of prevalence at monthly `dt` near threshold (exp 10) |
| Age-specific fertility/mortality (WPP) | | Age-structured ODE | ⬜ | |
| No transmission from children under 10 | | Age-structured FOI | ⬜ | Requires age structure |
| Progression by time since infection (`k_asy`, `k_non`) | **High** | Poisson-mixture latent classes (exact: class m with probability Poisson(m; r₀/k), progressing at rate m·k; 7–8 classes) | ✅ `TBExact` · ❌ `tbsim.TB` | `TBExact` ÷ ODE prevalence, set B (tb_LAI_TPT natural history): 1.00 / 1.01 ± 0.01 at weekly / monthly `dt` (40 reps). `tbsim.TB` ÷ ODE: 1.18 / **1.83**, because it applies the start-of-step hazard for the whole step, overstating early progression by 10–19%. Set A (tbsim defaults, near threshold): `TBExact` 0.99 / 0.94 ± 0.02, `tbsim.TB` 1.30 / 2.64 (exp 08) |
| Reinfection of latent people (clock reset) | **High** | Reinfected latents are reassigned across the classes with the Poisson weights (exact) | ⬜ | Only matters with time-since-infection progression |
| Treatment (coverage, adherence, efficacy, failure, relapse) | | Treatment and relapse-pending compartments | ⬜ | ABM uses fixed delays (6-month treatment, 1.5-year relapse) |
| Time-varying `sym_dead`, efficacy, coverage | | Time-varying parameters | ⬜ | |
| Household network | Low | Effective `beta` / household-structured ODE | ⬜ | Likely ⚠️ or ❌ |

## Plan

Experiments 01–11 are done; results are in the findings log below.
- **The baseline is consistent:** with exact step sampling (`TBExact`), the ABM matches the ODE. The 40-rep rerun of exp 05 settled the remaining question, so the population-size experiment considered earlier isn't needed. `TBExact` takes the same run time as `tbsim.TB` at every `dt` (exp 06).
- **`TBExact` now dates new infections within the step** (`infection_timing='within'`, exp 07), which removes exp 05's INFECTION artefact. Near-threshold shortfalls of up to 6% (exps 05, 08, 09) are caused by the per-step infection probability (exp 10). One residual is still unexplained: about 3.5% at `k` = 0 with weekly `'within'` timing. Set B, the reference, matches to ±1%.

**Decision (user, 2026-10-09): run every feature experiment with both `tbsim.TB` and `TBExact`.** Each comparison then separates two things:
- **`TBExact` vs the ODE:** how well the ODE can match the feature, without step ordering getting in the way.
- **`tbsim.TB` vs `TBExact`:** how much step ordering adds as features are added. At baseline it is 0.95 / 0.84 / 0.67 of the ODE prevalence at daily / weekly / monthly `dt`. Each feature experiment reports whether that gap grows.

This means `TBExact` must be extended alongside the ODE for each feature that changes the natural history. Fixing tbsim upstream remains an option to revisit later (see [exp 05's SUMMARY](experiments/05_exact_transitions/SUMMARY.md)).

**`dt`s for feature experiments:** monthly (tb_LAI_TPT's setting) and weekly. Add daily only if monthly and weekly differ substantially.

**Parameters from exp 08 on (user decision, 2026-10-12):**
- **Set B, the reference:** tb_LAI_TPT's best-fit natural history, with `k_asy = k_non` = 5.846 and reinfection risk ratio π = 0.24 for everyone in CLEARED. The pre-1995 `sym_dead` is used.
- **Set A, a bridge to exps 01–07:** tbsim defaults + `k`.
- In both, `beta` is recalibrated to 500/100k. See [exp 08's README](experiments/08_progression_tsi/README.md).

Next:

1. **Decided (user, 2026-10-12): no changes to tbsim; `TBExact` is the corrected ABM.** At monthly `dt` with `k` = 5.846, `tbsim.TB` gives 1.8× (set B) to 2.6× (set A) the ODE prevalence, because it applies the start-of-step hazard for the whole step (exp 08). `TBExact` will accumulate the fixes as features are added.
   - tb_LAI_TPT will need recalibrating with `TBExact` to be consistent with the ODE; that's deferred.
   - `k` is not recalibrated: measured exactly, set B already fits the data-derived progression targets (within 0.25 SD), better than tb_LAI_TPT's monthly `tbsim.TB` measurement does.
   - `TBExact` should be re-profiled with latent classes (exp 06 covered constant rates only).
   - Open: whether to keep running `tbsim.TB` in every experiment or only at key points.
2. **Done (exp 11): count incidence from transitions.** `TBExact`'s `incidence_event_py` matches the ODE's incidence per prevalent case to within 0.3%. The snapshot count (`incidence_py`) misses 3.3% of episodes at monthly `dt` and 0.75% at weekly. The transition-based count is the standard for ODE↔ABM comparisons (user decision, 2026-10-13). Notifications will need the same approach once treatment is added. Deferred: whether calibration targets should exclude brief subclinical episodes. If so, define the rule by duration in both models (see exp 11's SUMMARY).
3. **Exp 12: latent reinfection (clock reset),** with `rr_reinfection_inf = rr_reinfection_non` = π. In the classes, a reinfected latent agent is reassigned with the Poisson weights, which is also exact.
4. **Decision needed: transmission in `TBExact` (exp 10).** The near-threshold shortfall is the ABM's per-step infection probability. With mass action and a linear probability (`TBExactMassAction(infection_prob='linear')`), `TBExact` matches the ODE to within 1–2% at monthly `dt` in both sets (A: 0.984 ± 0.010; B: 1.000 ± 0.006). tb_LAI_TPT's random network (`dur=0`) behaves like mass action, so this could become `TBExact`'s default transmission. The household network, if ever matched, would still need network transmission.
5. **Later features with memory need the equilibrium start extended again:** age-specific rates need an initial age distribution, and fixed-duration treatment and relapse delays need initial time in state.

## Open issues

Differences between `TBExact` and the ODE that are not fully explained. They are listed in [exp 10's SUMMARY](experiments/10_transmission_step/SUMMARY.md), with the history in the experiments cited.

| Issue | Size | Status |
|---|---|---|
| `k` = 0, weekly `dt`, `'within'` timing: prevalence low (exps 07, 10) | 0.966 ± 0.007 of the ODE (320 reps pooled), against 0.995 predicted by a deterministic version of the step (z ≈ −4) | **Unexplained.** Ruled out: the network, per-step sampling, the random-number source, population size, and the per-step infection probability. The same setting with `'start'` timing matches, as do monthly `'within'` and set A weekly. Doesn't affect set B or monthly runs |
| CLEARED slightly low (exp 08) | 0.998–0.999 of the ODE, z = −2 to −5, even where prevalence matches | Not examined; small |

## Repo layout

```
tbode/          Python package: the ODE (ode.py), matched-ABM builders and analyzers (abm.py),
                tbsim.TB with exact per-step transitions (exact.py), comparison helpers (compare.py),
                per-step flux checks against expm(Q·dt) (flux.py)
experiments/    One folder per experiment (NN_name/): README.md (question, design), run.py, outputs, SUMMARY.md
tests/          Regression checks (e.g. TBODE with all features off reproduces tbsim's TB_ODE)
```

## Setup

```bash
uv sync                                       # creates .venv with editable ../tbsim and ../starsim
uv run python experiments/01_baseline/run.py
```

## Findings log

| # | Question | Result | Experiment |
|---|---|---|---|
| 1 | How well does the simplest matched ABM (`tbsim.TB`, monthly `dt`) reproduce the ODE? | Poorly: 10–35% below the ODE in prevalence and latent infection, at every `beta` tested (9, 12, 20) | [01_baseline](experiments/01_baseline/) |
| 2 | Does a smaller ABM time step close the gap? | Largely, well above threshold: within ~3% at weekly and daily `dt` at `beta` = 12 and 20. But not entirely: near threshold even daily `dt` is 5% low (exp 03). Exp 05 later showed the whole gap is a time-step effect, namely step ordering | [02_timestep](experiments/02_timestep/) |
| 3 | With `beta` calibrated to 500/100k and the ABM started at the ODE equilibrium, does the ABM stay there, at daily, weekly and monthly `dt`? | No: it drifts down at every `dt`, to 0.95 ± 0.02 / 0.84 / 0.67 of the ODE prevalence at daily / weekly / monthly `dt`. Even daily is significantly low (z = −2.9). The calibration sits near threshold (`beta` = 9.92 vs threshold 8.59), so prevalence is ~9× as sensitive as `beta`; the gaps are 0.5% / 1.8% / 3.8% in effective `beta` | [03_baseline_calibrated](experiments/03_baseline_calibrated/) |
| 4 | Which per-step flows cause the gap: transmission or natural-history transitions? | Natural history. Transmission matches its per-step expectation, `1 − e^(−x)`; exp 10 later showed that expectation is slightly below the ODE's flow. With the step ordering, forward movers carry on within a step and backward movers stop, which shortens time infectious | [04_flux_check](experiments/04_flux_check/) |
| 5 | Does exact per-step sampling (`TBExact`) remove the gap? | Yes: prevalence 1.00 / 1.00 / 0.97 ± 0.02 of the ODE at daily / weekly / monthly `dt` (40 reps). Daily `TBExact` (1.00) vs daily `tbsim.TB` (0.95) shows even the daily gap is step ordering. Per-step flows match `expm(Q·dt)` (sum of z² 18–19 over 14 flows, vs 691–2,688 for `tbsim.TB`) | [05_exact_transitions](experiments/05_exact_transitions/) |
| 6 | How much longer does `TBExact` take to run than `tbsim.TB`? | No longer: 0.98–1.01× the run time at every `dt`, within rep-to-rep noise. Both take ~10.5 ms per step at 50k agents, so run time scales with the number of steps (daily ≈ 30× monthly). About 43% of each step is the TB module | [06_profiling](experiments/06_profiling/) |
| 7 | Does dating each new infection within its step (`TBExact`, `infection_timing='within'`) remove exp 05's INFECTION artefact? | Yes: monthly INFECTION goes from 0.89 to 0.98 ± 0.02 of the ODE, and relative to prevalence from 0.91 to 1.00 (weekly: 0.98 to 1.00). Monthly prevalence, incidence and mortality are unchanged. At weekly `dt`, `'within'` came out 0.96 ± 0.01 vs 1.00 for `'start'` (160 reps each). Diagnostics found no mechanism or bug: per-step flows are exact; a deterministic version of the step predicts 0.999; two diagnostic variants agree within noise. The gap is probably mostly chance, but up to ~3% of prevalence near threshold can't be ruled out | [07_infection_timing](experiments/07_infection_timing/) |
| 8 | With progression by time since infection (`k` = 5.846), does the ODE (Poisson-mixture latent classes) match the ABM? | `TBExact`: yes. Prevalence 1.00 / 1.01 ± 0.01 of the ODE at weekly / monthly `dt` for set B (tb_LAI_TPT natural history); set A (near threshold) 0.99 / 0.94 ± 0.02. Daily `tbsim.TB` cohorts match the exact progression curve, which confirms the latent classes independently. `tbsim.TB`: no. 1.18 / 1.83× the ODE (set B) and 1.30 / 2.64× (set A), because the start-of-step hazard overstates early progression (cohort +10–19% at monthly `dt`). The step-ordering gap reverses and grows | [08_progression_tsi](experiments/08_progression_tsi/) |
| 9 | Is the near-threshold `TBExact` shortfall (set A, monthly: 0.94 of the ODE vs 0.995 predicted by a deterministic version of the step) a finite-population effect? | No. It stays about 5% at 25k / 50k / 100k / 200k agents (0.947 / 0.924 / 0.964 / 0.941 ± ~0.01). A constant fits (χ² = 7 on 3 df); 1/N doesn't (χ² = 60). Rep-to-rep SD shrinks with N, but the mean shortfall doesn't. Cause found in exp 10 | [09_population_size](experiments/09_population_size/) |
| 10 | Is the near-threshold shortfall caused by the ABM's per-step infection probability (`1 − e^(−x)`, with a susceptible pool fixed at the start of the step) rather than the ODE's flow (x)? | Yes. Set A monthly: mass action with `1 − e^(−x)` gives 0.941 ± 0.010, exactly as predicted, and the network gives the same (0.947), so the network adds nothing. With a linear probability: 0.984 ± 0.010 (A) and 1.000 ± 0.006 (B). Exp 05's monthly 0.97 is explained (rerun: 0.980 ± 0.010 vs 0.997 predicted). Still unexplained: `k` = 0 weekly `'within'`, 0.966 ± 0.007 pooled vs 0.995 | [10_transmission_step](experiments/10_transmission_step/) |
| 11 | Is `TBExact`'s monthly incidence shortfall (about 3% below prevalence) a snapshot artefact, and does counting from transitions fix it? | Yes. Counting each agent's within-step visits to ASYMPTOMATIC/SYMPTOMATIC (augmented chain, own random stream, dynamics bit-identical) gives incidence equal to prevalence ÷ ODE in every arm. Incidence per prevalent case is 1.002–1.003 ± 0.001 of the ODE at monthly `dt` and 1.000 at weekly. The snapshot misses 3.2–3.3% of episodes at monthly `dt` and 0.75% at weekly, whatever the parameters | [11_incidence_counting](experiments/11_incidence_counting/) |
