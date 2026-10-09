# tbsim-ode-testing

Building an ODE version of TBsim that is reasonably consistent with the agent-based `tbsim.TB` as it is used in practice (e.g. `tb_LAI_TPT`), and documenting which ABM features can and cannot be matched by an ODE.

## Approach

Start from the simplest matched case: single-strain `tbsim.TB` vs an ODE, with constant demographic rates, mass-action mixing and no treatment. Then add ABM features one at a time, extending the ODE (`tbode.TBODE`) to match each one where possible. For each feature, record whether the ODE can match it, how, and what gap remains.

## Consistency by feature

### Step ordering in `tbsim.TB` (affects every row)

`tbsim.TB` itself differs from the ODE in how it advances time, separately from any feature. Within each step, it applies transitions in a fixed order (INFECTION → NON_INFECTIOUS → ASYMPTOMATIC → SYMPTOMATIC). Agents moving forward can carry on through several states in one step, while agents moving backward stop for the step. That cuts time infectious, which acts like a lower `beta`:
- 0.5% / 1.8% / 3.8% less effective `beta` at daily / weekly / monthly `dt`;
- that is, 5% / 16% / 33% less prevalence at the near-threshold baseline calibration.

This is the whole ABM↔ODE gap at baseline. `tbode.TBExact`, a version of `tbsim.TB` that samples each agent's end-of-step state from the exact one-step probabilities `expm(Q·dt)`, removes it (exp 03–05), at no extra run time (exp 06). Feature experiments run both `tbsim.TB` and `TBExact`, to show how much the step ordering adds as complexity grows (see Plan).

### Features

Each row is an ABM feature of the tb_LAI_TPT setup, added one at a time on top of the rows above it. Where the two ABM versions differ, the status is given for each.

Status: ⬜ not started · 🟡 in progress · ✅ matched · ⚠️ approximate only · ❌ not achievable

| Feature (ABM) | Priority | ODE approach | Status | Notes |
|---|---|---|---|---|
| Baseline: natural history, single strain, constant birth and death rates, no treatment, `ss.RandomNet(dur=0)` | — | `TBODE` (generalised `TB_ODE`), `beta` tuned to 500/100k. The random network, which gives fresh contacts each step, becomes mass action, with the per-contact ABM `beta` set by `beta_edge = -ln(1 - beta·dt/K)/dt`, where K = `n_contacts` is the mean number of contacts | ✅ `TBExact` · ⚠️ `tbsim.TB` | ABM ÷ ODE prevalence at daily / weekly / monthly `dt`: `TBExact` 1.00 / 1.00 / 0.97 ± 0.02 (40 reps); `tbsim.TB` 0.95 / 0.84 / 0.67 (step ordering, above). Transmission is exact in both (exp 01–05) |
| Age-specific fertility/mortality (WPP) | | Age-structured ODE | ⬜ | |
| No transmission from children under 10 | | Age-structured FOI | ⬜ | Requires age structure |
| Progression by time since infection (`k_asy`, `k_non`) | **High** | Latent sub-compartments by time since infection | ⬜ | |
| Reinfection of latent people (clock reset) | **High** | Move reinfected latents back to the first sub-compartment | ⬜ | Only matters with time-since-infection progression |
| Treatment (coverage, adherence, efficacy, failure, relapse) | | Treatment and relapse-pending compartments | ⬜ | ABM uses fixed delays (6-month treatment, 1.5-year relapse) |
| Time-varying `sym_dead`, efficacy, coverage | | Time-varying parameters | ⬜ | |
| Household network | Low | Effective `beta` / household-structured ODE | ⬜ | Likely ⚠️ or ❌ |

## Plan

Experiments 01–06 are done; results are in the findings log below. The baseline is consistent: with exact step sampling (`TBExact`), the ABM matches the ODE. The 40-rep rerun of exp 05 settled the remaining question, so the population-size experiment considered earlier isn't needed. `TBExact` takes the same run time as `tbsim.TB` at every `dt` (exp 06).

**Decision (user, 2026-10-09): run every feature experiment with both `tbsim.TB` and `TBExact`.** Each comparison then separates two things:
- **`TBExact` vs the ODE:** how well the ODE can match the feature, without step ordering getting in the way.
- **`tbsim.TB` vs `TBExact`:** how much step ordering adds as features are added. At baseline it is 0.95 / 0.84 / 0.67 of the ODE prevalence at daily / weekly / monthly `dt`. Each feature experiment reports whether that gap grows.

This means `TBExact` must be extended alongside the ODE for each feature that changes the natural history. Fixing tbsim upstream remains an option to revisit later (see [exp 05's SUMMARY](experiments/05_exact_transitions/SUMMARY.md)).

Next:

1. **Progression by time since infection (`k_asy`/`k_non`)**, the first high-priority feature, followed by **latent reinfection**. *Not started.*
   - **`TBExact` needs extending first.** It assumes constant rates, so it doesn't support `k` > 0 yet. With `k` > 0, each agent's exit hazard from INFECTION depends on its time since infection, so exact sampling needs that hazard integrated over the step for each agent: a transition matrix per agent, or agents grouped by time since infection. It should be re-profiled afterwards (exp 06 covered constant rates only).
   - **Time step.** Monthly `TBExact` is adequate at baseline (0.97 ± 0.02 of the ODE), but with `k` ≈ 5.85/yr the progression hazard falls by a factor of e about every 2 months. Infections should be timed within the step rather than at its start, and monthly compared directly with weekly `dt` ([exp 05](experiments/05_exact_transitions/SUMMARY.md)).
   - **The equilibrium start needs extending.** `init_from_ode` assigns each agent a state from the ODE equilibrium proportions. That's a complete description only while rates have no memory. It sets `ti_infected = 0` for all infected agents, so with `k` > 0 every latent agent would look newly infected, causing a burst of early progression (tb_LAI_TPT has `k` ≈ 5.85/yr). The ABM's initial times since infection will need drawing from the ODE's time-since-infection sub-compartments. The same applies to later features with memory: age-specific rates need an initial age distribution, and fixed-duration treatment and relapse delays need initial time in state.

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
| 4 | Which per-step flows cause the gap: transmission or natural-history transitions? | Natural history. Transmission is exact. With the step ordering, forward movers carry on within a step and backward movers stop, which shortens time infectious | [04_flux_check](experiments/04_flux_check/) |
| 5 | Does exact per-step sampling (`TBExact`) remove the gap? | Yes: prevalence 1.00 / 1.00 / 0.97 ± 0.02 of the ODE at daily / weekly / monthly `dt` (40 reps). Daily `TBExact` (1.00) vs daily `tbsim.TB` (0.95) shows even the daily gap is step ordering. Per-step flows match `expm(Q·dt)` (sum of z² 18–19 over 14 flows, vs 691–2,688 for `tbsim.TB`) | [05_exact_transitions](experiments/05_exact_transitions/) |
| 6 | How much longer does `TBExact` take to run than `tbsim.TB`? | No longer: 0.98–1.01× the run time at every `dt`, within rep-to-rep noise. Both take ~10.5 ms per step at 50k agents, so run time scales with the number of steps (daily ≈ 30× monthly). About 43% of each step is the TB module | [06_profiling](experiments/06_profiling/) |
