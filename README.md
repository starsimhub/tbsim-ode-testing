# tbsim-ode-testing

Building an ODE version of TBsim that is reasonably consistent with the agent-based `tbsim.TB` as it is used in practice (e.g. `tb_LAI_TPT`), and documenting which ABM features can and cannot be matched by an ODE.

**How this README is organised:**

- **[Consistency by feature](#consistency-by-feature):** where each ABM feature stands.
- **[Plan](#plan):** what to do next, in order, and which row of the table each step moves forward.
- **[Findings log](#findings-log):** one row per experiment, the record of what each found.
- **Elsewhere:**
  - modelling conventions decided along the way, written to carry into a new repo or tb_LAI_TPT: [CONVENTIONS.md](CONVENTIONS.md);
  - process decisions and established facts: [CLAUDE.md](CLAUDE.md);
  - each experiment's design and results: its folder under `experiments/`.

## Approach

Start from the simplest matched case: single-strain `tbsim.TB` vs an ODE, with constant demographic rates, mass-action mixing and no treatment. Then add ABM features one at a time, extending the ODE to match each one where possible. For each feature, record whether the ODE can match it, how, and what gap remains.

## Versions and tracks

| Track | tbsim | starsim | Environment and code | Experiments | Status |
|---|---|---|---|---|---|
| v0.10 | v0.10.3 @ `b92d7e8` | v3.6.1 @ `3521adb4` | repo root: `.venv`, `tbode/` | 01–11 | Frozen (2026-11-03): rerunnable, not extended |
| v4 | `origin/v4` @ `dcfb522` (unreleased) | `origin/v4` @ `80c81f53` (4.0.0rc1 + 21) | `v4/`: `v4/.venv`, `v4/tbode4/` | 12 on | Active |

- **Pinned commits:** each environment installs tbsim and starsim from exact git commits (`pyproject.toml`, `uv.lock`), not the editable `../tbsim` and `../starsim` checkouts that tb_LAI_TPT shares. Bump the v4 pins on purpose, then rerun `v4/tests/test_anchor.py`.
- **Recording versions:** each experiment's README has a "Track" line, and each run writes `versions*.json` beside its outputs.
- **Between tracks:** the environments can't import each other, so results pass between them as saved files. New experiments use v0.10 only as the "tb_LAI_TPT as run" reference (v0.10's `tbsim.TB` at monthly `dt`), to see how v4 changes the original ABM's behaviour, not to test v0.10 ABM↔ODE consistency (user, 2026-11-05).
- **Cliff's v4 redo** (2026-10-27; scripts at the top of `v4/`, write-up [`v4/v4_findings.qmd`](v4/v4_findings.qmd)) reran the key findings of exps 01–11 under v4. Exp 13 reproduced it exactly at our pinned versions, with more reps, so the v4 numbers here are exp 13's.

**What v4 changes.** `tbsim.TB` is declared as states and transitions, so one model runs as an ODE (`solver='ode'`) or as agents (`solver='abm'`). Its ODE equals our v0.10 ODE, `TBODE`, to 1e-5 (exp 12).

- **No step-ordering gap:** its agent solver draws each state's exits jointly from the competing rates, which removes v0.10's ABM↔ODE gaps without anything like `TBExact` (see the table).
- **Two timing effects remain** (exp 13). Neither needs action unless a calibration target is a snapshot of latent infection, or the dynamics change fast.
  - **Latent snapshots:** end-of-step snapshots of latent infection read high at monthly `dt` (+9% in set A, +3–4% in B; about half that at weekly, gone at daily). The mean over a step is right.
  - **Growth lag:** growth from a seed lags the ODE by 0.1–0.3 years.

**Where v4's model differs from the v0.10 ABM, and what we're doing about it:**

| Difference | v0.10 / tb_LAI_TPT | v4 | Status |
|---|---|---|---|
| Front-loaded progression | Hazard r·e^(−kτ) (`k_asy`, `k_non`) | Two latent states, recent and remote (`k_latent`) | **Tentative: keep v4's two states, with remote progression = φ × recent** (user, 2026-11-05; [options note](proposals/02_latent_progression/options.md), [CONVENTIONS.md](CONVENTIONS.md) 2) |
| Latent reinfection | Resets the clock | Remote → recent; recent unchanged; non-infectious TB not reinfected | Set every relative risk to π (CONVENTIONS.md 3). In a cohort, the two-state rule adds 83–93% as much progression as tb_LAI_TPT's model solved exactly (continuous time), and 75–84% as much as tb_LAI_TPT as run (v0.10, monthly `dt`): about 0.5% of cumulative progression at set B's force of infection (exp 17) |
| Treatment | Fixed 6-month course; adherence × efficacy; failure back to the previous state | Completion is a rate; competing completion and failure; failure → symptomatic. Relapse still a fixed 1.5 years | A treatment proposal, to do |
| Per-contact `beta` on a network | Matched `beta_edge` (ours) | `beta/K`: 4% fewer infections at monthly `dt`, about 25% less prevalence near threshold. Since starsim `5ba36ebb` (newer than our pins), the ODE of a network applies the same correction: calibrated there, it gives exactly our matched `beta_edge` (exp 18) | **Decided: matched `beta_edge`, in our own code only** (user, 2026-11-05; [proposal 01](proposals/01_network_beta/proposal.md), CONVENTIONS.md 1). Open: switch to the network ODE. Neither handles reinfection at a reduced risk: on a network, set B's agents run 1.20× the ODE's prevalence at monthly `dt` (exp 18) |

## Consistency by feature

Each row is an ABM feature of the tb_LAI_TPT setup, added one at a time on top of the rows above it. Status is for each ABM against its ODE: v0.10's `tbsim.TB` and our `TBExact` (against `TBODE`), and v4's agents (against v4's ODE). Values are ABM ÷ ODE prevalence (1.00 = match), mean ± SE across reps; details and rep counts are in the experiments cited.

Status: ⬜ not started · 🟡 in progress · ✅ matched · ⚠️ approximate only · ❌ not achievable

| Feature (ABM) | Priority | ODE approach | v0.10 status | v4 status | Notes |
|---|---|---|---|---|---|
| Baseline: natural history, single strain, constant birth and death rates, no treatment, `ss.RandomNet(dur=0)` | — | `TBODE` (generalised `TB_ODE`); v4's `tbsim.TB(solver='ode')` equals it. `beta` tuned to 500/100k. The random network becomes mass action, with the matched per-contact `beta_edge` | ✅ `TBExact` · ⚠️ `tbsim.TB` | ✅ | **v0.10** (daily / weekly / monthly): `TBExact` 1.00 / 1.00 / 0.97 ± 0.02; `tbsim.TB` 0.95 / 0.84 / 0.67, from step ordering (exps 03–05, 10). **v4:** set A 0.99 / 1.01 / 1.000 ± 0.009, set B 0.99 (weekly, monthly); on the network, matched `beta_edge` 1.03 ± 0.02, v4's `beta/K` 0.74 (exps 12–13) |
| Age-specific fertility/mortality (WPP) | | Age-structured ODE | ⬜ | ⬜ | v4 has age stratification with aging (`sim.stratify`) and country demography (`ss.people`) for both solvers |
| No transmission from children under 10 | | Age-structured force of infection | ⬜ | ⬜ | Needs age structure; in v4, `ss.modify(..., where=age < 10)` |
| Progression by time since infection | **High** | v0.10: Poisson-mixture latent classes (exact r·e^(−kτ)). v4: two states with remote progression = φ × recent, chosen tentatively | ✅ `TBExact` · ❌ `tbsim.TB` | ✅ | **v0.10**, set B (weekly / monthly): `TBExact` 1.00 / 1.01 ± 0.01; `tbsim.TB` 1.18 / **1.83**, because it applies the start-of-step hazard for the whole step (exp 08). **v4:** classes and `k_latent` both match across solvers, 0.99–1.01 (exp 13). **The two-state structure** meets the progression targets within tb_LAI_TPT's priors, but its latent hazard floors at φ·r₀ (exp 16). For φ ≤ 0.03 its cumulative progression is within 2–6% of tb_LAI_TPT's curve B as tb_LAI_TPT measures it (end-of-step snapshots; exp 16). Counted from transitions, it's 4–15% below tb_LAI_TPT, mostly from v0.10's step bias (exp 17) |
| Reinfection of latent people | **High** | v0.10: reassign across the classes (exact clock reset). v4: remote → recent | ⬜ | ✅ | Only matters with time-since-infection progression. **v4**, set B, two-state model at φ = 0.003–0.03: agents ÷ ODE prevalence 0.991 ± 0.006, pooled over monthly and weekly `dt`, the same as without latent reinfection; 1.004 ± 0.004 pooled over 6,250–200,000 agents. **Against tb_LAI_TPT's clock reset,** in a cohort at λ = 0.10 per year: latent reinfection adds 83–93% as much progression by 30 years as tb_LAI_TPT's model solved exactly (continuous time), and 75–84% as much as tb_LAI_TPT as run (v0.10, monthly `dt`, which overstates progression after each infection). The difference is about 0.5% of cumulative progression. Progression after each (re)infection comes later (exp 17) |
| Treatment (coverage, adherence, efficacy, failure, relapse) | | Treatment and relapse-pending compartments | ⬜ | ⬜ | v0.10 uses fixed delays; v4 a completion rate (see "Where v4's model differs") |
| Time-varying `sym_dead`, efficacy, coverage | | Time-varying parameters | ⬜ | ⬜ | v4: `ss.timeseries` with `ss.modify`, for both solvers |
| Household network | Low | Effective `beta` / household-structured ODE | ⬜ | ⬜ | Likely ⚠️ or ❌ |
| Drug resistance (tb_LAI_TPT's 2 drugs, RIF and BDQ) | Low | v0.10: `TwoStrainODE`. v4: `tbsim.TB(drugs=...)`, strains as a dimension of the states, for both solvers | ⬜ | ⬜ | Superinfection of latent infections is latent reinfection, so the latent structure must carry over; v4's two states do (`check_strains.py` in the options note's folder). Aim: v0.10 `TBResistant`'s behaviour where possible (to revisit). v4 has no strain counts |

**Step ordering in v0.10's `tbsim.TB`** (affects every v0.10 row). Within each step it applies transitions in a fixed order (INFECTION → NON → ASY → SYM). Agents moving forward can carry on through several states, while agents moving backward stop for the step. That cuts time infectious: 0.5% / 1.8% / 3.8% less effective `beta` at daily / weekly / monthly `dt`, or 5% / 16% / 33% less prevalence near threshold.

- **`tbode.TBExact`** samples each agent's end-of-step state from the exact one-step probabilities `expm(Q·dt)`, which removes this gap at no extra run time (exps 03–06).
- **The rest of the gap is transmission:** a per-step infection probability `1 − e^(−x)` with the susceptible pool fixed at the start of the step, about x/2 below the ODE's flow (exp 10). Mass action with a linear probability fixes it.
- **v4 doesn't need `TBExact`.** Its agent solver draws each state's exits jointly from the competing rates, so the order transitions are declared in doesn't matter and the step-ordering gap is gone. v4's own agents match its ODE at daily, weekly and monthly `dt` (exp 13; baseline row above). `TBExact` remains only as v0.10's corrected ABM.

## Plan

### Next, in order

1. **Treatment, next** (user, 2026-11-06; the treatment row).
   - **First:** decide a treatment proposal, not yet written (a fixed course with outcome probabilities, and failure back to the previous state, against v4's completion rate).
   - **Then:** add it to the ODE and agents, and recalibrate with mortality ÷ incidence brought back. Notifications will need a transition-based count, like incidence.
2. **The remaining rows, in table order:** age structure and demography, child transmission, time-varying parameters, then the low-priority household network and drug resistance.
   - **Starting runs:** starting agents at an equilibrium is simple while rates have no memory. Age structure needs an initial age distribution, and fixed delays (treatment, relapse) need an initial time in state.
   - **Drug resistance:** v4's resistance agents against its ODE in our matched setup, and a check of what tb_LAI_TPT's `TBResistant` setup uses that v4 doesn't have.

### Open decisions

- **Latent structure:** option B is tentative. What's open, and what would reopen it, is in the [options note](proposals/02_latent_progression/options.md#decision-tentative): φ (swept for now, calibrated with all the targets eventually), face validity of B's underlying behaviour against tb_LAI_TPT's, and reinfection of recent infections. Latent reinfection gives no reason to reopen B: in a cohort, its rule changes cumulative progression by about 0.5% at set B's force of infection (exp 17).
- **Progression targets, how curve B is counted:** tb_LAI_TPT's `FirstPassage` reads curve B (ever reached asymptomatic TB) from end-of-step snapshots. At monthly `dt` that misses about 11% of entries, people who reach asymptomatic TB and leave it within the step. Counted from transitions, tb_LAI_TPT's ProgressionY1 is 0.060, not 0.054 (exp 17). New calibrations should probably count from transitions, as for incidence (CONVENTIONS.md 5); to decide.
- **Upstream:** whether to propose anything for latent progression (proposal 02). Network `beta` stays in our code (proposal 01).
- **Network `beta` after Starsim's fix ([exp 18](experiments/18_network_ode/SUMMARY.md)):**
  - **What the fix does:** Starsim's ODE of a network now gives the same model as our matched `beta`, so our override is no longer needed. The convention could become "calibrate and run the ODE on the network".
  - **What neither handles:** on a network, both under-transmit reinfection at a reduced risk relative to agents. Set B's agents run 1.20× the ODE's prevalence at monthly `dt`, 1.05× at weekly.
  - **Decided (user, 2026-11-06):** the ODE recalibration uses homogeneous mixing, where agents match the ODE (exps 13–17) and no `beta` conversion is needed.
  - **To send:** a per-transition factor in Starsim ([proposal 03](proposals/03_network_reinfection/proposal.md), draft for review).
  - **Also open:** whether to move the v4 pins to these commits (the anchor test passes).
- **Treatment:** a proposal, not yet written, before the treatment row.

### Parameters

- **Set A:** tbsim defaults + `k` = 5.846, near threshold.
- **Set B:** tb_LAI_TPT's best-fit natural history (row 0 of its posterior) + `k` = 5.846, the reference.
- **`beta`:** recalibrated to 500/100k in each.
- **tb_LAI_TPT** stays on v0.10 for now (user, 2026-11-03). New projects will probably use v4, starting from tb_LAI_TPT's parameters, which need recalibrating for the v4 structure: its calibration partly compensates for v0.10's monthly-step biases.

### v0.10 track (frozen 2026-11-03)

Exps 01–11 are complete. `TBExact` matches the ODE at baseline and with time-since-infection progression; v0.10's `tbsim.TB` doesn't. The decision not to change tbsim (2026-10-12) stands for v0.10. Open at the freeze:

- **Infection timing in `TBExact`** (`'start'` vs `'within'`): v4 has its own timing (end of step), whose effects are the timing effects above.
- **Checking that latent classes don't change other infection dynamics** (time in INFECTION, time to progression, by outcome): applies to whichever latent structure is used.
- **Recalibrating tb_LAI_TPT with `TBExact`:** deferred, probably superseded by v4.

**Unexplained differences between `TBExact` and the ODE,** not pursued unless they reappear in v4 (details in [exp 10's SUMMARY](experiments/10_transmission_step/SUMMARY.md)):

| Issue | Size | Status |
|---|---|---|
| `k` = 0, weekly `dt`, `'within'` timing: prevalence low (exps 07, 10) | 0.966 ± 0.007 of the ODE (320 reps), against 0.995 predicted by a deterministic version of the step | Unexplained. Ruled out: the network, per-step sampling, the random-number source, the per-step infection probability, and `'within'`'s generation-interval lag. Doesn't affect set B or monthly runs |
| CLEARED slightly low (exp 08; also v4, exps 13 and 17) | 0.998–0.999 of the ODE at 50k agents, z = −2 to −5 | Not explained by the step structure. In v4 the shortfall shrinks with population size (0.6% at 6,250 agents, 0.1% at 200,000; exp 17), so it's most likely a finite-population effect. It doesn't affect prevalence or mortality, so it isn't pursued further (user, 2026-11-06) |

## Repo layout

```
CONVENTIONS.md  The modelling conventions decided here, to carry into a new repo or tb_LAI_TPT
tbode/          v0.10 track: the ODE (ode.py), matched-ABM builders and analyzers (abm.py), tbsim.TB with exact
                per-step transitions (exact.py), comparison helpers (compare.py), per-step flux checks (flux.py), versions.py
tests/          v0.10 regression checks (e.g. TBODE with all features off reproduces tbsim's TB_ODE)
v4/             v4 track: its own pyproject.toml and .venv. At the top: Cliff's v4 redo (common.py, run_*.py, analyze.py,
                v4_findings.qmd, results/summary.md)
  tbode4/       matched v4 sims as ODE or agents, and the network-beta mapping both ways (sims.py); latent-structure
                calculations (latent.py); a standalone copy of the v4 ODE with first-entry incidence, for calibration
                (fastode.py); versions.py
  tests/        test_anchor.py: the v4 ODE reproduces TB_ODE and TBODE (rerun after bumping the v4 pins);
                test_fastode.py: fastode reproduces the v4 ODE (rerun after bumping the v4 pins);
                test_latent.py: the latent-structure calculations
proposals/      Proposed changes to tbsim and starsim, and the decisions on them (some kept in our own code instead),
                numbered like experiments (NN_name/), with the scripts and outputs behind them
experiments/    One folder per experiment (NN_name/), numbered across both tracks: README.md (question, design, track),
                run script(s), outputs, versions*.json, SUMMARY.md
```

## Setup

```bash
uv sync                                         # v0.10 track: .venv, pinned tbsim v0.10.3 / starsim v3.6.1
uv run python experiments/01_baseline/run.py

cd v4 && uv sync                                # v4 track: v4/.venv, pinned tbsim and starsim v4
v4/.venv/bin/python -c "import starsim as ss; ss.compile()"   # once: compiles the kernels (~2 minutes)
v4/.venv/bin/python experiments/12_v4_setup/run.py            # from the repo root
```

## Findings log

| # | Track | Question | Result | Experiment |
|---|---|---|---|---|
| 1 | v0.10 | How well does the simplest matched ABM (`tbsim.TB`, monthly `dt`) reproduce the ODE? | Poorly: 10–35% below the ODE in prevalence and latent infection, at every `beta` tested (9, 12, 20) | [01_baseline](experiments/01_baseline/) |
| 2 | v0.10 | Does a smaller ABM time step close the gap? | Largely, well above threshold: within ~3% at weekly and daily `dt` at `beta` = 12 and 20. But near threshold even daily `dt` is 5% low (exp 03). Exp 05 later showed the whole gap is step ordering | [02_timestep](experiments/02_timestep/) |
| 3 | v0.10 | With `beta` calibrated to 500/100k and the ABM started at the ODE equilibrium, does the ABM stay there? | No: it drifts down to 0.95 ± 0.02 / 0.84 / 0.67 of the ODE prevalence at daily / weekly / monthly `dt`. The calibration sits near threshold (`beta` = 9.92 vs 8.59), so prevalence is ~9× as sensitive as `beta`; the gaps are 0.5% / 1.8% / 3.8% in effective `beta` | [03_baseline_calibrated](experiments/03_baseline_calibrated/) |
| 4 | v0.10 | Which per-step flows cause the gap: transmission or natural history? | Natural history. With the step ordering, forward movers carry on within a step and backward movers stop, which shortens time infectious. Transmission matches its per-step expectation, `1 − e^(−x)` | [04_flux_check](experiments/04_flux_check/) |
| 5 | v0.10 | Does exact per-step sampling (`TBExact`) remove the gap? | Yes: 1.00 / 1.00 / 0.97 ± 0.02 of the ODE at daily / weekly / monthly `dt` (40 reps); even the daily gap was step ordering. Per-step flows match `expm(Q·dt)` | [05_exact_transitions](experiments/05_exact_transitions/) |
| 6 | v0.10 | How much longer does `TBExact` take to run? | No longer: 0.98–1.01× `tbsim.TB`'s run time; both ~10.5 ms per step at 50k agents | [06_profiling](experiments/06_profiling/) |
| 7 | v0.10 | Does dating each new infection within its step (`infection_timing='within'`) remove exp 05's INFECTION artefact? | Yes: monthly INFECTION goes from 0.89 to 0.98 ± 0.02 of the ODE; prevalence, incidence and mortality unchanged. Weekly `'within'` came out 0.96 ± 0.01 vs 1.00 for `'start'`, with no mechanism found | [07_infection_timing](experiments/07_infection_timing/) |
| 8 | v0.10 | With progression by time since infection (`k` = 5.846), does the ODE (Poisson-mixture latent classes) match the ABM? | `TBExact`: yes, 1.00 / 1.01 ± 0.01 at weekly / monthly `dt` for set B; set A 0.99 / 0.94 ± 0.02. `tbsim.TB`: no, 1.18 / 1.83× (set B) and 1.30 / 2.64× (set A), because the start-of-step hazard overstates early progression (cohort +10–19% at monthly `dt`) | [08_progression_tsi](experiments/08_progression_tsi/) |
| 9 | v0.10 | Is set A's near-threshold `TBExact` shortfall (0.94 monthly) a finite-population effect? | No: about 5% at 25k–200k agents. Cause found in exp 10 | [09_population_size](experiments/09_population_size/) |
| 10 | v0.10 | Is the shortfall caused by the per-step infection probability (`1 − e^(−x)`, susceptibles fixed at the start of the step)? | Yes: mass action with `1 − e^(−x)` gives 0.941 ± 0.010 as predicted, and the network the same. A linear probability gives 0.984 ± 0.010 (A) and 1.000 ± 0.006 (B). Still unexplained: `k` = 0 weekly `'within'`, 0.966 ± 0.007 | [10_transmission_step](experiments/10_transmission_step/) |
| 11 | v0.10 | Is `TBExact`'s monthly incidence shortfall a snapshot artefact? | Yes: counted from transitions, incidence per prevalent case is 1.002–1.003 ± 0.001 of the ODE at monthly `dt`. Snapshots miss 3.2–3.3% of episodes at monthly `dt` and 0.75% at weekly | [11_incidence_counting](experiments/11_incidence_counting/) |
| — | v4 | Cliff's v4 redo (2026-10-27): do the key findings of exps 01–11 hold under v4? | The gap is gone without `TBExact`: agents ÷ ODE 0.98–1.04 at weekly / monthly `dt`; one step matches `expm(Q·dt)` to 0.0005; classes and `k_latent` match across solvers. Network: `beta/K` 0.78, matched 1.04. Abbreviated reps, unpinned versions; reproduced in exp 13 | [v4/v4_findings.qmd](v4/v4_findings.qmd) |
| 12 | v4 | Does tbsim v4's ODE reproduce ours, and how should v4 agents be matched to it? | The v4 ODE equals `TB_ODE` and `TBODE` (≤ 3.2e-5) and calibrates to the same `beta`, 9.917338. `RandomNet(n_contacts=10)` has mean degree 10. v4's `beta/K` gives 0.956 ± 0.003 of the ODE's infections per step at monthly `dt`; the matched `beta_edge` 0.997 ± 0.003. `init=` starts agents at the equilibrium. Agents 2.5–5× faster than v0.10 | [12_v4_setup](experiments/12_v4_setup/) |
| 13 | v4 | Do Cliff's findings hold at pinned versions, with more reps and a daily step? | Yes, exactly. Set A 0.99 ± 0.02 / 1.01 ± 0.02 / 1.000 ± 0.009 at daily / weekly / monthly `dt` (Cliff's monthly +4% was chance); set B 0.99. Progression: classes and `k_latent` 0.99–1.01. Latent snapshots +9% (A) / +3–4% (B) at monthly `dt`. Incidence per prevalent case within 0.3%. Network: `beta/K` 0.74, matched 1.03 | [13_v4_pinned_rerun](experiments/13_v4_pinned_rerun/) |
| 14 | v4 | Across tb_LAI_TPT's posterior, is v4's two-state model without remote progression close to r·e^(−kτ), and does its reinfection rule match a clock reset? | Yes, with `k_latent` refitted (≈ 0.95 k): within 0.85% [0.56%, 1.03%] of the exact model (max 2.0%); reinfection within 0.0004. But that version is effectively single-state (clearance raised to about `inf_cle` + 0.95k), so it was set aside. Superseded by exp 16 | [14_posterior_progression](experiments/14_posterior_progression/) |
| 15 | v4 | With remote progression (φ × recent), can the two-state model fit each posterior sample's exact-model targets? | Up to φ ≈ 0.03, as well as the exact model; worse beyond. Superseded by exp 16, because the exact model isn't what tb_LAI_TPT ran | [15_two_state_fit](experiments/15_two_state_fit/) |
| 16 | v4 + v0.10 | Can the two-state model (φ swept) meet tb_LAI_TPT's progression targets within its priors, and how do its dynamics compare with tb_LAI_TPT as run? | Targets met within 1σ for ≥ 92% of posterior weight up to φ = 0.01 with the cascade fixed (74% at 0.03); ≥ 96% up to φ = 0.1 if clearance may rise (tb_LAI_TPT as run: 47%). For φ ≤ 0.03, cumulative progression within 2–6% of tb_LAI_TPT as run. Latent hazard floors at φ·r₀ instead of declining to zero; 27–35% of years 1–10 progression from latent infection, against 12.5% | [16_two_state_vs_lai_tpt](experiments/16_two_state_vs_lai_tpt/) |
| 17 | v4 + v0.10 | With latent reinfection on, do v4's agents match its ODE for the two-state model, and how does the two-state rule change reinfection's effects compared with tb_LAI_TPT as run? | Question 1: yes. Set B, φ = 0.003–0.03: agents ÷ ODE prevalence 0.991 ± 0.006 pooled, as without latent reinfection; no dependence on population size (6,250–200,000 agents; 1.004 ± 0.004); the cleared share has a small finite-population shortfall (≤ 0.6%). Recent latent infection reads ×1.18–1.33 high in monthly end-of-step snapshots. Question 2, cohorts, 12 sets: the two-state rule adds 83–93% as much progression as tb_LAI_TPT's model solved exactly (continuous time), and 75–84% as much as tb_LAI_TPT as run (v0.10, monthly `dt`), about 0.5% of cumulative progression at λ = 0.10. Its progression after each (re)infection comes later (88–91% within a year, against 96%). tb_LAI_TPT's snapshot curve B misses about 11% of entries to asymptomatic TB | [17_latent_reinfection](experiments/17_latent_reinfection/) |
| 18 | v4 (starsim `bc2704b2`, tbsim `efeca56`) | Does Starsim's new ODE of a network make our per-contact `beta` convention unnecessary? | It replaces it: calibrated on the network ODE, the per-contact `beta` equals our matched one exactly, so it's the same model. But both apply the susceptibles' factor to reinfection, which agents transmit at rr·b per contact. Reinfection is 7–8% higher in agents than in the ODE at monthly `dt`, so set B's prevalence is 1.204 ± 0.010 of the ODE's (weekly 1.052 ± 0.009); set A matches. The reinfection factor predicts the gap exactly (1.201) | [18_network_ode](experiments/18_network_ode/) |
| 19 | v4 | Incidence as the first entry to asymptomatic TB per episode: does the ODE match agents, how much does it differ from every entry, and can a standalone copy of the ODE be used for calibration? | Agents match the ODE per prevalent case within 0.4% (0.996–1.002 ± 0.002, sets A and B, monthly and weekly). Re-entries are about half of all entries (46% A, 50% B), so prevalence ÷ incidence is 1.82 years by first entry, 0.9–1.0 by every entry. The standalone ODE (`tbode4.fastode`) reproduces v4's within 1e-14 at equilibrium, 30 times faster. v4's existing `outputs()` misses remote → asymptomatic | [19_v4_incidence](experiments/19_v4_incidence/) |
| 20 | v4 | Calibrate the ODE without treatment (two-state latent model, homogeneous mixing) to tb_LAI_TPT's applicable targets, with IMIS | Mortality ÷ incidence is out of reach without treatment (prior 0.18–0.46 against 0.088–0.144), so it's dropped until treatment is added. The other five targets are all reproduced (prevalence 307 [233–391]/100k against 310; prevalence ÷ incidence 1.60 against 1.57). `beta` 19.4 (10.7–34.1) and π 0.34 (0.21–0.66) are constrained; most other parameters stay broad, and φ isn't identified (re-identification: flat likelihood). Untreated mortality per case is 0.23 (0.15–0.34) | [20_ode_calibration](experiments/20_ode_calibration/) |
| 21 | v4 | Do v4's agents match its ODE across exp 20's calibrated posterior (12 drawn sets)? | Yes. Averaged over the sets: prevalence 0.997 ± 0.007 (monthly) and 0.987 ± 0.007 (weekly); the share asymptomatic, first-entry incidence per prevalent case, prevalence ÷ incidence and mortality ÷ incidence within 0.3%, every set within about 2 SE. Recent latent infection reads ×1.20–1.67 in monthly snapshots, as predicted. Weekly prevalence pooled over every v4 arm is 0.993 ± 0.003 (monthly 1.001 ± 0.003), most likely noise | [21_calibrated_comparison](experiments/21_calibrated_comparison/) |
