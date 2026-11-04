# tbsim-ode-testing

Building an ODE version of TBsim that is reasonably consistent with the agent-based `tbsim.TB` as it is used in practice (e.g. `tb_LAI_TPT`), and documenting which ABM features can and cannot be matched by an ODE.

## Approach

Start from the simplest matched case: single-strain `tbsim.TB` vs an ODE, with constant demographic rates, mass-action mixing and no treatment. Then add ABM features one at a time, extending the ODE to match each one where possible. For each feature, record whether the ODE can match it, how, and what gap remains.

## Versions and tracks

The work runs on two versions of tbsim, kept apart:

| Track | tbsim | starsim | Environment and code | Experiments | Status |
|---|---|---|---|---|---|
| v0.10 | v0.10.3 @ `b92d7e8` | v3.6.1 @ `3521adb4` | repo root: `.venv`, `tbode/` | 01–11 | Frozen (2026-11-03): rerunnable, not extended |
| v4 | `origin/v4` @ `dcfb522` (unreleased) | `origin/v4` @ `80c81f53` (4.0.0rc1 + 21) | `v4/`: `v4/.venv`, `v4/tbode4/` | 12 on | Active |

How the tracks are kept separate:
- **Pinned commits.** Each environment installs tbsim and starsim from exact git commits, recorded in its `pyproject.toml` and `uv.lock`. It doesn't use the editable `../tbsim` and `../starsim` checkouts, which tb_LAI_TPT shares and which move.
- **Bumping v4.** Bump the v4 pins on purpose, then rerun `v4/tests/test_anchor.py`, which checks the v4 ODE against our v0.10 ODE.
- **Recording versions.** Each experiment's README has a "Track" line, and each run writes the exact versions to `versions*.json` beside its outputs (`tbode.save_versions` / `tbode4.save_versions`).
- **Moving results between tracks.** The two environments can't import each other, so results pass between them as saved files.
- **No new v0.10 arms.** New experiments compare against saved v0.10 results rather than rerunning them (user, 2026-11-03).

**Cliff's v4 redo (2026-10-27).** Cliff reran the key findings of exps 01–11 under v4.
- **Where it is:** the scripts are at the top of `v4/` (`common.py`, `run_*.py`, `analyze.py`), the write-up is [`v4/v4_findings.qmd`](v4/v4_findings.qmd), and the tables are in [`v4/results/summary.md`](v4/results/summary.md).
- **Versions:** it ran on tbsim `v4` @ `7b2e644` and starsim `v4-oct27` @ `fe6c6f9`, both with uncommitted changes. So its exact versions can't be reconstructed.
- **Reps:** cut for time, to 12–40 per arm, with no daily step.
- **Rerunning it:** **exp 13 reran it at our pinned versions and reproduced it exactly.** Every rep with one of Cliff's seeds gives his result. Exp 13 also added reps (40 per arm, 200 where close) and a daily step, so the v4 numbers below are exp 13's.

**What v4 changes.** `tbsim.TB` is declared as states and transitions, so one model runs as an ODE (`solver='ode'`) or as agents (`solver='abm'`). Its ODE is our v0.10 ODE, `tbode.TBODE`, to 1e-5 or better (exp 12). So "v4 agents ÷ v4 ODE" can be set directly beside v0.10's "`tbsim.TB` ÷ `TBODE`" and "`TBExact` ÷ `TBODE`".

The agent solver draws each state's exits jointly from the competing rates, so the order transitions are declared in doesn't matter. Cliff's redo found, and exp 13 confirmed, that this removes v0.10's gaps without anything like `TBExact`:
- **Steady state.** Agents started at the ODE equilibrium stay there. Set A prevalence is 0.99 ± 0.02 / 1.01 ± 0.02 / 1.000 ± 0.009 of the ODE at daily / weekly / monthly `dt`, against 0.95 / 0.84 / 0.67 for v0.10's `tbsim.TB`. Set B is 0.989 ± 0.006 / 0.99 ± 0.02 (weekly / monthly).
- **One step.** One step from each state matches `expm(Q·dt)` to within 0.0005.
- **Incidence.** Incidence counted from transitions is within 0.3% of the ODE per prevalent case.

What remains:
- **Latent infection reads high in end-of-step snapshots** at monthly `dt` (+9% in set A, +3–4% in set B; gone at daily), because infections are dated at the end of their step.
- **Growth from a seed lags the ODE** by 0.1–0.3 years.

**Where v4's model differs from the v0.10 ABM.** We intend to propose changes upstream for each of these:
1. **Front-loaded progression.** v0.10's `k_asy`/`k_non` (progression rates declining as e^(−kτ) with time since infection τ) are replaced by `k_latent`: two latent states, recent and remote, with constant rates.
   - This isn't consistent with the earlier work. We'll propose an alternative, but which one is still open (see Plan).
   - **Both options are exact under both solvers in v4.** Cliff's redo declares the Poisson-mixture classes directly in v4 (`TBClasses` in `v4/common.py`), and agents match the ODE with them. `k_latent` with no progression from the remote state also matches across solvers. With `k_latent` = 5.846 it lands close to the classes on the progression targets: ProgressionY1 0.049 vs 0.051, ProgressionRatio 0.676 vs 0.674, against targets of 0.0523 / 0.663. So the choice is about the model's shape and reinfection, not about matching the ODE.
   - **With `k_latent` refitted to 5.533,** the two-state model is within about 1% of the exact e^(−kτ) model's cohort progression at every time from 3 months to 50 years. Its reinfection rule (remote → recent) adds the same progression as a clock reset, to 0.0001 ([options note](proposals/02_latent_progression/options.md)).
   - **Pitfall:** `inf_non_remote=None` means "the same as `inf_non`". So by default only progression to asymptomatic TB is front-loaded.
2. **Latent reinfection** only moves remote infections back to recent. Reinfecting a recent latent does nothing, where v0.10 reset its clock, and reinfection of NON_INFECTIOUS is gone.
3. **Treatment.**
   - Completion is a rate (`delta`, exponential with a 6-month mean). v0.10's `Tx` and tb_LAI_TPT use a fixed 6-month course.
   - Success comes from competing completion and failure rates, rather than from adherence × efficacy drawn at the start.
   - Failure always goes to SYMPTOMATIC, where v0.10 returned people to their previous state.
   - Relapse is still a fixed 1.5-year delay.
4. **Per-contact `beta` on a network.** v4 uses `beta_ode / K`, which gives 4% fewer infections than the ODE at monthly `dt` (exp 12). Near threshold that is 22% less prevalence (Cliff's redo). Its ODE approximation of a network assumes that mapping. We use the matched `beta_edge`, as in v0.10. The user is raising this with Cliff directly ([proposal](proposals/01_network_beta/proposal.md)).

## Consistency by feature

### Step ordering in v0.10's `tbsim.TB` (affects every row)

`tbsim.TB` v0.10 differs from the ODE in how it advances time, separately from any feature. Within each step, it applies transitions in a fixed order (INFECTION → NON_INFECTIOUS → ASYMPTOMATIC → SYMPTOMATIC). Agents moving forward can carry on through several states in one step, while agents moving backward stop for the step. That cuts time infectious, which acts like a lower `beta`:
- 0.5% / 1.8% / 3.8% less effective `beta` at daily / weekly / monthly `dt`;
- that is, 5% / 16% / 33% less prevalence at the near-threshold baseline calibration.

This is nearly the whole ABM↔ODE gap at baseline. The rest is transmission: each agent's per-step infection probability is `1 − e^(−x)` with a susceptible pool fixed at the start of the step, a fraction of about x/2 below the ODE's flow (exp 10). That matters near threshold at monthly `dt` (up to 6% of prevalence in set A) and is fixed by mass action with a linear probability.

`tbode.TBExact` is a version of `tbsim.TB` that samples each agent's end-of-step state from the exact one-step probabilities `expm(Q·dt)`. It removes the step-ordering gap (exps 03–05) at no extra run time (exp 06). From exp 07 on, `TBExact` also dates each new infection within its step. That corrects end-of-step snapshots but makes the generation interval half a step too long.

v4's agent solver has neither ordering problem: it draws each state's exits jointly (Cliff's redo, exp 13; see above). It keeps the `1 − e^(−x)` transmission step, and it dates infections at the end of the step.

### Features

Each row is an ABM feature of the tb_LAI_TPT setup, added one at a time on top of the rows above it. Status is given for each ABM: v0.10's `tbsim.TB` and our `TBExact` (against `TBODE`), and v4's agents (against v4's ODE). The v4 statuses come from exp 12 and exp 13. Exp 13 reproduces Cliff's redo at pinned versions, with more reps.

Status: ⬜ not started · 🟡 in progress · ✅ matched · ⚠️ approximate only · ❌ not achievable

| Feature (ABM) | Priority | ODE approach | v0.10 status | v4 status | Notes |
|---|---|---|---|---|---|
| Baseline: natural history, single strain, constant birth and death rates, no treatment, `ss.RandomNet(dur=0)` | — | `TBODE` (generalised `TB_ODE`); in v4, `tbsim.TB(solver='ode')`, which equals it. `beta` tuned to 500/100k. The random network, which gives fresh contacts each step, becomes mass action, with the per-contact ABM `beta` set by `beta_edge = -ln(1 - beta·dt/K)/dt`, where K = `n_contacts` is the mean number of contacts | ✅ `TBExact` · ⚠️ `tbsim.TB` | ✅ (exps 12–13) | v0.10, ABM ÷ ODE prevalence at daily / weekly / monthly `dt`: `TBExact` 1.00 / 1.00 / 0.97 ± 0.02 (40 reps); `tbsim.TB` 0.95 / 0.84 / 0.67 (step ordering, above). Per-step transmission is `1 − e^(−x)` in both, about 1% of prevalence at monthly `dt` near threshold (exp 10). v4, agents ÷ ODE prevalence, mean ± SE (exp 13; 40 reps, 200 where marked): set A 0.99 ± 0.02 / 1.01 ± 0.02 / 1.000 ± 0.009 (200) at daily / weekly / monthly `dt`; set B 0.989 ± 0.006 (200) / 0.99 ± 0.02 at weekly / monthly. On the network at monthly `dt`, set A: matched `beta_edge` 1.03 ± 0.02, `beta/K` 0.74 ± 0.02. Exp 12: the ODE equals `TBODE`, and with the matched `beta_edge` infections per step are 0.997 ± 0.003 of the ODE's |
| Age-specific fertility/mortality (WPP) | | Age-structured ODE | ⬜ | ⬜ | v4 has age stratification with aging (`sim.stratify`) and country demography (`ss.people`) for both solvers |
| No transmission from children under 10 | | Age-structured FOI | ⬜ | ⬜ | Requires age structure; in v4, `ss.modify(..., where=age < 10)` |
| Progression by time since infection (`k_asy`, `k_non`) | **High** | v0.10: Poisson-mixture latent classes (exact: class m with probability Poisson(m; r₀/k), progressing at rate m·k; 7–8 classes). v4: to be decided (see Plan) | ✅ `TBExact` · ❌ `tbsim.TB` | ✅ classes · ✅ `k_latent` (exp 13); which structure to use is open | v0.10, `TBExact` ÷ ODE prevalence, set B (tb_LAI_TPT natural history): 1.00 / 1.01 ± 0.01 at weekly / monthly `dt` (40 reps). `tbsim.TB` ÷ ODE: 1.18 / **1.83**, because it applies the start-of-step hazard for the whole step, overstating early progression by 10–19%. Set A (tbsim defaults, near threshold): `TBExact` 0.99 / 0.94 ± 0.02, `tbsim.TB` 1.30 / 2.64 (exp 08). v4, agents ÷ ODE prevalence at weekly / monthly `dt` (exp 13; 40 reps, 200 where marked): set B with the classes 0.993 ± 0.004 (200) / 1.003 ± 0.008; with `k_latent` 1.01 ± 0.01 / 1.01 ± 0.01; set A with the classes 1.00 ± 0.02 (monthly). v4's own `k_latent` isn't the same model as v0.10's `k`; refitted, it is within ~1% of it ([options note](proposals/02_latent_progression/options.md)) |
| Reinfection of latent people (clock reset) | **High** | Reinfected latents are reassigned across the classes with the Poisson weights (exact) | ⬜ | ⬜ | Only matters with time-since-infection progression. v4's reinfection differs (see above). Not in Cliff's redo |
| Treatment (coverage, adherence, efficacy, failure, relapse) | | Treatment and relapse-pending compartments | ⬜ | ⬜ | v0.10 ABM uses fixed delays (6-month treatment, 1.5-year relapse); v4 uses a completion rate (see above) |
| Time-varying `sym_dead`, efficacy, coverage | | Time-varying parameters | ⬜ | ⬜ | v4: `ss.timeseries` with `ss.modify`, for both solvers |
| Household network | Low | Effective `beta` / household-structured ODE | ⬜ | ⬜ | Likely ⚠️ or ❌ |

## Plan

### v4 track (active)

**Done: exp 12 (complete), the v4 ODE and how to set up agents.**
- The v4 ODE equals `TBODE` and tbsim's `TB_ODE`, and calibrates to the same `beta` (9.917338).
- `RandomNet`'s mean degree is still `n_contacts`.
- Agents need the matched `beta_edge`, not v4's `beta/K`.
- `init=` starts agents at the equilibrium exactly.
- Agents run 2.5–5× faster than in v0.10.

**Done: Cliff's v4 redo (2026-10-27).** It reran the key findings of exps 01–11 under v4. The results are summarised under "Versions and tracks" above and in the features table, and in full in [`v4/v4_findings.qmd`](v4/v4_findings.qmd):
- the baseline from the equilibrium (exps 03, 05);
- one step against `expm(Q·dt)` (exps 04–05);
- progression with the classes and with `k_latent` (exp 08);
- incidence from transitions (exp 11);
- growth from a seed (exp 07's timing question);
- the network `beta` (exp 10);
- run time (exp 06).

**Done: exp 13 (complete), Cliff's redo at pinned versions.** We reran Cliff's scripts unchanged in `v4/.venv`, with 40 reps per arm, 200 where the first 40 were about 2 SE out, and a daily step for set A (user decision, 2026-11-03).
- **Reproduced:** every rep with one of Cliff's seeds gives his result.
- **Baseline and progression** (classes and `k_latent`) match the ODE to within 1% at daily, weekly and monthly `dt`.
- **The +4% at monthly `dt` was chance:** set A is 1.000 ± 0.009 with 200 reps.

**Done: [options note on latent progression and reinfection](proposals/02_latent_progression/options.md)** (a comparison, no recommendation). The candidates:
- Poisson-mixture classes (exact);
- v4's two states, as they are and with `k_latent` refitted;
- a chain of recent stages;
- a non-exponential dwell time.

These are compared on shape, reinfection, exactness under both solvers, and size. Main numbers, for set B:
- `k_latent` refitted to 5.533 (from 5.846) puts the two-state model within ~1% of the exact e^(−kτ) model's cohort progression at every time.
- Its reinfection rule adds the same progression as a clock reset, to 0.0001, even at 5× set B's force of infection.

**Next:**
1. **Decision (user): the latent structure and reinfection rule to propose for v4.** The note ends with the questions this leaves open:
   - whether ~1% from e^(−kτ) is consistent enough;
   - whether that holds across tb_LAI_TPT's posterior, not just its best fit;
   - late reactivation;
   - reinfection of recent infections;
   - size once the model is stratified.
2. **Exp 14: latent reinfection, agents vs ODE,** for the structure chosen, with `rr_reinfection_inf` = π. Both candidate rules should be exact under both solvers by construction. This waits for decision 1.
3. **Upstream proposals** (`proposals/`, numbered):
   - **01, per-contact `beta` on networks:** a single [proposal](proposals/01_network_beta/proposal.md), which the user is taking to Cliff directly. It covers Starsim's ODE approximation of a network, and tbsim's docs and `ode_compare.py`. Rerunning tbsim's `ode_compare.py` setting at monthly `dt` showed its network shortfall is the per-contact step effect (0.59 ± 0.05 of the ODE with `beta/K`, 0.99 ± 0.06 matched).
   - **02, latent progression and reinfection:** the options note so far; it becomes a proposal after decision 1.
   - **03, treatment:** a fixed course with outcome probabilities, and failure back to the previous state. To do.
4. **Timing effects to keep in mind** (exp 13).
   - **Latent snapshots:** end-of-step snapshots of latent infection read high at monthly `dt` (+9% in set A, +3–4% in B; about half that at weekly, gone at daily), while the mean over a step is right.
   - **Growth lag:** growth from a seed lags the ODE by 0.1–0.3 years.

   Neither needs action unless a calibration target is a snapshot of latent infection, or the dynamics change fast.

**Parameters.** Sets A and B carry over (exp 08's README):
- **A:** tbsim defaults + `k` = 5.846;
- **B:** tb_LAI_TPT's best-fit natural history + `k` = 5.846.

tb_LAI_TPT itself stays on v0.10 for now (user, 2026-11-03), but new projects will probably use v4, with tb_LAI_TPT's parameters as a starting point. Its calibration partly compensates for v0.10's monthly-step biases, so `beta` is recalibrated to 500/100k as before. The progression targets need rechecking under whichever latent structure is chosen.

**Later rows.** Starting agents at an equilibrium is simple while rates have no memory. Age structure needs an initial age distribution, and fixed delays (treatment, relapse) need initial time in state.

### v0.10 track (frozen 2026-11-03)

Exps 01–11 are complete; results are in the findings log, and the decisions they led to are in `CLAUDE.md`. Their main conclusions:
- `TBExact` matches the ODE at baseline and with time-since-infection progression.
- v0.10's `tbsim.TB` doesn't. Step ordering costs up to a third of prevalence, and the start-of-step hazard gives up to 2.6× prevalence at monthly `dt`.

The decision not to change tbsim (2026-10-12) stands for v0.10.

Items open at the freeze, and where they go:
- **Infection timing in `TBExact`** (`'start'` vs `'within'`). v4 has its own timing (end of step); Cliff's redo measured its effects (Plan item 4).
- **Linear mass-action transmission as `TBExact`'s default.** Moot in v4: its agents match the ODE without it (Cliff's redo).
- **Checking that latent classes don't change other infection dynamics** (time in INFECTION, time to progression, by outcome). This applies to whichever latent structure is chosen.
- **Recalibrating tb_LAI_TPT with `TBExact`.** Deferred, and probably superseded by v4 for new projects.
- **Re-profiling `TBExact` with latent classes.** Dropped.

## Open issues (v0.10 track)

Differences between `TBExact` and the ODE that are not fully explained. They are listed in [exp 10's SUMMARY](experiments/10_transmission_step/SUMMARY.md), with the history in the experiments cited. Not pursued further unless they reappear in v4.

| Issue | Size | Status |
|---|---|---|
| `k` = 0, weekly `dt`, `'within'` timing: prevalence low (exps 07, 10) | 0.966 ± 0.007 of the ODE (320 reps pooled), against 0.995 predicted by a deterministic version of the step (z ≈ −4) | **Unexplained.** Ruled out: the network, per-step sampling, the random-number source, the per-step infection probability, and `'within'`'s generation-interval lag (the deterministic step includes it). Population size was tested only in a different arm (exp 09). The same setting with `'start'` timing matches, as do monthly `'within'` and set A weekly. Doesn't affect set B or monthly runs |
| CLEARED slightly low (exp 08) | 0.998–0.999 of the ODE, z = −2 to −5, even where prevalence matches | Not explained by the step structure: the deterministic step (exp 10's `step_map.py`, `'within'`, `1 − e^(−x)`) predicts 0.999–1.001 (checked 2026-10-13). Otherwise not examined; small |

## Repo layout

```
tbode/          v0.10 track: the ODE (ode.py), matched-ABM builders and analyzers (abm.py), tbsim.TB with exact
                per-step transitions (exact.py), comparison helpers (compare.py), per-step flux checks (flux.py), versions.py
tests/          v0.10 regression checks (e.g. TBODE with all features off reproduces tbsim's TB_ODE)
v4/             v4 track: its own pyproject.toml and .venv. At the top: Cliff's v4 redo (common.py, run_*.py, analyze.py,
                v4_findings.qmd, results/summary.md)
  tbode4/       matched v4 sims as ODE or agents (sims.py), versions.py
  tests/        test_anchor.py: the v4 ODE reproduces TB_ODE and TBODE (rerun after bumping the v4 pins)
proposals/      Drafts of changes to propose upstream (tbsim, starsim), numbered like experiments (NN_name/): one draft
                issue per repo, with the scripts and outputs behind them
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
| 2 | v0.10 | Does a smaller ABM time step close the gap? | Largely, well above threshold: within ~3% at weekly and daily `dt` at `beta` = 12 and 20. But not entirely: near threshold even daily `dt` is 5% low (exp 03). Exp 05 later showed the whole gap is a time-step effect, namely step ordering | [02_timestep](experiments/02_timestep/) |
| 3 | v0.10 | With `beta` calibrated to 500/100k and the ABM started at the ODE equilibrium, does the ABM stay there, at daily, weekly and monthly `dt`? | No: it drifts down at every `dt`, to 0.95 ± 0.02 / 0.84 / 0.67 of the ODE prevalence at daily / weekly / monthly `dt`. Even daily is significantly low (z = −2.9). The calibration sits near threshold (`beta` = 9.92 vs threshold 8.59), so prevalence is ~9× as sensitive as `beta`; the gaps are 0.5% / 1.8% / 3.8% in effective `beta` | [03_baseline_calibrated](experiments/03_baseline_calibrated/) |
| 4 | v0.10 | Which per-step flows cause the gap: transmission or natural-history transitions? | Natural history. Transmission matches its per-step expectation, `1 − e^(−x)`; exp 10 later showed that expectation is slightly below the ODE's flow. With the step ordering, forward movers carry on within a step and backward movers stop, which shortens time infectious | [04_flux_check](experiments/04_flux_check/) |
| 5 | v0.10 | Does exact per-step sampling (`TBExact`) remove the gap? | Yes: prevalence 1.00 / 1.00 / 0.97 ± 0.02 of the ODE at daily / weekly / monthly `dt` (40 reps). Daily `TBExact` (1.00) vs daily `tbsim.TB` (0.95) shows even the daily gap is step ordering. Per-step flows match `expm(Q·dt)` (sum of z² 18–19 over 14 flows, vs 691–2,688 for `tbsim.TB`) | [05_exact_transitions](experiments/05_exact_transitions/) |
| 6 | v0.10 | How much longer does `TBExact` take to run than `tbsim.TB`? | No longer: 0.98–1.01× the run time at every `dt`, within rep-to-rep noise. Both take ~10.5 ms per step at 50k agents, so run time scales with the number of steps (daily ≈ 30× monthly). About 43% of each step is the TB module | [06_profiling](experiments/06_profiling/) |
| 7 | v0.10 | Does dating each new infection within its step (`TBExact`, `infection_timing='within'`) remove exp 05's INFECTION artefact? | Yes: monthly INFECTION goes from 0.89 to 0.98 ± 0.02 of the ODE, and relative to prevalence from 0.91 to 1.00 (weekly: 0.98 to 1.00). Monthly prevalence, incidence and mortality are unchanged. At weekly `dt`, `'within'` came out 0.96 ± 0.01 vs 1.00 for `'start'` (160 reps each). Diagnostics found no mechanism or bug: per-step flows are exact; a deterministic version of the step predicts 0.999; two diagnostic variants agree within noise. The gap is probably mostly chance, but up to ~3% of prevalence near threshold can't be ruled out | [07_infection_timing](experiments/07_infection_timing/) |
| 8 | v0.10 | With progression by time since infection (`k` = 5.846), does the ODE (Poisson-mixture latent classes) match the ABM? | `TBExact`: yes. Prevalence 1.00 / 1.01 ± 0.01 of the ODE at weekly / monthly `dt` for set B (tb_LAI_TPT natural history); set A (near threshold) 0.99 / 0.94 ± 0.02. Daily `tbsim.TB` cohorts match the exact progression curve, which confirms the latent classes independently. `tbsim.TB`: no. 1.18 / 1.83× the ODE (set B) and 1.30 / 2.64× (set A), because the start-of-step hazard overstates early progression (cohort +10–19% at monthly `dt`). The step-ordering gap reverses and grows | [08_progression_tsi](experiments/08_progression_tsi/) |
| 9 | v0.10 | Is the near-threshold `TBExact` shortfall (set A, monthly: 0.94 of the ODE vs 0.995 predicted by a deterministic version of the step) a finite-population effect? | No. It stays about 5% at 25k / 50k / 100k / 200k agents (0.947 / 0.924 / 0.964 / 0.941 ± ~0.01). A constant fits (χ² = 7 on 3 df); 1/N doesn't (χ² = 60). Rep-to-rep SD shrinks with N, but the mean shortfall doesn't. Cause found in exp 10 | [09_population_size](experiments/09_population_size/) |
| 10 | v0.10 | Is the near-threshold shortfall caused by the ABM's per-step infection probability (`1 − e^(−x)`, with a susceptible pool fixed at the start of the step) rather than the ODE's flow (x)? | Yes. Set A monthly: mass action with `1 − e^(−x)` gives 0.941 ± 0.010, exactly as predicted, and the network gives the same (0.947), so the network adds nothing. With a linear probability: 0.984 ± 0.010 (A) and 1.000 ± 0.006 (B). Exp 05's monthly 0.97 is explained (rerun: 0.980 ± 0.010 vs 0.997 predicted). Still unexplained: `k` = 0 weekly `'within'`, 0.966 ± 0.007 pooled vs 0.995 | [10_transmission_step](experiments/10_transmission_step/) |
| 11 | v0.10 | Is `TBExact`'s monthly incidence shortfall (about 3% below prevalence) a snapshot artefact, and does counting from transitions fix it? | Yes. Counting each agent's within-step visits to ASYMPTOMATIC/SYMPTOMATIC (augmented chain, own random stream, dynamics bit-identical) gives incidence equal to prevalence ÷ ODE in every arm. Incidence per prevalent case is 1.002–1.003 ± 0.001 of the ODE at monthly `dt` and 1.000 at weekly. The snapshot misses 3.2–3.3% of episodes at monthly `dt` and 0.75% at weekly, whatever the parameters | [11_incidence_counting](experiments/11_incidence_counting/) |
| — | v4 | Cliff's v4 redo (2026-10-27): do the key findings of exps 01–11 still hold under v4? | The gap is gone without `TBExact`. Agents ÷ ODE prevalence at weekly / monthly `dt`: set A 0.98 / 1.04 ± 0.02, set B 0.98 / 1.00 (± 0.01–0.02). One step matches `expm(Q·dt)` to within 0.0005. Progression with Poisson-mixture classes (`TBClasses`) or `k_latent` matches across solvers (set B 0.98–1.02). Incidence per prevalent case is 0.996–1.006. What remains: latent snapshots read high at monthly `dt` (+3 to +13%), and growth lags by 0.1–0.3 years. On the network, `beta/K` gives 0.78 and the matched rate 1.04 ± 0.04. Abbreviated runs (12–40 reps, no daily step), at unpinned versions with uncommitted changes. Reproduced exactly at pinned versions, with more reps, in exp 13 | [v4/v4_findings.qmd](v4/v4_findings.qmd) |
| 12 | v4 | Does tbsim v4's ODE reproduce ours, and how should v4 agents be matched to it? | The v4 ODE equals `TB_ODE` (2.8e-5) and `TBODE` (≤ 3.2e-5 relative after the first year), and calibrates to the same `beta`, 9.917338. `RandomNet(n_contacts=10)` still has mean degree 10 with fresh edges each step. v4's per-contact `beta_ode/K` gives 0.956 ± 0.003 of the ODE's infections per step at monthly `dt` (per-edge saturation, predicted 0.960); the matched `beta_edge` gives 0.997 ± 0.003. `init=` starts agents at the equilibrium exactly. Agents run at 1.6–4 ms per step at 50k, 2.5–5× faster than v0.10 | [12_v4_setup](experiments/12_v4_setup/) |
| 13 | v4 | Do Cliff's v4 findings hold at pinned versions, with 40 reps per arm (200 where close) and a daily step? | Yes, and they reproduce exactly: every rep with one of his seeds gives his result. Agents ÷ ODE prevalence: set A 0.99 ± 0.02 / 1.01 ± 0.02 / 1.000 ± 0.009 at daily / weekly / monthly `dt`. Cliff's monthly +4% was chance (his 40 seeds 1.043; seeds 40–199 0.990). Set B 0.989 ± 0.006 / 0.99 ± 0.02 (weekly / monthly). Progression: classes 0.993 ± 0.004 / 1.003 ± 0.008, `k_latent` 1.01 / 1.01 (set B). Latent snapshots +9% (A) and +3–4% (B) at monthly `dt`, gone at daily. Incidence per prevalent case within 0.3%. Network: `beta/K` 0.74, matched 1.03. 1.3–2.4 ms per step at 50k agents | [13_v4_pinned_rerun](experiments/13_v4_pinned_rerun/) |
