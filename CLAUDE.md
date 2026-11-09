# CLAUDE.md

Guidance for Claude Code when working in this repository.

## Purpose

The end goal is an **ODE version of TBsim that is reasonably consistent with the ABM as the user actually runs it**, plus a clear account of which ABM features an ODE can match exactly, which only approximately, and which not at all. The approach is incremental: start from the simplest matched comparison and add one feature at a time. `README.md` holds the per-feature consistency table. Keep it short, and update it as results come in.

**Priorities (from the user):** progression by time since infection (`k_asy`/`k_non`) and latent reinfection matter most; the ODE should reproduce them at least approximately. The household network matters least; it's acceptable to conclude it can't be matched, as long as we describe the gap. Drug resistance (multiple strains) is low priority too (user, 2026-11-04), except that the latent-progression and reinfection structure must carry over to it, since superinfection of a latent infection is latent reinfection.

This is **not** the tbsim repo. Do not edit `../tbsim` unless asked. New ODE code belongs here; it can be proposed upstream later.

**`CONVENTIONS.md`** lists the modelling conventions decided here, written to be carried into a new repo or tb_LAI_TPT. When a decision changes, update it as well as README.

## Upstream files (in `../tbsim`)

- `tbsim/compartmental/lshtm_ode.py`: `default_pars`, `TB_ODE` (odeint; the reference), `TB_SS` (same equations, Euler steps in Starsim). Both share `default_pars`, have identical equations, and were last updated 2026-05-01. They predate the ABM's time-since-infection progression, latent reinfection, and `init_prev_active` (all added July 2026).
- `tbsim/tb.py`: agent-based `TB`. See `define_pars`, `set_prognoses`, `progression_rates`, `step_transitions` and `step_bookkeeping`.
- `tbsim/interventions/treatments.py`: `Tx` (adherence → success; success → possible relapse) and `TxDelivery`.
- `tests/test_compartmental.py`: existing ODE tests and the two-strain ABM↔ODE validation. A good template for matched runs.

**tbsim v4** (`origin/v4`, on starsim v4) is a rewrite. `tbsim.TB` is declared as states and transitions, so one model runs as an ODE (`solver='ode'`) or as agents (`solver='abm'`). The interventions are products and deliveries (`tbsim.care`, `tbsim.regimen`). `TB_ODE` remains as the reference, and `TB_SS`, `Tx` and `TxDelivery` are gone.

Facts from exp 12:

- The v4 ODE equals `TBODE`.
- Agents start at the equilibrium with `init=`.
- On `RandomNet`, use the matched `beta_edge` (`tbode4.beta_edge_from_ode`, the default in `tbode4.make_sim`), not v4's `beta/K`, which is 4% low at monthly `dt`.

**Cliff's v4 redo (2026-10-27).**

- **Where it is:** the scripts at the top of `v4/` (`common.py`, `run_*.py`) and the write-up `v4/v4_findings.qmd`.
- **What it found:** it reran the key findings of exps 01–11. v4 agents match the ODE without `TBExact`, and progression works both with Poisson-mixture classes declared in v4 (`TBClasses`) and with `k_latent`. It didn't cover latent reinfection.
- **Caveats:** abbreviated reps, and unpinned versions with uncommitted changes.
- **Using it:** its code runs in the pinned `v4/.venv`. Build on it rather than redoing it.
- **Exp 13 reproduced it exactly at the pinned versions,** with 40–200 reps and a daily step. `TB_V4_FOLDER` redirects its `results/` to an experiment folder, inherited by worker processes.
- **Seeds:** v4's draws are keyed by seed, so reps with the same seed and unchanged code repeat a run. Extra reps need new seeds.

Where v4's model differs from the v0.10 ABM (front-loading via `k_latent`, latent reinfection, treatment as a rate, network `beta`), see README "Versions and tracks". Decisions so far (user, 2026-11-05):

- **Network `beta`:** the matched per-contact `beta` (proposal 01's preferred option) stays in our own code (`tbode4.beta_edge_from_ode`), not in base Starsim or TBsim, for now. Don't push it upstream.
- **Homogeneous mixing by default (user, 2026-11-08):** v4 agents run without a network unless a feature needs real contacts (household network, contact-targeted interventions). Converted `RandomNet(10, dur=0)` matches it (exp 22) but runs about 2× slower with the same rep-to-rep noise, and needs rr′ recomputed whenever `beta` or `dt` changes.
- **Network reinfection risks (user, 2026-11-08):** on `RandomNet(dur=0)`, also convert each `rr_reinfection_*` (`tbode4.rr_edge_from_ode`), with `trans_asymp` unchanged, until Starsim's ODE of a network handles reinfection (proposal 03, sent to Cliff). Rerun exp 22's `diagnose.py` after any Starsim bump. CONVENTIONS.md 1.
- **Latent progression:** tentatively, v4's two latent states with remote progression = φ × the recent rates (option B in `proposals/02_latent_progression/options.md`), for now.
  - Not the version with no remote progression, which is effectively single-state.
  - φ is to be calibrated eventually, but may only be identifiable with all the targets, burden data included. Until then, sweep a couple of reasonable values (about 0–0.03, up to 0.1).
  - Reopen it if the underlying behaviour lacks face validity against tb_LAI_TPT's.
  - With strains, match v0.10's multistrain model where possible: `rr_reinfection_inf` = `rr_superinfection_latent` = `rr_superinfection_non` = π.
  - Criteria: fit to the data matters much more than consistency with tb_LAI_TPT as run.
  - Whether to propose anything upstream is still open.

## Target ABM configuration (`../tb_LAI_TPT`)

The user's real ABM setup is in `../tb_LAI_TPT/run_one_sim.py`, which calls `utils.make_sim` / `utils.build_sim`. Calibrated parameters are in `results/hm_real_20260728/trajectory_selection/collapsed_samples_eps1.2.csv` (row 0 is the best fit). Main features:

- **Time:** monthly `dt`, 1950–2050.
- **Demography:** `ss.Pregnancy` and `ss.Deaths` driven by India WPP age-specific rates (quarterly `dt`). The initial age distribution comes from WPP, and the population grows.
- **Networks:**
  - `HouseholdNetFixed`: DHS India household sizes 1–13, households split over time via a time-varying `prob_move_out`.
  - `ss.RandomNet(n_contacts=10, dur=0)`: fresh contacts every step, so close to mass action.
  - Each network has its own `beta`: `{'random': beta, 'household': beta * rel_beta_hh}`.
  - Age-mixing pools exist in the code but are commented out.
- **Child transmission:** `NoChildTransmission` sets `rel_trans = 0` for ages under 10.
- **Time-since-infection progression:** `k_asy = k_non ≈ 5.85` per year in the best fit, i.e. strongly front-loaded with an e-folding time of about 2 months. This applies only to exits from `INFECTION` (`inf_non`, `inf_asy`).
- **Reinfection:** `rr_reinfection_rec = rr_reinfection_treat = rr_reinfection_cleared = pi` (calibrated, ≈ 0.24). `rr_reinfection_inf` / `rr_reinfection_non` keep their defaults, so they also equal `pi`.
- **Treatment:**
  - `SymptomaticTxDelivery` treats a time-varying fraction of symptomatic cases, with separate interventions for under-15 and 15+ (`rel_tx_covg_u15`). Coverage is 0 before 1995.
  - `TxVaryEfficacy` has adherence and time-varying efficacy, and a fixed 6-month course.
  - Failure returns the person to their prior state.
  - A success relapses with probability `p_relapse`, after a fixed 1.5 years, moving to `SYMPTOMATIC` (from `CLEARED`, `INFECTION` or `NON_INFECTIOUS`).
- **Mortality:** time-varying `sym_dead` (`TimeVaryingSymDead`).
- **Initial conditions:** `SeedInfections` seeds several states at once; `init_prev = 0.2`.
- **Outputs only:** `over_detection` affects notifications output, not dynamics.

## Prior work (`../tb_LAI_TPT/ode_testing`, May 2026)

Earlier ODE↔ABM work aimed at using an ODE burn-in and handing the result to the ABM. **The folder was deleted** in tb_LAI_TPT commit `bb19486` (2026-10-28); recover files with `git show bb19486^:ode_testing/<file>` there. Its `beta_ode_from_tb` uses the unmatched mapping (K·b, see `CONVENTIONS.md`), so don't reuse it.

- `tb_ode_demog.py`: `TB_ODE_demog`, which adds time-varying crude birth and death rates and lets N vary (TB deaths are not replaced). It is not age-structured, and its docstring says it predates the ABM's treatment/relapse refactor.
- `run_ode_vs_abm.py`: `beta_ode_from_tb()` translates the per-edge monthly betas plus `n_contacts` and average household size into an ODE `beta`. It compares the models before treatment (to 1979).
- `run_ode_handoff_test.py`: seeds the ABM from the ODE state at 1960 or 1970.

Reuse what is useful, but check it against current tbsim first.

## Established facts

- Reinfection of `INFECTION` / `NON_INFECTIOUS` leaves the state unchanged and only resets `ti_infected` (`set_prognoses`). It does inflate `incidence_kpy` / `new_infections` counts.
- `ti_infected` is only read in `progression_rates`, which only applies to `INFECTION` exits. So reinfecting latent people affects dynamics only when `k_asy` or `k_non` > 0. Reinfecting `NON_INFECTIOUS` people never affects state dynamics, because they cannot return to `INFECTION`.

- On `ss.RandomNet(n_contacts=K)` (starsim 3.6.1), K is the **mean degree** (~K/2 edges per agent, transmission both ways). Measured directly; the plugin docs and the class docstring can mislead here. The matched ABM `beta` comes from `tbode.beta_edge_from_ode`.
- `ss.Births` / `ss.Deaths` take rates per 1000 (`rate_units=1e-3`): use `ss.peryear(1000 * mu)`.
- A `tbsim.TB` step runs transmission first, then transitions that apply immediately in the fixed order INFECTION → NON → ASY → SYM. Agents can pass through several states in one step and can leave INFECTION in the step they are infected.

- **Baseline calibration** (user-chosen): single target of 500/100k (range 200–800) active ASY+SYM prevalence, all ages; only `beta` is tuned (`tbode.calibrate_beta`, ODE equilibrium); everything else fixed. With tbsim default rates, `beta` = 9.917, which is near threshold: prevalence is ~10× as sensitive as `beta`, and the approach to equilibrium is slow. Start ABMs at the ODE equilibrium (`init_from_ode`) rather than burning in. `init_from_ode` assigns states only, with `ti_infected = 0`. That's enough while rates have no memory, but it must be extended feature by feature: time since infection for `k` > 0, age, and time in treatment or relapse-pending.
- **Also report:** incidence (first ASY entry per episode, tb_LAI_TPT definition), mortality, % asymptomatic, and progression (ProgressionY1 / ProgressionRatio; targets 0.0523 / 0.663).
- **Both ABM versions in every feature experiment (user decision, 2026-10-09):**
  - Run `tbsim.TB` and `TBExact` side by side, and report `TBExact` ÷ ODE (can the ODE match the feature?) and `tbsim.TB` ÷ `TBExact` (how much does step ordering add?), compared with the baseline gap.
  - Extend `TBExact` alongside `TBODE` for each feature that changes the natural history.
- **Parameters (user decision):** use tbsim's default natural history with the 500/100k `beta` calibration. The goal is to compare ABM and ODE versions with similar dynamics and complexity to tb_LAI_TPT, *not* to reproduce the exact tb_LAI_TPT calibration in the ODE. (For reference: tb_LAI_TPT's calibrated rates only make sense with `k_asy = k_non ≈ 5.85`; without it, 42% of a latent cohort progresses within 10 years.)
- **Step-ordering bias (exp 03–05):**
  - Standard `tbsim.TB` lets agents carry on forward within a step, and stops agents moving backward, cutting time infectious. This is the whole cause of the ABM↔ODE gap: 0.5% / 1.8% / 3.8% of effective `beta` at daily / weekly / monthly `dt`, or 5% / 16% / 33% of prevalence near threshold.
  - `tbode.TBExact` (exact `expm(Q·dt)` sampling; constant hazards only) removes it: prevalence is 1.00 / 1.00 / 0.97 ± 0.02 of the ODE with 40 reps.
  - Transmission (`beta_edge_from_ode`) matches its per-step expectation, `1 − e^(−x)` with x = `beta`·w·`dt`. But that is a fraction of about x/2 below the ODE's flow, because the susceptible pool is fixed at the start of the step (exp 10). Near threshold at monthly `dt`, this costs up to 6% of prevalence (set A). Mass action with a linear per-step probability (`TBExactMassAction(infection_prob='linear')`) matches the ODE to within 1–2%. The random network behaves like mass action.
- **Incidence: use `incidence_event_py` (`TBExact`, counted from transitions), not the snapshot-based `incidence_py`.** The snapshot misses 3.3% of episodes at monthly `dt` and 0.75% at weekly (exp 11). The event count matches the ODE's incidence per prevalent case to within 0.3%. Notifications will need a transition-based count too once treatment is added. **Decision (user, 2026-10-13): the transition-based count is the standard for ODE↔ABM comparisons.** **Decision (user, 2026-11-05): count only the first entry to ASY per episode**, as tb_LAI_TPT does. In v4, use `tbode4.fastode` (ODE) and exp 19's `Incidence` analyzer (agents). `v4/common.py`'s `outputs()` counts every entry, and misses remote → ASY, so it isn't incidence by this definition. Re-entries are about half of all entries (exp 19). Deferred: whether calibration targets should exclude brief subclinical episodes. If so, define the rule by duration in both models, not by step boundaries (see exp 11's SUMMARY).
- **Network `beta` (exp 18):**
  - **Starsim's fix:** since starsim `5ba36ebb`, newer than the v4 pins, Starsim's ODE of a network gives the same model as our matched `beta_edge`.
  - **The gap both share:** neither matches agents for reinfection at a reduced risk. Set B on `RandomNet(10)` runs at 1.20× the ODE's prevalence at monthly `dt`, from 7–8% more reinfection in agents.
  - **Not affected:** homogeneous mixing, as in exps 13–17.
  - **How v4 agents transmit on a network (exp 22):** per edge per step, the transition's rate for its most infectious source saturates (1 − e^(−rr·b·D)). Everything else multiplies that chance linearly: other sources' coefficients such as `trans_asymp`, `rel_trans`, `rel_sus`, and modifiers. So converting b and each `rr_reinfection_*`, with `trans_asymp` unchanged, makes agents on `RandomNet(dur=0)` reproduce the homogeneous ODE (1.000 ± 0.006 monthly over 12 sets).
  - **The exp 18 environment:** tbsim `efeca56` works with starsim up to `bc2704b2`, not later. That environment is `experiments/18_network_ode/.venv`; run with `PYTHONPATH=v4`.
- **tb_LAI_TPT's curve B (`cohort_progression.FirstPassage`) is read from end-of-step snapshots, so it misses agents who reach asymptomatic TB and leave it within a step** (v0.10 can move latent → ASY → NON in one step). At monthly `dt` that's about 11% of entries: ProgressionY1 is 0.060 counted from transitions, against 0.054 from snapshots and 0.052 for tb_LAI_TPT's model solved exactly in continuous time (exp 17, 12 sets). The snapshot loss roughly offsets v0.10's monthly overstatement of progression at 1 year, not at 10. When comparing with tb_LAI_TPT's progression, say which count is used. Exp 17's `run_lai_tpt.py` (`TBConstantFOI`) records every move.
- **Use `TBExact(infection_timing='within')` from exp 08 on (under review).** It dates each new infection uniformly within its step (exp 07), which removes the half-step INFECTION artefact.
  - **Trade-off, found 2026-10-13:** the force of infection is taken from the start-of-step state. So `'within'` dates infections on average half a step after the infectiousness that caused them, making the generation interval half a step too long. `'start'` gets it right on average, but its end-of-step INFECTION snapshots read low.
  - Neither choice affects the equilibrium. Neither changes how many progress: the latent classes shift with the infection date.
  - **Open at the v0.10 freeze:** which timing to use once runs change over time. v4 has its own timing (README, "What v4 changes").

  Still unexplained: `k` = 0 with weekly `'within'` timing runs about 3.5% low (0.966 ± 0.007 vs 0.995 predicted), with or without the network (exps 07, 10). Other near-threshold shortfalls are explained by transmission (exp 10). Compare `TBExact` with a deterministic version of its step (`experiments/10_transmission_step/step_map.py`), using the `'exp'` infection-probability form for network transmission. Runs with the same seed are uncorrelated after 100 years, so pairing by seed doesn't reduce noise.

- **Time-since-infection progression in the ODE: Poisson-mixture latent classes (user decision, 2026-10-12).**
  - With r₀ = `inf_non` + `inf_asy` and `k_asy = k_non = k`, a progression hazard of r₀·e^(−kτ) gives exactly the same progression times as a mixture of classes. A new infection joins class m with probability Poisson(m; r₀/k); class m progresses at constant rate m·k; class 0 never progresses.
  - Clearance and death act identically on every class, so the model is exact, not an approximation. Latent reinfection (clock reset) is also exact: a reinfected agent leaves its class and is reassigned with the same Poisson weights.
  - 4–5 classes are enough at `k` = 5.85 (checked against the exact time-varying cohort, to 1e-6). A chain of time-since-infection stages converges slowly (1-month stages give 4–7% too little progression) and is not used.
  - `TBExact` uses the same classes. Check both independently against `ode_cohort_progression` (exact time-varying hazard) and against `tbsim.TB` at daily `dt`.
  - **Exp 08: `TBExact` matches the ODE** (set B: 1.00 / 1.01 ± 0.01 at weekly / monthly `dt`). **`tbsim.TB` doesn't:** it evaluates the hazard at the start of each step and applies it for the whole step, overstating cohort progression by 10–19% at monthly `dt` (3–5% weekly; exact at daily). Prevalence is then 1.8× (set B) to 2.6× (set A) the ODE at monthly `dt`, and 1.2–1.3× at weekly. tb_LAI_TPT's calibration compensates for this bias. **Decision (user, 2026-10-12): no changes to tbsim.** `TBExact` is the corrected ABM and accumulates the fixes. tb_LAI_TPT will need recalibrating with it to be consistent with the ODE; that's deferred. `k` is not recalibrated: set B, measured exactly, already fits the data-derived progression targets (within 0.25 SD).
- **Parameter sets for the progression experiment (user decision, 2026-10-12):** A = tbsim defaults + `k` = 5.85 (links back to exps 01–06); B = tb_LAI_TPT best-fit natural history (row 0, pre-1995 `sym_dead` = 0.263) + `k` = 5.85. In both, `beta` is recalibrated to 500/100k. B is the reference set for later experiments. B uses tb_LAI_TPT's reinfection risk ratios: `rr_reinfection_cleared = rr_reinfection_rec = rr_reinfection_treat` = π (0.24). Latent reinfection stays off until it is the feature under test (exp 12: `rr_reinfection_inf = rr_reinfection_non` = π). A keeps tbsim's defaults, so near-threshold sensitivity is retained.

## Environment

- **Two tracks, two environments (user decision, 2026-11-03).** See README "Versions and tracks".
  - **v0.10 (frozen, exps 01–11):** the repo root's `.venv` and `tbode/`, with tbsim v0.10.3 @ `b92d7e8` and starsim v3.6.1 @ `3521adb4`.
  - **v4 (active, exps 12 on):** `v4/.venv` and `v4/tbode4/`, with tbsim `origin/v4` @ `dcfb522` and starsim `origin/v4` @ `80c81f53`.

  Run v4 scripts from the repo root with `v4/.venv/bin/python`. If `VIRTUAL_ENV` points elsewhere (e.g. tb_LAI_TPT's), unset it.

- **Both environments pin tbsim and starsim to exact git commits** in `pyproject.toml` and `uv.lock`. Never point them back at the editable `../tbsim` / `../starsim` checkouts: tb_LAI_TPT shares those, and they move.
  - Bump the v4 pins only on purpose, and rerun `v4/tests/test_anchor.py` after each bump.
  - tbsim v4 still reports `__version__ = '0.10.3'`, so only the commit identifies it.
- **Record versions.** Every experiment README has a "Track" line, and every run writes `versions*.json` via `tbode.save_versions` / `tbode4.save_versions`. The two environments can't import each other, so results cross between them as saved files (e.g. exp 12's `reference_v010.json`).
- **v0.10 arms in new experiments (rule revised by the user, 2026-11-05):**
  - **Allowed as the "tb_LAI_TPT as run" reference:** v0.10's `tbsim.TB` at tb_LAI_TPT's settings (monthly `dt`), to see how v4's behaviour changes the original ABM's. Exps 16 and 17 do this, in the pinned v0.10 environment.
  - **Not for v0.10 ABM↔ODE consistency:** that track is frozen. No new `TBExact` or `TBODE` work; compare against saved v0.10 results.
  - tb_LAI_TPT stays on v0.10 for now. New projects will probably use v4, starting from tb_LAI_TPT's parameters.
- **Starsim v4 compiles its kernels once** (`ss.compile()`, ~2 minutes).
- ABM cost (exp 06): about 10.5 ms per step at 50k agents on one core, the same for `tbsim.TB` and `TBExact` and at every `dt`. A 100-year run takes ~13 s monthly, ~1 min weekly, ~6 min daily; parallel runs are somewhat slower each. Runs are on the user's laptop: use at most 4 cores, via `sc.parallelize(..., ncpus=tbode.ncpus)` over seeds.

## Conventions

- **One folder per experiment** under `experiments/`, e.g. `experiments/01_baseline/`:
  - `README.md`, written before running: the question, the design, and the tbsim commit
  - the run script(s)
  - outputs, including the raw results behind every figure (e.g. per-seed time series and ODE runs, saved with `sc.save`), so figures can be redrawn without rerunning
  - `SUMMARY.md`, written after running: the findings
- After each experiment, update `README.md`, keeping its three parts distinct and without repeating results:
  - **the consistency table:** each row's status, with headline numbers only;
  - **the Plan:** what to do next, in order, naming the row each step moves forward, plus open decisions. No "Done" summaries of experiments;
  - **the findings log:** one concise row per experiment.
  - **Decisions** go in `CONVENTIONS.md` (modelling conventions), `CLAUDE.md` (process and facts) or the relevant proposal; the README links to them.
- **Keep each `README.md` and `SUMMARY.md` a coherent, current document, not a log.** When an experiment gains new arms or reps, or a finding changes:
  - rewrite the affected sections so the document reads as one account of the current results;
  - don't append "update" or "superseded" sections;
  - merge new results into the existing tables;
  - say how the runs came about in one short note, and keep superseded numbers out of the main text;
  - make each experiment's `README.md` design describe every arm actually run;
  - update other documents that quote the changed numbers: other experiments' summaries and the findings log.
- **Make each experiment's documents readable on their own.** When a `README.md` or `SUMMARY.md` mentions a method, class or result from another experiment (e.g. `TBExact`, step ordering), explain it in a sentence or two there, as well as linking to it.
- **Make tables self-explanatory:**
  - say what each value is (e.g. "ABM ÷ ODE, 1.00 = match");
  - state which statistic is used (mean ± SE or median) and how many reps;
  - use one statistic consistently within a table.
- Change **one** thing at a time.
- Run several ABM reps (seeds) and show the mean ± SE and a 10–90% interval against the ODE. Use enough reps that the SE is small relative to the gap being tested; near threshold, the rep-to-rep spread is large (about 14% of prevalence at 50k agents).
- Compare proportions, not counts. Use the same outputs on both sides: state proportions, prevalence (`ASYMPTOMATIC + SYMPTOMATIC`), TB mortality, and notifications once treatment is on. Don't compare ABM infection-incidence outputs naively, because they include reinfections.
- Shared code goes in the `tbode/` package, not copied between experiments:
  - `ode.py`: `TBODE`, extended feature by feature, with each feature behind an option defaulting to off
  - `abm.py`: `make_abm` and the `StateCounts` analyzer
  - `compare.py`: outputs, summaries and plots
- With all options off, `TBODE` must reproduce tbsim's `TB_ODE`; keep a regression test for this in `tests/`.
- Experiments are numbered by question, not by table row; one row may need several experiments.
- **Number an experiment or proposal only when its folder is created** (user, 2026-11-06). Refer to planned ones by description (e.g. "a treatment proposal, not yet written"), not by a reserved number, so numbers follow the folders without gaps.

## Style

- Follow the Starsim style guide; use the Starsim-AI skills when writing Python. Max line length 200.
- On a `tbsim.Sim`, use `sim.get_tb()`. Pass TB parameters directly as kwargs.
- **Writing documents** (READMEs, SUMMARYs, proposals, CONVENTIONS.md):
  - **Concise:** say each thing once, and cut what doesn't help the reader.
  - **Plain language, for a scientific audience:** define terms and symbols where they first appear, and don't use the repo's shorthand without explaining it.
  - **Standing on their own:** describe the current design or findings, not how the document got there. Don't refer to earlier versions or drafts of it ("now", "no longer", "revised", "as before"), and give decision dates only where they matter. This extends "a coherent, current document, not a log" under Conventions.
- In Markdown, always leave a blank line before a list, including one that follows a lead-in line such as "**Design:**" or a nested list that follows continuation text. Some renderers otherwise merge the list into the paragraph above it.
