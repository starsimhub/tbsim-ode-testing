# CLAUDE.md

Guidance for Claude Code when working in this repository.

## Purpose

The end goal is an **ODE version of TBsim that is reasonably consistent with the ABM as the user actually runs it**, plus a clear account of which ABM features an ODE can match exactly, which only approximately, and which not at all. The approach is incremental: start from the simplest matched comparison and add one feature at a time. `README.md` holds the per-feature consistency table. Keep it short, and update it as results come in.

**Priorities (from the user):** progression by time since infection (`k_asy`/`k_non`) and latent reinfection matter most; the ODE should reproduce them at least approximately. The household network matters least; it's acceptable to conclude it can't be matched, as long as we describe the gap.

This is **not** the tbsim repo. Do not edit `../tbsim` unless asked. New ODE code belongs here; it can be proposed upstream later.

## Upstream files (in `../tbsim`)

- `tbsim/compartmental/lshtm_ode.py`: `default_pars`, `TB_ODE` (odeint; the reference), `TB_SS` (same equations, Euler steps in Starsim). Both share `default_pars`, have identical equations, and were last updated 2026-05-01. They predate the ABM's time-since-infection progression, latent reinfection, and `init_prev_active` (all added July 2026).
- `tbsim/tb.py`: agent-based `TB`. See `define_pars`, `set_prognoses`, `progression_rates`, `step_transitions` and `step_bookkeeping`.
- `tbsim/interventions/treatments.py`: `Tx` (adherence → success; success → possible relapse) and `TxDelivery`.
- `tests/test_compartmental.py`: existing ODE tests and the two-strain ABM↔ODE validation. A good template for matched runs.

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

Earlier ODE↔ABM work aimed at using an ODE burn-in and handing the result to the ABM:

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
  - Transmission (`beta_edge_from_ode`) is exact.

## Environment

- This repo has its own uv environment (`uv sync`; run with `uv run python ...` or `.venv/bin/python`), with editable `../tbsim` and `../starsim`. Record both commits for each experiment. Starting point: tbsim v0.10.3 @ `b92d7e8`, starsim v3.6.1 @ `3521adb4`.
- ABM cost (exp 06): about 10.5 ms per step at 50k agents on one core, the same for `tbsim.TB` and `TBExact` and at every `dt`. A 100-year run takes ~13 s monthly, ~1 min weekly, ~6 min daily; parallel runs are somewhat slower each. Runs are on the user's laptop: use at most 4 cores, via `sc.parallelize(..., ncpus=tbode.ncpus)` over seeds.

## Conventions

- **One folder per experiment** under `experiments/`, e.g. `experiments/01_baseline/`:
  - `README.md`, written before running: the question, the design, and the tbsim commit
  - the run script(s)
  - outputs, including the raw results behind every figure (e.g. per-seed time series and ODE runs, saved with `sc.save`), so figures can be redrawn without rerunning
  - `SUMMARY.md`, written after running: the findings
- After each experiment, update the consistency table, the Plan section (what has and hasn't been run, open decisions) and the findings log in `README.md`.
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

## Style

- Follow the Starsim style guide; use the Starsim-AI skills when writing Python. Max line length 200.
- On a `tbsim.Sim`, use `sim.get_tb()`. Pass TB parameters directly as kwargs.
