# 30: Do v4's agents match its ODE with TB, treatment and the corrected demography?

**Status:** design, not run.
**Track:** v4 (`v4/.venv`; tbsim `origin/v4` @ `619c296`, starsim `origin/v4` @ `d83f4f87`).
**Run on:** zebra.

## Why this experiment

Exp 28 compared v4's agents with v4's ODE with TB, treatment and India's fertility and mortality from the UN World Population Prospects (WPP) 2024. Ratios per prevalent case matched within 1%, but population size and the age pattern didn't:

- agents' population was 6–7% below the ODE's from 1990 to 2010;
- prevalence by age group differed by up to 10% with WPP's 5-year age groups in the ODE, and up to 5% with single years;
- all-age prevalence was 1–3% low in agents before 2025.

Exp 29 explained the population gap fully, without TB, and found settings that remove it ([CONVENTIONS.md 10](../../CONVENTIONS.md#10-demography-births-ages-and-the-starting-population)):

- **Agents:** `ss.Pregnancy` with fertility shifted younger by its mean gestation (0.769 years), so births come at WPP's ages; the steady-state share of women pregnant at the start; and births as female as the ODE's (48.63%).
- **The ODE:** age groups of a quarter or an eighth of a year, since aging by constant-rate exits between groups spreads a cohort's ages (single years overstate population by 2.2% by 2000); and WPP's split by sex within each age group at the start, instead of Starsim's default of half of each age female.

With these, agents' population was within 0.1% of the ODE's (extrapolated to exact aging) through 2050, and within about 1% in each broad age group's share.

Exp 28 also started both solvers from an equilibrium without ages, the same TB state shares at every age, which took the ODE 33–77 years to settle and overlapped its first comparison window. This experiment starts from a state with ages instead.

## Question

With the corrected demography and a starting state with ages, do agents match the ODE with exp 28's TB model and treatment: in all-age prevalence, ratios per prevalent case, and prevalence by age group?

## Design

Exp 28's stage 2, with two changes: the corrected demography, and the starting state.

- **The model:** exp 27's, `tbode4.TBTreat`:
  - v4's `tbsim.TB` with the two-state latent model, every reinfection risk at π, and homogeneous mixing;
  - a fixed 6-month treatment course from the end of 1994, with coverage, efficacy and TB mortality (`sym_dead`) over time as calibrated in exp 24;
  - the same rates at every age (no child-specific rates or mixing).
- **Parameter sets:** exp 28's two, sets 0 and 6 from exp 24's posterior (one from each stretch of its ridge; set 0 front-loads progression strongly, `k_latent` = 11.1 per year, and set 6 weakly, 2.4). Not recalibrated with demography.
- **Demography:** WPP 2024 as bundled with Starsim, interpolated linearly by year; deaths by age and sex (`ss.death`) in both solvers; `mu` = 0, and TB deaths aren't replaced.
  - **Agents:** exp 29's P2b: the shifted fertility table, `PregnancyStart`, and `p_female` = 48.63%.
  - **The ODE:** `ss.birth` by age and sex, WPP's split by sex within each age group (`joint_sex()`), and the age-group width chosen below. The rates stay on WPP's groups.
- **The starting state, with ages:** the ODE's untreated steady state at 1950's demography, by age and sex.
  - **How:** run the ODE with eighth-year ages, 1950's fertility and mortality frozen, no treatment and `sym_dead` at its pre-1995 value, from exp 28's starting state, for 200 years (exp 28's B2 settled within 77; checked by the largest change in any group's share over the last 20 years). Average each age and sex group's TB state shares to single years of age, weighted by the frozen run's people, and start both solvers from that table, on WPP's 1950 population (`init=` with `ss.by(age, sex)`, every state but susceptible, which takes the rest).
  - **Why single years:** the ODE takes a table as coarse as its groups or coarser, so one table serves every width. Agents need age groups no coarser than the table's bins, so agents are stratified by single years of age; their groups only report results and read the table. Two Starsim behaviours found in the quick test, both avoided this way: agents given a table finer than their age groups start with wrong shares when several states are given (most adults susceptible, with WPP's groups and an eighth-year table), and they refuse a share by group for every state, susceptible included ("the initial values sum to 0").
  - **What it isn't:** with frozen rates the population's own age distribution moves towards the stable one, and the force of infection with it. So the shares are those of the frozen model's age distribution, not exactly an equilibrium at WPP's 1950 ages, and the runs settle a little after 1950. Both solvers start from the same state, so agents ÷ ODE is still a fair comparison; the ODE alone shows how far it moves before 1995.
  - **For incidence:** the share of non-infectious people already asymptomatic this episode, from the frozen run's NON0 (non-infectious, not yet asymptomatic this episode; tracked as in exp 28), for the ODE's first-entry incidence (NON0 at the start) and agents' `Incidence` flag (drawn at random among non-infectious people, the same share at every age).
- **The quick test** (laptop, before running):
  - **The ODE** starts exactly at the table at every width (largest difference in any group and state about 5e-16).
  - **Agents** start at it within sampling noise in every state and age band up to 85 (1.5 million agents); at 85+, a few cells of 1–5 expected people are 0.
  - **Agents' population with TB** is 0.9997 ± 0.0018 of the eighth-year ODE's in 1990 (8 reps of 100,000; births, all deaths and TB deaths each within SE). Stratifying agents by single years leaves their demography unchanged, and so does the TB module without infection.
  - **`tbode4`'s demography** (moved from exp 29) reproduces exp 29's fertility table, ages and female share exactly.
- **The ODE's age groups:** single, quarter- and eighth-year groups, all affordable with TB (in the quick test, 1950–1998 at weekly steps took 11, 15 and 22 s). Quarter and eighth give the ODE extrapolated to exact aging (2 × eighth − quarter), as in exp 29; single years link back to exp 28's arm S. The rates stay on WPP's groups.
- **Agents:** monthly and weekly steps, 20 reps each per set, new seeds; 100,000 agents at the start. Incidence and the pitfall checks are counted with `tbode4.Incidence` and `tbode4.TxChecks`. Results by age group are counted in our own analyzer (`AgeCounts`) with 1e-3 years added to ages (agents' 32-bit ages otherwise put ±1.7% bands in 5-year groups; exp 29).
- **Years:** 1950–2041, as in exp 28.
- **Size:** 2 sets × 2 steps × 20 reps = 80 agent runs, about 52 s each monthly and 168 s weekly in exp 28 (2.4 core-hours), plus 12 ODE runs of up to about a minute each and two frozen runs. On zebra, in parallel.

## Measured

As exp 28's stage 2, with the ODE at each width:

- **Agents ÷ ODE:**
  - active prevalence (asymptomatic + symptomatic) per person, all ages and in five age groups (0–14, 15–29, 30–49, 50–64, 65+), and the share asymptomatic;
  - per prevalent case: incidence (first entry to asymptomatic TB per episode, plus relapses; CONVENTIONS.md 5), notifications, the number on treatment, and TB mortality (from symptomatic TB and on failing courses);
  - recent and remote latent infection, and the cleared share;
  - the population's size and age shares.

  Each as a yearly mean, then over three windows: 1950–1994 (before treatment), 1995–2024 (the scale-up) and 2025–2040. Mean ± SE across the 20 reps, for each set and step.
- **The starting state:** the 1950 shares by age and sex in each solver against the intended ones; and, in the ODE, how far prevalence by age group moves before 1995 (against exp 28's 33–77 years from an age-independent start).
- **Against exp 28:** the same ratios side by side, to show what the corrected demography and the start change.
- **Run time** of the ODE at each width, with TB.
- **Expected effects, reported but not counted as mismatches:**
  - recent latent infection reads high in monthly end-of-step snapshots (×1.6 in exp 28; exp 21);
  - at monthly steps, treatment starts about 9% low in 1995 in agents (time-varying rates applied half a step earlier; exp 27), and agents drift by about 1% from the ODE's starting state.

**Figures** (raw outputs saved with `sc.save`), for each set, agents' reps and mean against the ODE, monthly and weekly:

- burden over time (prevalence, incidence, TB mortality and notifications per 100,000);
- prevalence by age group over time;
- the population's size and age shares;
- agents ÷ ODE by output and window, with SE, next to exp 28's.

## What would answer the question

- **TB matches with demography** if agents ÷ ODE is within about 2 SE of 1, and within about 1%, for all-age prevalence, each ratio per prevalent case and population size, in the two later windows at weekly steps (as in exps 27 and 28); and prevalence in each broad age group is within about 2%.
- **If prevalence by age still differs** while population and age shares match, the cause is in how TB states move with age, not in demography: compare the ODE's widths to see whether its aging is the cause.
- **If the before-1995 window still falls short,** compare with the ODE's own movement from the starting state, to see whether it's settling.
- **Then the demography row** of the consistency table can move from ⚠️, and recalibration with demography is next.

## Not covered

- **Recalibration** with demography (next): the sets' fits to exp 24's targets move with it.
- **Child transmission, child-specific rates and age mixing,** and lower treatment coverage for children.
- **WPP's sex ratio at birth:** the ODE's births take its starting female share (CONVENTIONS.md 10).
- **Migration.**

## Files

- `run.py`: `check` (the quick test of the starting state), `timing`, `start` (the starting states → `start.obj`), `tb` (the ODE at each width → `ode.obj`, then the agents → `runs.obj`); summaries and figures to add once the runs are back. Every run's raw outputs are saved with `sc.save`, with `versions.json`.
- `quick/`: the quick test of every step (`EXP30_QUICK=1`).
- `SUMMARY.md`: findings, written after running.
- Uses exp 27's sets and treatment schedules, exp 28's demography and outputs, `tbode4.demog` (exp 29's shifted fertility table, `PregnancyStart` and `joint_sex()`, moved there), and `tbode4`'s `TBTreat`, `Incidence`, `TxChecks`, `ode_incidence` and `fastode`.
