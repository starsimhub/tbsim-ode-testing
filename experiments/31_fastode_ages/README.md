# 31: A fast standalone ODE with ages and demography, for calibration

**Status:** run (findings in [SUMMARY.md](SUMMARY.md)).
**Track:** v4 (`v4/.venv`; tbsim `origin/v4` @ `619c296`, starsim `origin/v4` @ `d83f4f87`).
**Run on:** the laptop (4 cores); every run here is cheap.

## Why this experiment

Exp 24 calibrated the ODE with treatment using `tbode4.fastode`, a standalone copy of v4's TB ODE with no ages. In it, births replace every death at μ = 1/70. One draw, from its untreated equilibrium in 1994 to 2040, took about 90 ms, and each IMIS run took 310,000 draws. India's demography raises the ODE's burden 1.2–1.6× at exp 24's sets (exp 28), so the calibration needs redoing with it.

The demography is as exps 29–30 matched it ([CONVENTIONS.md 10](../../CONVENTIONS.md#10-demography-births-ages-and-the-starting-population)): WPP 2024 fertility and mortality by age and sex from 1950, and a starting state with ages. Through Starsim's ODE, with the quarter-year age groups that convention uses, one draw would take about 15 s:

- **the start:** a 200-year run at 1950's demography, frozen, for each draw's own starting shares by age (up to 11 s);
- **1950–2041 at monthly steps:** 7 s (exp 30, on zebra).

That's about 11 days per IMIS run on 4 cores. A standalone copy with ages, written like `fastode`, should be much faster. This is possible because TB's rates are the same at every age and mixing is homogeneous: every TB flow is the same operation on each age group, and ages are linked only through the force of infection (one number per step) and births.

**Monthly steps are accurate enough for calibration in Starsim's quarter-year ODE.** In exp 30's saved runs, its monthly and weekly runs differ by at most 0.2% in prevalence, incidence, TB mortality and notifications, 1960–2040, at both sets. The targets' σ are about 10%.

## Question

Can a standalone ODE with ages and demography reproduce v4's ODE (exp 30), and its exact-aging limit, quickly enough for IMIS? With exact aging, which time step is accurate enough?

## Design

### The model

`tbode4.fastode_age`, a new module; `fastode` itself is unchanged. The TB model is exp 30's:

- **Natural history:** v4's `tbsim.TB` with the two-state latent model (CONVENTIONS.md 2), every reinfection risk at π (convention 3), and homogeneous mixing. The rates are the same at every age.
- **Treatment:** CONVENTIONS.md 9's design. Coverage, efficacy and `sym_dead` vary over time as in exps 24 and 30, and there is over-detection of notifications.
- **Counting:** non-infectious TB is split into NON0 and NON1, as in `fastode`, for first-entry incidence (convention 5).

Its state is an array of TB state × sex × age group. Demography follows exp 30's ODE:

- **Deaths:** WPP's death rates by age and sex, interpolated linearly between WPP's decades (`tbode4.WPP`), at every TB state alike. `mu` = 0, and TB deaths aren't replaced.
- **Births:** WPP's fertility applied to women by their WPP age group. Births go into susceptible in the first age group.
- **Female share at birth:** the starting population's female share (48.63%). Starsim's ODE uses this share with `joint_sex()` (CONVENTIONS.md 10), so the copy keeps it; WPP's sex ratio at birth is for later.
- **Ages:** groups up to 100, with an open 100+ group.

It has two ways of aging, as two functions.

**Option 1, `run_rate()`: a copy of v4's ODE.**

- Quarter-year age groups, with aging out of each group at rate 1 per group width, as Starsim's ODE does. The rates stay on WPP's groups.
- It's a system of ODEs solved with an adaptive solver (`solve_ivp`), as `fastode.run` is.
- Time-varying inputs (WPP rates, coverage, efficacy, `sym_dead`) are functions of continuous time.
- The fixed 6-month course is 12 Erlang stages, as in `fastode` and exp 24 (within 0.1% of the exact fixed course's prevalence; exp 23).
- **Purpose:** to check the copy's TB and demography equations against Starsim's ODE, apart from how it ages people.

**Option 2, `run_cohort()`: exact aging, for calibration.**

- **Groups one step wide:** age groups are one time step h wide. At the end of each step, every group moves up by one group, so a cohort keeps its exact age. This removes the spread of ages that constant-rate aging gives (exp 29).
- **TB within a step:** the generator Q is the matrix of TB transition rates, about 35×35 with treatment. The step applies its matrix exponential, expm(Q·h), to every age and sex group at once. Q is the same at every age, so this is one matrix exponential per step.
- **Force of infection:** taken at mid-step by a predictor-corrector. First predict the end of the step with the force of infection at its start. Then redo the step with the mean of the start's and the predicted end's force of infection. This makes the error proportional to h² instead of h. The start-of-step version is also run, to show what the corrector buys.
- **Inputs over time:** time-varying inputs are taken at mid-step, as v4's ODE applies them.
- **Deaths:** each group survives the step with probability exp(−m·h), where m is its death rate at mid-step.
- **Rates by age over a step:** during a step, a group's people age from [g·h, (g+1)·h) to [(g+1)·h, (g+2)·h), spending exactly half their time in each. So each group takes the mean of its own death and fertility rates and the next group's. With WPP's age edges on the grid, that's exact. The rate at the group's lower edge alone would apply each change with age half a step late: population 1.2% low by 2020 at monthly steps.
- **Births:** births during the step enter the first group at its end. They are counted from the mean of the women at the step's start and end. They survive half a step at the first group's death rate, and get half a step of the step's TB transitions from susceptible (expm(Q·h/2)), since they were born on average half a step before its end.
- **The predicted end of the step** counts these births among the living, the force of infection's denominator. Without them the mid-step force of infection is too high by about half a step's births (0.2% at monthly steps), a first-order error.
- **The fixed course:** a conveyor of whole steps, as in Starsim's ODE. With K the course in whole steps, there are K + 1 slots. People starting treatment in a step enter slot 0, every slot moves along one at the step's end, and the last slot is released over the next step: half at its start, half at its end, standing in for Starsim's even release over the step. A course then lasts K steps on average. Failing courses face `sym_dead` in every slot. K is 2, 6, 12 and 26 at quarter-year, monthly, half-month and weekly steps.
- **The steps:** quarter-year, monthly, half-month (1/24 year) and weekly. The weekly step is 1/52 year, not Starsim's 7 days, so that WPP's age edges fall on the grid (which the rates by age above need) and the course is exactly 26 steps.
- **Incidence:** first entries to asymptomatic TB, counted from each step's transitions, plus relapses. NON0 is tracked by age, so this needs no approximation. Exp 30's `ode_incidence` decays NON0 at its mean death rate over ages.

### The start

Every run starts in 1950 from a table: each TB state's share in each single year of age and sex, put on WPP's 1950 population. Exp 30 made this table with a run of Starsim's ODE:

- **The frozen run:** start from exp 27's untreated equilibrium without ages, and run 200 years with 1950's demography frozen, no treatment, and `sym_dead` at its 1995 value, until the shares by age settle.
- **The table:** average the frozen run's final shares to single years of age, weighted by each group's people.

The table is the TB epidemic each parameter set would have reached by 1950, so it differs between parameter sets. A higher `beta`, for example, means more people latently infected at every age. Exp 30 made tables only for sets 0 and 6. In the calibration, every IMIS draw has different parameters, so every draw needs its own frozen run to make its own table. Running each through Starsim's ODE would take up to 11 s, so the calibration will do the frozen runs in `fastode_age`. This experiment checks those frozen runs separately from the 1950–2041 dynamics:

- **Dynamics (arms C, X, X0):** start from exp 30's saved tables (`start.obj`) for sets 0 and 6. Differences from exp 30's runs then come from the dynamics alone.
- **Building the table (arm S):** make the table with `fastode_age`'s frozen run and compare it, and the 1950–2041 run started from it, with exp 30's.
- **Burn-in length (arm S):** 200, 100 and 50 years. Each is judged by the largest change in any group's share over its last 20 years, and by its effect on the targets. A shorter burn-in makes every draw cheaper. Exp 28 found the ODE took 33–77 years to settle from a start without ages.
- **No endemic state:** a draw whose frozen run dies out, with prevalence below 1e-9, has no endemic state. The calibration gives it zero likelihood.

### Arms

Exp 24's sets 0 and 6 (exp 27's sets file), with treatment, 1950–2041:

| Arm | Purpose | Model | Start | Compared with |
|---|---|---|---|---|
| C ("copy") | Are the TB and demography equations right? Option 1 ages people as Starsim's ODE does, so any difference is an error in the equations | Option 1, quarter-year groups | Exp 30's table | Exp 30's quarter-year Starsim ODE, weekly and monthly steps (`ode.obj`) |
| X ("exact aging") | Does option 2 reach exact aging, and how small can its step be before it drifts? This picks the step for the calibration | Option 2 at quarter-year, monthly, half-month and weekly steps, force of infection at mid-step | Exp 30's table | Exp 30's exact-aging reference; weekly option 2; and option 2's own limit, extrapolated from its half-month and weekly runs |
| X0 ("start of step") | Is the mid-step force of infection worth its cost? It doubles the work per step, so drop it if the simpler version is as accurate at the chosen step | Option 2 at the same steps, force of infection at the start of each step | Exp 30's table | The same as X |
| S ("start") | Can `fastode_age` make the starting table itself, as the calibration will for each draw, and how short can the frozen run be? | `fastode_age`'s frozen run (option 2, the step chosen from X) for 200, 100 and 50 years, then 1950–2041 from the table it makes | Its own table | Exp 30's table and its runs from it |
| P ("prior") | Does the chosen step hold across the prior, not only at two posterior sets? IMIS starts from prior draws, including fast rates (`k_latent` up to 12 per year) and high `beta`, where a coarse step does worst | The whole draw as the calibration would run it (the 200-year frozen run, then 1950–2041), all at one step: quarter-year, monthly and weekly, at 20 prior draws with an untreated endemic equilibrium without ages (exp 27's starting point) | Its own table (arm S's method) | Weekly (no Starsim runs exist at these draws) |

- **Exp 30's exact-aging reference:** 2 × eighth-year − quarter-year, from Starsim's weekly ODE runs. It stands in for exact aging (exps 29–30), but carries Starsim's weekly-step error and the extrapolation's own error. Exp 30 found eighth-year groups 0.2–0.3% from it.
- **What arm C should show:** without treatment (1950–94), the copy should be closer to Starsim's weekly run than to its monthly run, since both Starsim runs approach the same continuous ODE as their step shrinks. From 1995, the 12 Erlang stages that stand in for the fixed course add up to about 0.1%.

### Measured

- **Exp 30's yearly outputs:** prevalence (asymptomatic + symptomatic) for all ages and in five age groups (0–14, 15–29, 30–49, 50–64, 65+), the share asymptomatic, and per person: incidence, notifications, the number on treatment, TB mortality, recent and remote latent infection, and the cleared share. Also population and the age groups' shares.
- **The calibration outputs, as exp 24 defines them:**
  - point values on 1 January: prevalence and the share asymptomatic in 2020, prevalence in 2040;
  - totals over calendar years: incidence, notifications × (1 + `over_detection`), and TB deaths;
  - the 10 targets' values;
  - prevalence at 15+, the survey's age range, as well as all ages.
- **Ratios:** each as model ÷ reference, by year and over exp 30's windows (1950–94, 1995–2024, 2025–40). All runs are deterministic, so there's no SE.
- **Against agents, reported only:** arms X and X0 at each step ÷ exp 30's saved agent runs for sets 0 and 6 (`runs.obj`; same starting tables; monthly and weekly steps, 20 reps for set 0 and 60 for set 6), by window, mean ± SE across reps. No new agent runs.
  - **As a sanity link:** `fastode_age` ÷ agents should look like exp 30's Starsim ODE ÷ agents. If it doesn't, something in the new model is off that the comparisons between ODEs missed.
  - **As a possible explanation:** X0 at monthly steps takes the force of infection at the start of each step, as monthly agents do. If it reproduces monthly agents' offset of about 1% from the ODE, that explains the slow settling found in exp 27 (README, "Known unexplained differences in v4").
  - **Not used to choose the step:** agents at monthly steps carry their own step artefacts of about 1% (exps 27, 30). The calibrated model should be the continuous-time one, which agents match at weekly steps.
- **Run time** per run on one laptop core, for each option and step, with and without the start, and the projected cost of a 310,000-draw IMIS run on 4 cores.

## What would answer the question

- **The copy is right** if arm C matches Starsim's weekly quarter-year ODE within about 0.1% in every output and age group before 1995, closer than Starsim's monthly run does, and within about 0.2% after 1995.
- **Option 2 is right** if, as the step shrinks, its outputs converge: differences fall by about 4× when h halves with the midpoint force of infection, and about 2× without it. The limit should be within about 0.3% of exp 30's exact-aging reference in population, prevalence (all ages and by age group) and ratios per prevalent case.
- **The step for calibration:** the coarsest step whose calibration outputs and yearly burden are within 0.2% of weekly option 2 at both sets and at every prior draw in arm P. 0.2% is about the size of the monthly-to-weekly gap in Starsim's ODE, and well under the targets' σ.
- **The start:** the shortest burn-in whose targets are within 0.05% of the 200-year start's.
- **Fast enough** if a draw (start plus 1950–2041) takes about 1 s or less on one core, about a day per IMIS run on the laptop.
- **If option 2 can't match exp 30's reference** while arm C matches Starsim, the difference is in the aging or in the step scheme. Compare option 2's limit with option 1 at eighth-year groups, extrapolated as exp 30 did.

## Next

The recalibration with demography, as its own experiment, using `fastode_age` at the chosen step. It will follow exp 24's IMIS settings and checks ([CONVENTIONS.md 8](../../CONVENTIONS.md#8-calibrating-the-ode-a-standalone-copy-at-equilibrium-with-imis)):

1. a coverage check;
2. re-identification;
3. the real targets;
4. a second seed, compared on parameter marginals as well as predictions.

Then agents against the calibrated ODE, as exps 21 and 27 did: v4's agents with exp 30's demography and start, at about 12 sets drawn from the posterior (from both seeds), at monthly and weekly steps, ÷ `fastode_age` at the same sets. This checks that agents run at the calibrated parameters reproduce the fit.

`v4/tests/test_fastode_age.py` keeps the check against v4's ODE, as `test_fastode.py` does for `fastode`: set 0 against a small saved reference from exp 30, at the tolerances found here. Rerun it after bumping the v4 pins.

## Not covered

- **The calibration** itself (next).
- **Child transmission, child-specific rates and lower coverage for children:** rates stay the same at every age.
- **WPP's sex ratio at birth, and migration.**
- **Networks:** mixing is homogeneous.

## Files

- `v4/tbode4/fastode_age.py`: the model (both options), its start and its calibration outputs.
- `run.py`: `copy` (arm C), `steps` (arms X and X0), `start` (arm S), `prior` (arm P), `summary` (tables, the agents comparison and `steps.png`; run times are recorded with each run). Each saves its raw outputs (`copy.obj`, `steps.obj`, `start.obj`, `prior.obj`) with `sc.save` and writes `versions.json`.
- `SUMMARY.md`: findings, written after running.
- `v4/tests/test_fastode_age.py`: option 1 against v4's ODE run live, option 2's convergence, and the starting table's round trip.
- Uses exp 30's `ode.obj`, `start.obj` and `runs.obj`, exp 27's sets, exp 24's priors (`calib.py`), and `tbode4.demog`.
