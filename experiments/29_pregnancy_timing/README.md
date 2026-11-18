# 29: Do agents' births match the ODE's when `ss.Pregnancy` uses fertility by age at delivery?

**Status:** run (findings in `SUMMARY.md`).
**Track:** v4 (`v4/.venv`; tbsim `origin/v4` @ `619c296`, starsim `origin/v4` @ `d83f4f87`).

## Why this experiment

Exp 28 compared v4's agents with v4's ODE using India's fertility and mortality from the UN World Population Prospects (WPP) 2024. Agents' population grew more slowly: 0.92 of the ODE's by 2000, with the ODE split by single years of age (D3 ÷ D2 in exp 28).

The likely cause is how `ss.Pregnancy` uses fertility rates.

- **WPP's age-specific fertility rate** is births per woman per year by the mother's age *at delivery*.
- **The ODE uses it that way.** `ss.birth(rate=ss.by(age, sex))` gives each woman births at her age group's rate.
- **`ss.Pregnancy` applies it at the woman's age *at conception*,** rescaled to women who aren't already pregnant. It delivers 32–42 weeks later (0.769 years on average), taking the rate of the year of delivery but not the age.

So agents' deliveries come about 0.77 years later in a mother's life than WPP's rates say. That lengthens each generation and slows growth. In exp 28, agents' realized fertility by age at delivery was 0.80–0.93 of the input at ages 15–24 and 1.0–2.8 of it at ages 35–49, while total fertility was 0.99–1.00. A rough estimate put this at about two thirds of the gap. Exp 28 didn't test it.

This experiment tests the fix: give `ss.Pregnancy` a fertility table shifted younger by the mean gestation, so a woman who conceives at age x takes WPP's rate for age x + 0.77. It also measures what remains, and checks the mechanism from the ODE's side. The results are meant to support a proposal to Starsim about how `ss.Pregnancy` reads fertility rates.

## Background

- **Exp 28's settings carry over:**
  - WPP 2024 rates as bundled with Starsim (every 10 years), interpolated linearly by year;
  - fertility at ages 10–54 in 5-year groups;
  - deaths by age and sex (`ss.death(rate=ss.by(age, sex))`) in both solvers;
  - WPP's 1950 population by age and sex at the start;
  - 1950–2050, monthly steps;
  - no TB.
- **Single years of age in the ODE.** Exp 28 found that WPP's 5-year groups distort the ODE's aging. Single years of age (0–99, 100+) fix most of it at about the same run time, so the reference ODE here uses them. The rates stay on WPP's groups.
- **Feasibility, checked before this design** (one rep, 50,000 agents, 1950–1970):
  - `ss.Pregnancy` accepts a fertility table on a quarter-year age grid.
  - With the shifted table, realized fertility by age at delivery over the 1960s is 0.96–1.02 of WPP's at ages 15–44, against 0.82–1.27 unshifted. Total fertility is 1.00 in both.
  - Population in 1970 ÷ 1950 is 1.609 shifted, against 1.592 unshifted.
  - Ages 45–54 are still above WPP's rate (1.12 and 1.29). Ages 10–14 are well below it (0.26), because nobody conceives before 10 − 0.77 years.
- **Other differences between `ss.Pregnancy` and the ODE's births,** which the shift doesn't remove:
  - **Pregnancy loss when the mother dies.** A fetus is lost if its mother dies while pregnant. At India's adult female mortality in the 1950s (about 1–2% a year), that's roughly 1% of births.
  - **The start.** With `burnin=True` (the default), women are pregnant at the start, at the gestational ages they'd have reached. In exp 28, agents' births in 1950 were 12% below the ODE's crude birth rate. `burnin_plan()` makes each woman pregnant with probability 1 − e^(−h·d), with h her fertility rate and d the mean gestation. A steady state has a share h·d pregnant, so this seeds too few pregnancies: about 10% too few at the peak rate (0.29 per year), consistent with a check of one run (3,189 pregnant at the start, against about 3,530 expected).
  - **Whole steps.** At monthly steps a pregnancy lasts a whole number of months.
- **Starsim's warning** that "pregnancy makes births, and so does the sim's ss.Births or ss.birth()" appears in every agent run, in exp 28 too. It's a false positive: the flow that holds `ss.death` has an `n_born` attribute. In a check, births equalled deliveries. Every run here checks that the population changes by births minus deaths at each step.

## Question

1. With fertility shifted younger by the mean gestation, does `ss.Pregnancy` give the same population size, age distribution and births as the ODE?
2. What remains, and how much comes from each of pregnancy loss when the mother dies, the start and whole steps?
3. From the ODE's side: does shifting the ODE's fertility *older* by the mean gestation reproduce exp 28's slower growth in agents?

## Design

| Arm | Solver | Fertility | Steps | What it tests |
|---|---|---|---|---|
| R | ODE, single years of age | WPP's rates by age at delivery (exp 28's D2) | Monthly | The reference |
| P0 | Agents | `ss.Pregnancy` as in exp 28's D3: WPP's rates at the age at conception, ages 10–54 | Monthly | Exp 28's gap, rerun here with new seeds |
| P1 | Agents | `ss.Pregnancy` with the shifted table: at conception age x, WPP's rate for age x + g, with g = 0.769 years (its mean gestation); quarter-year age grid; conception from age 10 − g to 55 − g | Monthly | The fix |
| P1w | Agents | As P1 | Weekly | Whole-month pregnancies at monthly steps |
| O1 | ODE, single years of age | WPP's rates shifted older: at age a, the rate for age a − g | Monthly | The mechanism from the ODE's side: births later in life without pregnancy states |
| R2, R4, R8 | ODE, half-, quarter- and eighth-year age groups | As R | Monthly | Whether R's aging has converged |
| P2 | Agents | As P1, with the steady-state share pregnant at the simulation's start: h·d (conception rate × mean gestation), overriding `burnin_plan()`'s 1 − e^(−h·d) | Monthly | Fix 2 |
| P3 | Agents | As P2, with each rate raised by e^(m·g), m women's death rate at age x + g/2 | Monthly | Offsetting pregnancies lost when the mother dies |
| P2s, P3s | Agents | As P2 and P3, starting with WPP's 1950 ages and half of each age women | Monthly | Whether the ODE's starting sex split explains the rest |
| R4x, R8x, P2b | ODE (quarter- and eighth-year ages) and agents | The ODE started from WPP's sex split within each age group; agents as P2, with births as female as the ODE's (48.63%) | Monthly | The fix to the ODE's start |

- **The arms were added in rounds,** each to explain what the previous ones left:
  - **R to O1:** the original design. P1 was still 4% below R.
  - **R2 and R4:** exp 28 had shown that the width of the ODE's age groups changes its population (by up to 7%, 5-year groups against single years), so part of P1 ÷ R could be R's own aging error.
  - **P2, P3 and R8:** the start and pregnancy losses were then estimated to account for about 1% of what was left; P2 and P3 test each, and R8 checks the extrapolation from R2 and R4 to exact aging.
  - **P3s and P2s:** P3's births per woman were right, but its births per person were 2% below the ODE's. That led to the ODE's starting sex split: Starsim's compartmental solver starts with half of each age women, where WPP's 1950 population for India has 48.6%.
  - **R4x, R8x and P2b:** the fix, a joint split by age and sex in `sim.stratify()`. The plan was to test it with births 50% girls and with an illustrative 48.5%, but the ODE's births then follow its starting female share (48.63%) whatever `ss.birth`'s `p_female` is, so only that version could be built.

  Earlier arms keep their seeds; `run.py more` runs the latest arms and merges them into the saved runs. The first round also found a bug in exp 28's realized-fertility calculation, which counted girls under 10 among women aged 10–14 (fixed); every arm then was rerun with the same seeds.
- **R8 keeps every step's results by group:** Starsim keeps them at every step only up to 100 MB per run (`Engine.group_bytes`), raised here.
- **Agents:** 100,000 at the start, growing to about 500,000 by 2050; 10 reps per arm, with new seeds.
- **Deaths:** exp 28's `ss.death` by age and sex in every arm. Pregnancy's own maternal and infant death rates stay at their default of 0.
- **The shifted table:**
  - built from the same interpolated WPP rates;
  - each quarter-year age takes the 5-year group's rate at age + g (a step function, as WPP's groups are);
  - one row per calendar year (`ss.Pregnancy` takes the year of delivery itself).
- **Size:** 80 agent runs, about 23 s each monthly and 90 s weekly; 2–3 s for each ODE run. Run on the laptop, on 4 cores.

## Measured

- **Population:**
  - size ÷ its 1950 value every 10 years;
  - the share in age groups 0–14, 15–29, 30–49, 50–64 and 65+;
  - for agents, mean ± SE across reps;
  - ratios between arms: each agent arm against R and against the finer ODEs, each fix against the arm before it, and the finer ODEs against R;
  - exact aging, extrapolated from the finer ODEs.
- **Births:**
  - crude birth and death rates per 1,000 by year;
  - realized fertility by mother's age at delivery (5-year groups) ÷ WPP's rate, and total fertility, by decade;
  - mean age at childbearing, realized against WPP's.
- **What remains in P1:**
  - **Pregnancy loss when the mother dies:** pregnancies ended this way per 1,000 births. `ss.Pregnancy`'s miscarriage and stillbirth results don't count these losses: they read 0 even at a death rate of 0.2 per year (checked before running). So they're counted by balance in each calendar year: women pregnant at its start, plus conceptions, minus deliveries, minus women pregnant at its end. Abortions and spontaneous losses are 0 here. In the check, the balance was exactly 0 without deaths, and 454 at 0.2 per year against about 488 expected.
  - **The start:** births in 1950 and 1951, P1 ÷ R.
  - **Whole steps:** P1w ÷ P1.
- **The balance check,** every agent run: population change against births minus deaths, every step.

## What would answer the question

- **The fix works** if P1 is within about 1% of R in population size through 2050, and within 1% of R's share in each age group. A possible exception is the first few years, if they're explained by the start.
- **The mechanism is confirmed from both sides** if O1 is within about 1% of P0 in population size.
- **If P1 is still more than 1% from R,** the measurements of what remains say which part is responsible. Pregnancy loss when the mother dies is real biology, so the proposal would leave it and describe it. The start and whole steps are things the proposal could address.
- **The proposal to Starsim** follows from the results. If P1 matches R, it asks for an option in `ss.Pregnancy` to read age-specific fertility at the expected age at delivery, or for its documentation to say that WPP-style rates are by age at delivery and need shifting. It would also report the false-positive warning. Drafted after this experiment, in its own folder under `proposals/`.

## Not covered

- **TB.** Once births match, the TB comparison comes with the next experiment: a starting state with ages for both solvers (the ODE's steady state at 1950's demography), on the consistency table's demography row.
- **tb_LAI_TPT's own setup** (Starsim's defaults, ages 15 and over; exp 28's D4). Unchanged here.
- **Migration,** which neither model includes.

## Files

- `run.py`: runs the arms and the summaries, with every run's raw outputs saved with `sc.save`, and `versions.json`.
- `SUMMARY.md`: findings, written after running.
