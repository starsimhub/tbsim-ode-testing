# 28: Do v4's agents match its ODE with India's age-specific fertility and mortality?

**Status:** run (findings in `SUMMARY.md`).
**Track:** v4 (`v4/.venv`; tbsim `origin/v4` @ `619c296`, starsim `origin/v4` @ `d83f4f87`).

## Why this experiment

This is the first step on the consistency table's "Age-specific fertility/mortality (WPP)" row. tb_LAI_TPT runs its ABM from 1950 with a growing population: births from India's age-specific fertility rates (`ss.Pregnancy`), deaths from India's mortality rates by age and sex (`ss.Deaths`), and an initial age distribution, all from the UN World Population Prospects (WPP) 2024. Every v4 experiment so far has used a constant population, with births replacing every death at μ = 1/70 per year and no ages. An ODE with age structure has to approximate aging by moving people between age groups, so this is a feature where the two solvers may not agree exactly.

The TB model is exp 27's: treatment, falling TB mortality and the calibrated natural history, which agents already match without demography. Demography is the one change.

## Background

- **Demography data in Starsim v4:** `starsim_data` 2.2.0 bundles WPP 2024 (Medium variant) for every country: crude birth rates, mortality (`mx`) by sex and 5-year age group (with 0 and 1–4 separate), age-specific fertility by 5-year group of mother (10–14 to 50–54), and population by age and sex. These come every ten years (1950, 1960, …, 2100). `ss.population('India', year=1950)` gives the starting age and sex distribution.
- **tb_LAI_TPT's data** are the same WPP 2024 Medium estimates, pulled from the WPP files (`pull_demographics.py`): annual mortality by sex and the same age groups, and annual fertility by single year of the mother's age (15–49). The ODE won't use finer groups than 5 years (user, 2026-11-12), so the bundled data are enough if the interpolation between decades doesn't matter. The input check below tests that.
- **What each solver can run** (checked with small test runs before this design):
  - **v4's ODE refuses `ss.Births` and `ss.Deaths`,** which add and remove agents. Its form of demography is `ss.birth(rate=...)` and `ss.death(rate=...)` in `demographics=`, with rates by age and sex given as `ss.by(...)` of `ss.timeseries`, in a sim stratified by age and sex (`sim.stratify(age=..., sex=...)`).
  - **`ss.Pregnancy` in the ODE is too large with TB:** its compartmental form tracks pregnancy and TB states jointly, about 20,000 compartments. Births at each woman's age-specific fertility rate, `ss.birth(rate=ss.by(age, sex))`, give the same fertility without pregnancy states. With exp 27's treated model, a 1950–2040 ODE run takes 4 s.
  - **Agents can't take births by age and sex** (`ss.birth(rate=ss.by(age, sex))` is for compartmental runs only), so they need `ss.Pregnancy`. They do take the ODE's deaths, `ss.death(rate=ss.by(age, sex))`.
  - **`ss.Pregnancy`'s defaults differ from the ODE's births** in two ways: it ignores fertility under 15 (`min_age` = 15), and with a country's bundled rates it uses the nearest decade's, not an interpolation. Losses, abortion, infertility and maternal and infant deaths are all 0 by default. Given the same rates interpolated by year, and ages 10–54, the only remaining difference is gestation: a birth comes about 9 months after conception.
  - **First test runs** (one rep each, 50,000 agents, no TB), population ÷ 1950: in 2000, the ODE 3.19, agents with Pregnancy set to match 2.99, with its defaults 2.91, WPP 3.05; in 2050, 5.02, 4.98, 4.78 and 4.85. Stage 1 is designed to explain these gaps.
- **How the ODE ages people:** it moves people out of each age group at a constant rate, set once at the start from the 1950 pyramid. With groups no wider than WPP's 5-year groups, that rate is 1 ÷ the group's width (1 per year for age 0, 1/4 for 1–4, 1/5 for each 5-year group). Time in a group is then exponential: the mean stay is right, but a cohort's ages spread (the 1950 birth cohort's model ages in 2000 have mean 50 but an SD of about 16 years, instead of all being 50). Agents age exactly. Mortality rises steeply with age, so the spread may change the population's size and age distribution.
- **The TB model:** exp 27's, `tbode4.TBTreat`: v4's `tbsim.TB` with the two-state latent model, every reinfection risk at π, homogeneous mixing, a fixed 6-month treatment course, and coverage, efficacy and `sym_dead` over time as calibrated in exp 24 (`ss.timeseries`, linear between exp 24's knots). Treatment starts at the end of 1994. Exp 24 calibrated it with births replacing every death at μ = 1/70; that calibration isn't redone, so with India's younger, growing population the fit to the targets will move. This experiment compares agents with the ODE, not either with data.
- **Earlier agents ÷ ODE results without demography (exp 27):** across 12 sets from exp 24's posterior, prevalence 0.997–1.008 by window and step, and ratios per prevalent case within 0.3%. Two monthly-step effects of about 1%: agents apply a time-varying rate as if half a step earlier than the ODE (treatment starts 9% low in 1995, within 1% by 2000), and monthly agents settle slowly from the ODE's starting state.

## Question

With India's WPP fertility and mortality, an initial age distribution and a growing population, do v4's agents reproduce v4's ODE?

1. **Demography alone:** population size, age distribution, crude birth and death rates, against each other and against WPP.
2. **With TB and treatment:** prevalence (all ages and by age group), the share asymptomatic, and incidence, notifications and TB mortality per prevalent case.

Two related questions for the ODE alone: how long does the model take to settle (burn in) from an age-independent start once it has demography, and how much does demography change the burden?

## Design

**Age groups:** WPP's mortality groups, 0, 1–4, 5–9, …, 95–99 and 100+ (22 groups), by sex, so that mortality and fertility need no averaging. Fertility applies to women aged 10–54.

**Demography inputs, for both solvers:** the bundled rates, interpolated linearly by year between decades.

- **ODE:** `ss.birth(rate=ss.by(age, sex))` at each age group's fertility rate (0 for men), and `ss.death(rate=ss.by(age, sex))` at `mx`.
- **Agents:** `ss.Pregnancy` with the same fertility rates as an annual table, `min_age` = 10 and `max_age` = 55, and the same `ss.death`.

**Input check, before running:** the bundled rates, interpolated, against tb_LAI_TPT's annual files (fertility grouped to 5 years), for each year 1950–2050. Reported as total fertility and life expectancy at birth by year, and the largest relative difference in any age group's rate.

### Stage 1: demography alone

No TB. 1950–2050, monthly steps, starting from WPP's 1950 population by age and sex.

| Arm | Solver | Fertility | Mortality | What it tests |
|---|---|---|---|---|
| D1 | ODE | `ss.birth` by age and sex | `ss.death` by age and sex | The reference |
| D2 | ODE, single years of age | As D1 | As D1 | How much of the gap is from aging between 5-year groups |
| D3 | Agents | `ss.Pregnancy`, set to match D1 | As D1 | The same model as D1, as agents, but for gestation |
| D4 | Agents | `ss.Pregnancy('India')`, defaults | `ss.Deaths('India')` | Starsim's built-in modules, as tb_LAI_TPT would call them |
| D5 | Agents | `ss.Pregnancy` with tb_LAI_TPT's annual file | `ss.Deaths` with tb_LAI_TPT's annual file | tb_LAI_TPT's demography, run in v4 |

- **Agents:** 100,000 at the start, 10 reps each, new seeds.
- **D3 against D1** is the consistency check. **D4 against D3** shows the effect of Pregnancy's defaults. **D5 against D4** is the data source.
- **D1's ages of the 1950 birth cohort** over time show the spread from aging between groups; D2 shows how much it matters.

### Stage 2: TB with demography

- **The model:** exp 27's, at 2 of its 12 sets (set 0 from exp 24's first IMIS seed and set 6 from the second, one from each stretch of the posterior's ridge), with `mu` = 0 and stage 1's demography (D1 in the ODE, D3 as agents). TB deaths are not replaced.
- **Start:** 1950, with everyone's TB state drawn from the set's untreated equilibrium at μ = 1/70 without ages (`fastode`, as in exp 27), in the same shares in every age and sex group, for both solvers. This isn't an equilibrium with ages, but both solvers start from it, so agents ÷ ODE is still a fair comparison. Children are infected as much as adults at the start; it takes decades for this to wash out of the young age groups. Everything but demography is constant until treatment starts at the end of 1994. The run ends in 2041.
- **The ODE:** v4's ODE, at monthly and weekly steps. First-entry incidence from `tbode4.ode_incidence`.
- **Agents:** at monthly and weekly steps, 20 reps each, new seeds. 100,000 agents at the start, growing to about 440,000 by 2040, when prevalence is about 240/100k, about 1,000 people with TB. Incidence and the pitfall checks are counted with `tbode4.Incidence` and `tbode4.TxChecks`.
- **Size:** 2 sets × 2 steps × 20 reps = 80 agent runs. Timed on one run first; expected under an hour on 4 cores.
- **Over-detection:** applied to both sides' notifications in the same way, so it cancels in agents ÷ ODE.
- **Flows:** dated by Starsim's current indexing, where index i is the step from t_i to t_(i+1) (exp 26).

### Stage 2, ODE only: burn-in and the effect of demography

Cheap arms (seconds each) for how long the model takes to settle with demography, and how much demography changes the burden. ODE at monthly steps, both sets, the same starting state as the main arm.

| Arm | Demography | Treatment | Years | What it shows |
|---|---|---|---|---|
| Main (above) | WPP, by year | From 1995 | 1950–2041 | The reference for the others |
| B1 | WPP, by year | None | 1950–2041 | The untreated course with demography: settling from the start and demographic change together |
| B2 | Frozen at 1950's rates | None | 1950–2150 | Burn-in alone: this model has a steady state, so it shows how long the age-independent start takes to settle |
| B3 | None: births replace deaths at μ = 1/70, no ages (exp 27's ODE) | From 1995 | 1950–2041 | How demography changes the burden's level and trend |
| S | WPP, by year, with single years of age (0–99, 100+) | From 1995 | 1950–2041 | How much of agents ÷ ODE by age comes from aging between 5-year groups (stage 1's D2, with TB) |

- **S** is the main arm with the population split by single year of age, at both sets and at monthly and weekly steps; the rates stay on WPP's groups, so each single year in a group has its group's rate. It's compared with the agents as the main arm is, and with the main arm. Starsim keeps results by group at every step only while they take at most 100 MB in a run (`Engine.group_bytes`), and every few steps otherwise; S at weekly steps takes more, so the limit is raised for S.
- **Untreated (B1, B2)** means no treatment and TB mortality (`sym_dead`) held at its value before 1995, so that only demography moves the model.
- **B2's steady state** includes the population's own: with 1950's rates held fixed, the age distribution moves from WPP's 1950 pyramid towards the stable one those rates imply. Its population keeps growing, but every share settles. The shares at 2150 are taken as the steady state; the run is checked to be flat by then.
- **B1 against B2** separates settling from the start, which both share, from demographic change, which only B1 has.
- **B3 against the main arm** is the effect of demography on the burden at these parameters. Its size will change when the model is recalibrated with demography; how long burn-in takes depends mainly on the model's slow processes (remote latent infection, aging, demographic change), so B2 should carry over.

### Measured

- **Stage 1:**
  - population size ÷ its 1950 value, and its share in age groups 0–14, 15–29, 30–49, 50–64 and 65+, every 10 years, against WPP's population;
  - crude birth and death rates per 1,000 by year, against WPP's;
  - for each arm, mean ± SE across reps, and D3 ÷ D1, D4 ÷ D3 and D5 ÷ D4.
- **Stage 2, agents ÷ ODE:**
  - active prevalence (asymptomatic + symptomatic) per person, all ages and in the same five age groups, and the share asymptomatic;
  - per prevalent case: incidence (first entry to asymptomatic TB per episode, plus relapses), notifications, the number on treatment, and TB mortality (from symptomatic TB and on failing courses);
  - recent and remote latent infection, and the cleared share;
  - the population's size and age shares, as in stage 1.

  Each as a yearly mean, then means over three windows: 1950–1994 (before treatment, settling from the age-independent start), 1995–2024 (the scale-up) and 2025–2040. Mean ± SE across the 20 reps, for each set and step.
- **Stage 2, ODE-only arms:**
  - **burn-in time** in B2: the first year after which each output stays within 1% of its steady state for good, for prevalence (all ages and each of the five age groups), recent and remote latent infection, incidence and TB mortality per 100,000, and the five age groups' population shares;
  - **B1 ÷ B2** over time, for the same outputs;
  - **the main arm ÷ B3** over time, for prevalence, incidence, TB mortality and notifications per 100,000, and by window.
- **Expected:**
  - **Recent latent infection reads high** in monthly end-of-step snapshots (×1.2–1.7 in exp 21); it's reported, not counted as a mismatch.
  - **Exp 27's two monthly-step effects** (above): treatment starts about 9% low in 1995, and monthly agents drift by about 1% from the ODE's starting state.
  - **Stage 1's gap in population,** carried into TB counts. Comparing per person and per prevalent case keeps it out of the TB ratios as far as it can be; the age groups show what's left.

**Figures:**

- **Stage 1:**
  - **Population size and age distribution over time:** population ÷ 1950, the five age groups' shares by year, and age pyramids (share by 5-year group and sex) in 1950, 1980, 2000, 2020 and 2050, for each arm against WPP;
  - **Crude birth and death rates** by year, for each arm against WPP.
- **Stage 2,** for each set, agents' reps, their mean and the ODE, at monthly and weekly steps, 1950–2040:
  - **Burden over time:** prevalence, incidence, TB mortality and notifications per 100,000 per year;
  - **Prevalence by age group** over time, per 100,000 in each of the five groups;
  - **Age distribution** of the TB runs' population over time, as in stage 1, so stage 1's gap can be followed into the TB results;
  - **Agents ÷ ODE** by output and window, with SE.
- **Stage 2, ODE-only arms,** for each set:
  - **Burn-in:** prevalence (all ages and by age group) and latent infection over time in B1 and B2, with B2's steady state and its burn-in times marked;
  - **Effect of demography:** prevalence, incidence, TB mortality and notifications per 100,000 over time, the main arm against B3.

The raw outputs behind each figure are saved with `sc.save`: per-rep yearly values for each arm and stage, and the ODE runs.

## What would answer the question

- **Demography matches** if D3 is within about 1% of D1 in population size and every age share through 2050. D1 against WPP shows what the ODE leaves out (migration, interpolation between decades); it doesn't need to be exact.
- **TB matches** if agents ÷ ODE is within about 2 SE of 1, and within about 1%, for prevalence and each ratio per prevalent case, in the two later windows at weekly steps, as in exp 27.
- **If demography doesn't match,** D2 shows whether the age groups are the cause. If they are, the row is ⚠️ for the ODE at 5-year groups, with the size of the gap. Gestation, the other difference between D1 and D3, should matter little: it delays births by about 9 months.
- **If the TB comparison doesn't match but demography does,** compare prevalence by age group to find where it comes from. Aging spreads a cohort's ages in the ODE, which mixes children infected at the start into adulthood sooner.
- **Burn-in:** B2's burn-in times say how early a run with demography must start, from an age-independent state, for each output to settle; if they are longer than the run before the targets (1950 to 2015–2020), a starting state with ages is needed. B1 ÷ B2 shows how much of B1's movement before 1995 is demographic change rather than settling.
- **D4 and D5 show what tb_LAI_TPT's way of calling the modules would change.** If D4 differs from D3 by more than about 1% of population, the reason goes in SUMMARY.md.

## Not covered

- **Treatment by age:** tb_LAI_TPT treats under-15s at a lower coverage (`rel_tx_covg_u15`). Here coverage is the same at every age.
- **Child transmission** (no transmission under 10, the table's next row), other age-specific TB rates and age mixing.
- **Migration:** WPP gives net migration, which tb_LAI_TPT doesn't use.
- **Recalibration:** with demography, the sets' fits to exp 24's targets will move.
- **An equilibrium start with ages:** both solvers start from the same age-independent state instead.

## Files

- `run.py`:
  - `inputs` runs the input check → `inputs.csv`, `inputs.log`, `inputs.png`;
  - `demog` runs stage 1 → `demog.obj`, `demog.csv`, `demog_population.png` (size, age shares and pyramids), `demog_rates.png`, `demog.log`;
  - `tb` runs stage 2 → `ode.obj`, `runs.obj`, `by_window.csv`, `pitfall_checks.csv`, `tb_burden_<set>.png` (prevalence, incidence, mortality, notifications), `tb_prev_by_age_<set>.png`, `tb_population_<set>.png`, `agents_vs_ode.png`, `tb.log`;
  - `burnin` runs the ODE-only arms → `burnin.obj`, `burnin_times.csv`, `burnin_<set>.png`, `demog_effect_<set>.png`, `burnin.log`;
  - `single` runs arm S (after `tb`) → `ode_single.obj`, `by_window_single.csv` (agents ÷ S), `single_vs_groups.csv` (S ÷ the main arm's ODE), `single_by_age_<set>.png`, `agents_vs_ode_single.png`, `single.log`, `single.out`;
  - `summary` redraws the summaries and figures from the saved runs;
  - every run writes `versions.json`.
- Uses exp 27's `run.py` (its sets, `sets.obj`, and its starting state and treatment schedules), `tbode4`'s `TBTreat`, `Incidence`, `TxChecks`, `ode_incidence` and `fastode`, and tb_LAI_TPT's `data/fertility_rates_India.csv` and `data/death_rates_India.csv` (read only).
- `SUMMARY.md`: findings, written after running.
