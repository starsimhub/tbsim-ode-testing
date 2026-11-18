# Proposal 06: births and the starting population, agents against the compartmental solver

**For:** Starsim v4 (`ss.Pregnancy` in `starsim/demographics.py`; the compartmental solver's starting population in `starsim/sim.py`).
**Status:** draft, not sent.
**Track:** starsim `origin/v4` @ `d83f4f87` (4.0.0rc4), tbsim `origin/v4` @ `619c296`.
**Evidence:** [exp 29](../../experiments/29_pregnancy_timing/SUMMARY.md) (agents with `ss.Pregnancy` against the compartmental solver's births); [exp 28](../../experiments/28_age_demography/SUMMARY.md) (where the gap was found).

## Summary

We compared Starsim's agents with its compartmental (ODE) solver, using India's fertility and mortality from the UN World Population Prospects (WPP) 2024, 1950–2050.

- **The models:**
  - **The compartmental solver** takes births as `ss.birth(rate=ss.by(age, sex))`.
  - **Agents** can't take that, so they use `ss.Pregnancy` with the same age-specific fertility rates.
  - **Deaths** are the same in both, `ss.death(rate=ss.by(age, sex))`.
- **The gap:** agents' population was 7.5% smaller by 2000.
- **What explains it:** all of it, by three things in Starsim and one in how finely the compartmental model is split by age. With all four dealt with, the two match within 0.1% through 2050. Change 5 also covers how the compartmental solver sets the sex of births, which it ties to the starting population.

| Part of the gap, India, 2000 | Size | Proposed change |
|---|---|---|
| `ss.Pregnancy` applies fertility rates at the age at conception, not delivery | 3.3% of population (5.5% by 2050) | 1 |
| `ss.Pregnancy` seeds too few pregnancies at the start (`burnin=True`) | 0.5% | 2 |
| The compartmental solver starts with half of each age female; WPP's 1950 population has 48.6% | 1.9% (2.1% by 2050) | 5 (a workaround exists) |
| Single-year age groups don't age people exactly in the compartmental solver | 2.2% | None: use finer groups |

Two more changes affect what users see, not the dynamics:

- **Change 3:** `ss.Pregnancy`'s loss results don't count pregnancies lost when the mother dies.
- **Change 4:** a warning about double-counted births is a false positive.

## The setting

- **Model:** demography only. Rates are WPP 2024 as bundled with Starsim (every 10 years), interpolated linearly by year:
  - deaths by age and sex, in both solvers;
  - fertility at ages 10–54 in 5-year groups;
  - India's 1950 population by age and sex at the start.
- **Agents:** 10 reps of 100,000 at the start (about 500,000 by 2050), monthly steps.
- **Compartmental runs:** age groups of a year, half, quarter and eighth of a year. The population changes by half as much at each halving of the width, so we extrapolate to exact aging.

## 1. Age-specific fertility at the age at delivery

**What happens now.** `conception_rate()` looks up the rate for the woman's current age, and for the year in which a pregnancy started now would end (`year = self.t.now('year') + self.mean_dur_pregnancy`). The year is shifted to delivery, but the age isn't. So a rate for mothers aged 20–24 at delivery is applied to women conceiving at 20–24, who deliver at 20.8–25.8. WPP's age-specific fertility rates, like most published ones, are by the mother's age at delivery.

**What it does:**

| India, WPP 2024, agents (10 reps) | As now | Rates shifted to age at delivery |
|---|---|---|
| Realized fertility ÷ WPP's rate, by age at delivery, ages 15–24 | 0.84–0.93 | 0.99–1.04 |
| Realized fertility ÷ WPP's rate, ages 35–49 | 1.08–2.71 | 0.98–1.03 |
| Mean age at childbearing, against WPP's | 0.77–0.81 years older | within 0.04 years |
| Total fertility ÷ WPP's | 0.99–1.00 | 1.00 |
| Population in 2000 | 1 | 1.034 (1.055 by 2050) |

The compartmental solver confirms the size of the effect from the other side. Shifting its fertility 0.77 years *older* lowers its population by 3.2% by 2000 and 5.0% by 2050.

**Proposal.** In `conception_rate()`, look the rate up at the expected age at delivery as well as the year:

```python
a = (age if idx is None else age[idx]) + self.mean_dur_pregnancy  # Rates by age at delivery, as WPP's are
```

Our workaround gives the same result: a table on a quarter-year age grid, in which a woman conceiving at age x gets the rate for age x + 0.769 (`shifted_table()` in exp 29's `run.py`).

**Choices for the maintainers:**

- **An option or the default.** Something like `fertility_age='delivery'` (rates by age at delivery, shifted) or `'conception'` (as now). We'd suggest `'delivery'` as the default for data tables, since published age-specific fertility is by age at delivery, but it changes existing results.
- **`min_age` and `max_age`.** These are ages at conception now. With rates by age at delivery, users giving WPP's groups (10–54 or 15–49) would expect them to be ages at delivery too, i.e. conception from `min_age` − d to `max_age` − d.
- **Documentation either way.** The docstring could say which age the rates are by. That alone would have saved us the experiment.

## 2. The share pregnant at the start

**What happens now.** `burnin_plan()` gives each woman who could conceive the chance `-np.expm1(-rate*dt*steps)` of being pregnant at the start, with `steps*dt` the mean gestation d rounded to whole steps. In a steady state, the share pregnant in an age group is its birth rate × gestation, h·d. The exponential form is lower: at h = 0.284 per year (India 1950, ages 20–24) and d = 0.769, it gives 0.196 against 0.218, 10% too few.

**What it does:** with change 1 applied, 3,189 women were pregnant at the start in one run of 100,000 agents. The compartmental solver's births per year × gestation give about 3,530. Agents' crude birth rate in 1950 was 0.90 of the compartmental solver's (mean of 10 reps).

**Proposal.** Use the steady-state share, with the unrounded mean gestation:

```python
d = ss.solvers.convert(ss.years(self.mean_dur_pregnancy), engine.sim.unit_base)
return np.minimum(rate*d, 1)*kept
```

Overriding `burnin_plan()` this way (`PregnancyStart` in exp 29's `run.py`) gave 3,489 women pregnant at the start (mean of 10 runs), and 0.5% more population from 2000. The compartmental form's `backfill_burnin()` may want the same check.

## 3. Losses at the mother's death, in the results

**What happens now.** `on_death()` ends the pregnancies of women who die, and counts each as a miscarriage or stillbirth by gestational age (`lose(pregnant, force=False)`, which adds to `engine.flowcounts`). The results don't show them.

**Reproduction:**

```python
sim = ss.Sim(diseases=ss.SIS(beta=0, init_prev=0), demographics=[ss.Pregnancy(fertility_rate=200), ss.death(rate=ss.peryear(0.2))],
             n_agents=20_000, start=2000, stop=2005, dt=ss.months(1), solver='abm', rand_seed=1)
sim.run()
r = sim.results.pregnancy  # miscarriage, stillbirth, abortion, antepartum: all 0
lost = r.P.values[0] + r.conception.values[:-1].sum() - r.delivery.values[:-1].sum() - r.P.values[-1]  # 454; about 488 expected; exactly 0 with no deaths
```

In exp 29's India runs, the four results summed to 0, against 47,131 pregnancies lost this way (6 per 1,000 births in the 1950s, falling to 0.5 by the 2040s).

**Proposal.** Count these losses in the `miscarriage` and `stillbirth` results, or in a separate result (e.g. `loss_maternal_death`), so the four results plus deliveries balance conceptions. We haven't traced why the flow counts added in `on_death()` don't reach the results; it may be that they're added after the step's results are recorded.

These losses don't need offsetting in the fertility rate. A woman who dies while pregnant drops out of the count of women as well as of births, so births per woman alive at delivery already equal the input rate (realized total fertility 1.000). Raising conception to offset them overshoots, to 1.006.

## 4. The warning about double-counted births

**What happens now.** `init_pre()` warns if any declared transition has an `n_born` attribute:

```python
if any(isinstance(mod, Births) for mod in sim.modules) or any(hasattr(tr, 'n_born') for tr in declared):
```

The flow that holds `ss.death` has a transition with `n_born` (set to `None`), so the warning appears in every sim with `ss.Pregnancy` and `ss.death`. Births equal deliveries in those sims, and the population changes by exactly births minus deaths at every step (exp 29, every run).

**Proposal:**

```python
... or any(getattr(tr, 'n_born', None) is not None for tr in declared):
```

## 5. The compartmental solver's starting split by sex

**What happens now.** For a compartmental run, `Sim` takes a named population's ages, but not its sexes: each age group is split by the sex stratification's shares, equal unless given (`init_people()`; "A compartmental run takes the size and age split of a population, not its agents"). Agents take both from the population.

**What it does:** WPP's 1950 population for India is 48.6% female, and 48.1% at ages 15–49. The compartmental solver starts at 50.0% female at every age, so it has 2.5% more women aged 15–49, and 2.5% more births at the start. That's 1.9% more population by 2000 and 2.1% by 2050.

- **The test:** agents started as the compartmental solver is, with WPP's ages and half of each age female, match the compartmental solver (extrapolated to exact aging) within 0.1% through 2050, with changes 1 and 2 applied.
- **Giving overall shares doesn't fix it:** `sim.stratify(sex={'female': 0.486, 'male': 0.514})` corrects the overall split but not the split by age, and it also makes every birth 48.6% female.

**A workaround exists, with a side effect.** Giving the sex split within each age group as a joint split starts the compartmental solver from WPP's population:

```python
sim.stratify(age=edges, sex=ss.levels(['female', 'male'], fractions=ss.by(table, dims=['age', 'sex'])))  # table: one row per age group, female and male shares summing to 1
```

But the compartmental solver's births are then as female as its starting population overall (48.63% here). It ignores `ss.birth`'s `p_female` (only the agent solver uses it), and the youngest age group's row makes no difference either, though a comment in `Groups.init()` says newborns take the youngest group's share. Without a joint split, births follow the sex stratification's shares (50% by default), again whatever `p_female` is. Agents given the same share of girls, through `ss.Pregnancy`'s `p_female`, match the compartmental solver within 0.1% (exp 29). The share matters: 48.6% girls instead of 50% is 2.1% less population by 2000 and 5.5% by 2050.

**Proposal:**

- **The default:** when a sim is stratified by age and sex and given a population with both (`ss.population(country, year)`), start the compartmental solver from the population's joint split by age and sex, as agents do.
- **Births' sex:** take it from `ss.birth`'s `p_female` in the compartmental solver too (or from a sex ratio at birth), independently of the starting split. WPP's sex ratio at birth isn't in `starsim_data` 2.2.0, but a newer version has it.
- **Documentation:** say how the compartmental solver sets the sex of births.

## Not proposed

- **Finer age groups.** With constant-rate aging, single-year groups put India's population 2.2% too high by 2000; quarter-year groups 0.55%, eighth-year 0.3%. That's for users to choose. A note in the documentation of age stratification would help: aging by constant-rate exits spreads a cohort's ages with variance (age × group width).
- **Fertility by age for agents without pregnancy** (`ss.birth(rate=ss.by(age, sex))` for agents, as the compartmental solver has it). That would let both solvers share one birth model when pregnancy itself doesn't matter. With changes 1 and 2, we haven't needed it.
