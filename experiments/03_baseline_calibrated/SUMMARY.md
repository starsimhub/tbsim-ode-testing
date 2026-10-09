# 03 Baseline calibrated to 500/100k: summary

## Calibration

- ODE `beta` = **9.917** gives 500/100k at equilibrium. 200/100k needs 9.193 and 800/100k needs 10.494.
- **This regime is just above the epidemic threshold.** The threshold is `beta` ≈ 8.59: below it, TB dies out. That has three consequences:
  - **Prevalence is very sensitive to `beta`.** Near threshold, prevalence grows roughly in proportion to how far `beta` is above 8.59, so a 1% change in `beta` changes prevalence by about 9%. For example, 1% less `beta` gives 454/100k instead of 500, and the whole 200–800/100k range spans only 14% in `beta`.
  - **Small transmission biases show up as large prevalence gaps.** That makes this a stringent test of consistency, but a gap that would be negligible far from threshold looks large here.
  - **Approach to equilibrium is slow.** Transmission barely exceeds what's needed to replace cases, so the pull back towards equilibrium is weak. From a 1% SYMPTOMATIC seed, the ODE takes about 250–300 years to settle. The ABMs are therefore started at the ODE equilibrium rather than burned in.
- tb_LAI_TPT's calibrated rates were not used: without `k`, they put 42% of a latent cohort into active disease within 10 years.

## Are the implied outputs reasonable? (ODE equilibrium, tbsim default natural history, no treatment)

| Output | Value | Comment |
|---|---|---|
| Prevalence (ASY+SYM), all ages | 500/100k | Target |
| Incidence (first ASY entry per episode) | 273/100k/yr | Prevalence/incidence = 1.83 years, a plausible untreated disease duration |
| TB mortality | 84/100k/yr | Mortality/incidence = 0.31. Plausible for no treatment; WHO India 2000, with treatment, was 99/100k |
| % asymptomatic | 50% | Below the India survey's ~65% (tb_LAI_TPT target). Set by ASY↔SYM rates, not `beta` |
| ProgressionY1 | 0.034 | Target 0.052: too little early progression without time-since-infection progression |
| ProgressionRatio | 0.48 | Target 0.66: same reason |
| Progression by 10 years | 7.1% | Between Ferebee (6.2%) and Sutherland (9.6%) |

## ABM vs ODE

The ABM (`tbsim.TB`, 50k agents) was started at the ODE equilibrium and run for 100 years. Each rep's outputs are averaged over the last 50 years. The table shows the mean across reps divided by the ODE equilibrium (1.00 = match), ± SE. Daily has 40 reps; weekly and monthly have 10 each.

| ABM `dt` (reps) | Prevalence, ABM ÷ ODE | Incidence, ABM ÷ ODE | Mortality, ABM ÷ ODE | % asymp., ABM ÷ ODE | Prevalence, 10–90% across reps (ODE = 500) | Effective-`beta` shortfall |
|---|---|---|---|---|---|---|
| 1 day (40) | 0.95 ± 0.02 | 0.95 ± 0.02 | 0.96 ± 0.02 | 1.00 | 427–542 | 0.52 ± 0.18% |
| 1 week (10) | 0.84 ± 0.03 | 0.84 ± 0.03 | 0.86 ± 0.03 | 1.00 | 374–468 | 1.75 ± 0.35% |
| 1 month (10) | 0.67 ± 0.04 | 0.66 ± 0.04 | 0.71 ± 0.05 | 0.99 | 229–391 | 3.80 ± 0.52% |

The effective-`beta` shortfall is how much lower the ODE `beta` would have to be to give the ABM's mean prevalence. It expresses each gap as a transmission bias, which removes the threshold amplification described above.

- **At every `dt`, the ABM drifts below the equilibrium and is still falling after 100 years.** The drift is smooth. With 40 daily reps, prevalence averaged over successive 25-year windows is 498, 493, 485 and 467/100k (SE 5–9). Earlier, with only 10 daily reps, the median seemed to hold at 500 until about 1960 and then fall; with 40 reps that turned out to be noise.
- **The gap shrinks as `dt` falls, but daily `dt` still doesn't fully match.** Daily is 2.9 SE below the ODE, so the gap isn't noise, although the ODE lies inside the daily 10–90% range. The flux check in [04_flux_check](../04_flux_check/), which compares per-step transition rates in the two models, identifies step ordering in `tbsim.TB` as the cause, and exp 05 confirms that it accounts for all of the gap (see below). Within a step, tbsim applies transitions in a fixed order (INFECTION → NON_INFECTIOUS → ASYMPTOMATIC → SYMPTOMATIC). An agent moving forward can carry on through several states in one step, while one moving backward stops for the step. Agents therefore spend too little time infectious, which acts like a slightly lower `beta`.
- **Cohort progression** is 1.00–1.02× the exact value at daily `dt`, 0.96–0.97× at weekly and 0.91–0.95× at monthly. Part of the shortfall comes from observing states only once per step.
- **% asymptomatic matches at every `dt`.** It depends only on the ASY↔SYM rates, not on transmission.

## Is the whole daily gap due to `dt`?

**Yes: it is step ordering.** The answer comes from [05_exact_transitions](../05_exact_transitions/). `tbode.TBExact` is a version of `tbsim.TB` that removes the step-ordering bias. It draws each agent's end-of-step state from the exact continuous-time transition probabilities over that step (`expm(Q·dt)`), rather than applying transitions one after another. With 40 reps, the same setup gives:
- `TBExact` prevalence: 1.00 ± 0.02 of the ODE at daily `dt`, 1.00 ± 0.02 at weekly and 0.97 ± 0.02 at monthly;
- against `tbsim.TB`'s 0.95 ± 0.02 at daily `dt`.

So there's no part of the gap that doesn't depend on `dt`, such as a finite-population effect near threshold.

The shortfall in `beta` grows less than proportionally with `dt`: 0.52 ± 0.18% (daily), 1.75 ± 0.35% (weekly), 3.80 ± 0.52% (monthly). It doesn't fit a simple proportional-to-`dt` form. One reason may be that the ABM is still drifting at year 100, so the last-50-year mean understates each run's eventual equilibrium, by different amounts at each `dt`.

## How the ABM starts at equilibrium

All ABM runs in exps 03–05, both `tbsim.TB` and `TBExact`, start from the same ODE equilibrium; the ABM never searches for its own.
- **The equilibrium** (`TBODE.equilibrium`) comes from integrating the ODE's population proportions for 2,000 years from a 1% SYMPTOMATIC seed, then polishing with a root solve.
- **It's assigned to agents by `init_from_ode`,** after `sim.init()` and before the run:
  - Each agent's state is drawn at random with the ODE's proportions as probabilities. The ODE's CLEARED, RECOVERED and TREATED compartments map to the ABM's CLEARED state, with the matching reinfection risk (`rr_reinfection`).
  - Infected agents get `ever_infected` set and a time of infection (`ti_infected`) of 0. For the incidence counter, the right share of NON_INFECTIOUS agents is marked as already active in this episode.
  - Each rep therefore starts with slightly different proportions (e.g. about 126 ± 11 ASYMPTOMATIC agents), which is part of the rep-to-rep spread.

**Why no burn-in is needed:** at baseline, an agent's future depends only on its current state. Rates are constant, and births and deaths don't depend on age. So the state proportions fully describe the equilibrium, and `TBExact` stays at it.

**When this stops being enough:** any feature that gives agents memory. With time-since-infection progression (`k` > 0), `ti_infected = 0` would make every latent agent look newly infected and cause a burst of early progression. Initial times since infection will need drawing from the ODE's time-since-infection sub-compartments. Age-specific rates and fixed-duration treatment and relapse delays will likewise need an initial age distribution and time in state.

## Notes

- The daily arm and the 30 extra daily reps were added after the first run, which had weekly and monthly only. The weekly and monthly arms were rerun with the same seeds and reproduced the earlier numbers exactly.
- **Figures:** `baseline_equilibrium.png` (10 reps per `dt`) and `baseline_equilibrium_more_reps.png` (40 daily reps). Their data are in `outputs.obj` and `outputs_more_reps.obj`. Tables: `summary.csv` (medians, 10 reps) and `summary_more_reps.csv` (means ± SE, used above).
- The daily incidence and mortality curves drop to 0 at 2000, because the year-2000 bucket contains only the final step, which records no events. This has a negligible effect on the 50-year means.

## Next

Find which flows cause the shortfall: [04_flux_check](../04_flux_check/).
