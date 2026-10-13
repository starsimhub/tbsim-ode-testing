# 05 Exact per-step transitions: summary

## Background

In exp 03, with `beta` calibrated so the ODE equilibrium is 500/100k, the standard ABM (`tbsim.TB`) drifted below the ODE. Exp 04 traced this to tbsim's step ordering. Within a step, tbsim applies transitions one state at a time (INFECTION → NON_INFECTIOUS → ASYMPTOMATIC → SYMPTOMATIC), so agents moving forward can carry on into further states in the same step while agents moving backward stop. The net effect is too little time infectious, which acts like a lower `beta`.

`tbode.TBExact` replaces those ordered transitions. For each agent, it draws the end-of-step state from the exact probabilities of the continuous-time model over one step, `expm(Q·dt)`, where Q holds the ODE's natural-history rates. Transmission and everything else are unchanged.

## Result: `TBExact` matches the ODE at every `dt`

The ABM was started at the ODE equilibrium and run for 100 years, with 50k agents and 40 reps per `dt`. Each rep's outputs are averaged over the last 50 years. The table shows the mean across reps divided by the ODE equilibrium (1.00 = match), ± SE; `tbsim.TB` rows are from exp 03.

| Model, `dt` (reps) | Prevalence, ABM ÷ ODE | Incidence, ABM ÷ ODE | Mortality, ABM ÷ ODE | % asymp., ABM ÷ ODE | Prevalence, 10–90% across reps (ODE = 500) |
|---|---|---|---|---|---|
| `TBExact`, 1 day (40) | **1.00 ± 0.02** | 1.00 ± 0.02 | 1.00 ± 0.02 | 1.00 | 432–567 |
| `TBExact`, 1 week (40) | **1.00 ± 0.02** | 0.99 ± 0.02 | 1.00 ± 0.02 | 1.00 | 423–608 |
| `TBExact`, 1 month (40) | **0.97 ± 0.02** | 0.94 ± 0.02 | 0.97 ± 0.02 | 1.00 | 386–554 |
| `tbsim.TB`, 1 day (40) | 0.95 ± 0.02 | 0.95 ± 0.02 | 0.95 ± 0.02 | 1.00 | 427–542 |
| `tbsim.TB`, 1 week (10) | 0.84 ± 0.03 | 0.84 ± 0.03 | 0.86 ± 0.03 | 1.01 | 374–468 |
| `tbsim.TB`, 1 month (10) | 0.67 ± 0.04 | 0.66 ± 0.04 | 0.71 ± 0.05 | 0.99 | 229–391 |

- **At daily and weekly `dt`, `TBExact` is indistinguishable from the ODE** in prevalence, incidence and mortality (all |z| < 0.5). Its trajectories stay flat around the equilibrium for 100 years (`exact_equilibrium.png`), whereas `tbsim.TB` drifts down throughout.
- **The daily `tbsim.TB` gap is step ordering.** At the same daily `dt`, `TBExact` gives 1.00 against `tbsim.TB`'s 0.95, a difference of 0.05 ± 0.02 (2.0 SE). This answers the open question in exp 03. There is no sign of a gap that doesn't depend on `dt`, such as a finite-population effect near threshold, so no population-size experiment is needed.
- **At monthly `dt`, `TBExact` is slightly low, but prevalence and mortality are within noise.** Prevalence is 0.97 ± 0.02 (z = −1.3) and mortality 0.97 ± 0.02 (z = −1.3). Incidence is 0.94 ± 0.02 (z = −2.8), which is mostly a measurement artefact (confirmed in [exp 11](../11_incidence_counting/SUMMARY.md): the snapshot count misses 3.3% of episodes at monthly `dt`):
  - Incidence is counted as the first entry into ASYMPTOMATIC or SYMPTOMATIC per episode, and the ABM only observes states at the end of each step.
  - So an agent who enters ASYMPTOMATIC and leaves again for NON_INFECTIOUS within the same month isn't seen. Unless it re-enters ASYMPTOMATIC later in the episode, that episode is never counted.
- **INFECTION is low at monthly `dt`, a timing artefact rather than lost time:** 0.89 ± 0.02 of the ODE, against 0.98 at weekly and 1.00 at daily.
  - The force of infection uses the start-of-step state, and newly infected agents then go through that step's transitions. So in the ABM every infection effectively happens at the start of the step, on average half a step earlier than in the ODE, where infections are spread through the step.
  - From then on, each agent spends the correct exponential time in INFECTION. But states are only observed at the end of each step, so each cohort is seen half a step further into its time in INFECTION.
  - At equilibrium that predicts an INFECTION stock of `exp(−r·dt/2)` times the ODE value, with r ≈ 2.1/yr the exit rate from INFECTION. That gives 0.997, 0.98 and 0.92 at daily, weekly and monthly `dt`, matching the observed values.
  - Only INFECTION is affected, because it's the only state whose inflow (infection) is lumped at the start of the step; `expm(Q·dt)` spreads every other inflow correctly through the step.
  - At equilibrium this changes no durations and no other state. Away from equilibrium, it slightly shortens the delay between infection and progression, by half a step.
- **10 vs 40 reps:** the first, 10-rep version of this experiment gave about 0.96 at weekly and monthly `dt`. With 40 reps, that shortfall turned out to be noise.

## Flux check (`run_flux.py`, 20 reps)

This repeats exp 04's per-step check for both models, with the same seeds, at weekly and monthly `dt`: 50k agents, 10 years from the ODE equilibrium, 20 reps.
- For each pair of states (A, B), the ABM's observed frequency of starting a step in A and ending it in B is divided by the exact probability over one step's length, `expm(Q·dt)[A, B]`.
- That exact probability includes every path from A to B within the step. For example, INFECTION → SYMPTOMATIC happens by passing through ASYMPTOMATIC.
- 1.00 means the ABM matches the continuous-time model the ODE represents.
- There's no daily arm. `TBExact` is exact at any `dt` by construction, so this check tests whether its sampling is implemented correctly. Monthly is the hardest test, because multi-step paths within a step are commonest there. At daily `dt` they are too rare to test. In exp 04, INFECTION → SYMPTOMATIC occurred only 5 times in 4 reps.

**Per-flow results** (ABM ÷ exact, ± SE):

| Flow (start of step → end of step) | `tbsim.TB` 1 week | `tbsim.TB` 1 month | `TBExact` 1 week | `TBExact` 1 month |
|---|---|---|---|---|
| **Forward, then on again in the same step** | | | | |
| INFECTION → SYMPTOMATIC (via ASY) | 1.77 ± 0.11 | 1.87 ± 0.06 | 0.94 ± 0.11 | 1.03 ± 0.06 |
| NON_INFECTIOUS → SYMPTOMATIC (via ASY) | 2.01 ± 0.06 | 1.80 ± 0.03 | 0.97 ± 0.06 | 1.04 ± 0.03 |
| ASYMPTOMATIC → DEAD (via SYM) | 1.18 ± 0.05 | 1.45 ± 0.04 | 0.94 ± 0.05 | 1.04 ± 0.04 |
| **Forward into an infectious state** | | | | |
| INFECTION → ASYMPTOMATIC | 0.99 ± 0.01 | 0.95 ± 0.01 | 0.98 ± 0.01 | 1.01 ± 0.01 |
| NON_INFECTIOUS → ASYMPTOMATIC | 0.98 ± 0.01 | 0.90 ± 0.01 | 1.00 ± 0.01 | 1.00 ± 0.01 |
| ASYMPTOMATIC → SYMPTOMATIC | 0.99 ± 0.01 | 0.96 ± 0.01 | 1.00 ± 0.01 | 1.00 ± 0.01 |
| **Backward** | | | | |
| SYMPTOMATIC → ASYMPTOMATIC | 1.03 ± 0.01 | 1.10 ± 0.01 | 1.00 ± 0.01 | 1.02 ± 0.01 |
| ASYMPTOMATIC → CLEARED (via NON) | never (≈69 expected) | never (≈268 expected) | 1.11 ± 0.12 | 1.03 ± 0.06 |
| SYMPTOMATIC → NON_INFECTIOUS (via ASY) | never (≈204 expected) | never (≈792 expected) | 0.90 ± 0.07 | 1.00 ± 0.03 |
| **Other direct transitions** (INFECTION → CLEARED or NON; NON → CLEARED; ASY → NON; SYM → DEAD) | 0.99–1.01 | 1.00–1.02 | 0.98–1.00 | 0.99–1.00 |

- **`tbsim.TB` shows the pattern explained in exp 04.** Agents moving forward carry on in the same step: two-step forward paths are about twice the exact value, and direct entries into ASY and SYM are correspondingly low. Agents moving backward stop: paths that need a backward move followed by another move never occur.
- **`TBExact` matches on every flow,** including the rare multi-step paths that `tbsim.TB` gets wrong or never produces.

**Summary statistic.** Each flow's deviation can be expressed as a z-score: (ratio − 1) ÷ SE, i.e. how many standard errors it is from the exact value. If the only differences are random, each z² averages about 1. The sum covers 14 flows: the table's rows except ASYMPTOMATIC → DEAD, plus the very rare three-step SYMPTOMATIC → CLEARED (via ASY and NON). With only random differences, it should be about 14, and below about 24 in 95% of cases.

| | 1 week | 1 month |
|---|---|---|
| `TBExact`: sum of z² (largest \|z\|) | 19.4 (2.5) | 17.8 (2.3) |
| `tbsim.TB`: sum of z² (largest \|z\|) | 691 (18.2) | 2,688 (28.9) |

`TBExact` is consistent with pure noise; `tbsim.TB` is far from it. For the rare multi-step flows, where the normal approximation is poor, `TBExact` also passes exact Poisson tests (p = 0.06–0.96). In a 4-rep pilot, two of those flows had p ≈ 0.003; that disappeared at 20 reps, so it was noise.
- Transmission is unchanged by `TBExact`: infections of susceptibles are 0.997 and 0.996 (±0.003) of the per-step expectation `1 − e^(−x)`, x = `beta`·w·`dt`. That expectation is itself a fraction of about x/2 below the ODE's infection flow ([exp 10](../10_transmission_step/SUMMARY.md)).

**Unit tests (`tests/test_exact.py`).**
- One monthly step with 200k agents spread over INFECTION, NON, ASY and SYM. End-state frequencies match an independently built `expm(Q·dt)`, including the split of CLEARED agents by pathway (identified by `rr_reinfection`): χ² p = 0.93, 0.63, 0.98, 0.67.
- A control test confirms that `tbsim.TB` fails the same check (p = 0.002 to 10⁻³²).
- `TBExact` raises an error if `k_asy` or `k_non` is > 0.

## Conclusion for the baseline

**The single-strain ODE and ABM are consistent once the ABM's natural history is stepped exactly.** With `TBExact`, they agree within about ±2% at daily and weekly `dt` and within about 3% at monthly, even in this near-threshold regime, where prevalence amplifies any transmission bias about 9×. The whole inconsistency in standard `tbsim.TB` is its step ordering. That biases the ABM towards less transmission:
- about 0.5% in effective `beta` at daily `dt`, 1.8% at weekly and 3.8% at monthly;
- that is, about 5%, 16% and 33% of prevalence here.

## Is monthly `dt` good enough for `TBExact`?

At baseline, yes.
- **Prevalence and mortality** are 0.97 ± 0.02 of the ODE at monthly `dt`, not significantly different from 1.00. The 95% interval is about 0.93–1.01, which allows at most about 0.7% bias in effective `beta`.
- **This calibration is the hardest case:** near threshold, any transmission bias shows up about 9× larger in prevalence.
- **The per-step transitions are exact at monthly `dt`** (flux check).
- **The two visible monthly deviations are measurement artefacts:** INFECTION (0.89) because of infection timing, and incidence (0.94), because brief visits to ASYMPTOMATIC are missed (3.3% of episodes; exp 11).

**What remains at monthly `dt`:**
- **Transmission is still split into steps.** The force of infection is fixed at its start-of-step value, and each infection is dated to the start of its step. Fixing the susceptible pool at the start of the step has a cost. Each agent's per-step infection probability is `1 − e^(−x)` (x = `beta`·w·`dt`), which is a fraction of about x/2 below the ODE's infection flow. So "transmission is exact" holds per step, but not against the ODE ([exp 10](../10_transmission_step/SUMMARY.md)). With this, a deterministic version of the step predicts 0.997 here. Exp 10 reran this arm with 160 reps and got 0.980 ± 0.010, consistent with the 0.97.
- **Incidence counted from monthly snapshots reads low.** `TBExact` now counts it from transitions (`incidence_event_py`, [exp 11](../11_incidence_counting/SUMMARY.md)), which matches the ODE.

**This will matter more with time-since-infection progression.** With `k` ≈ 5.85/yr (tb_LAI_TPT), the progression hazard falls by a factor of e about every 2 months. Dating infections to the start of the step shifts them by about half a month on average, a quarter of that decay time, right where progression is fastest. That experiment should:
- draw each infection's time within the step;
- compare monthly with weekly `dt` directly.

## Options for the user

1. **Fix it upstream in tbsim.** Either sample the end-of-step state from `expm(Q·dt)`, as `TBExact` does, or make a smaller change: draw every agent's transition from its start-of-step state, with at most one transition per step and no carrying on. That still has an O(dt) error, but it no longer favours forward moves.
2. **Keep `tbsim.TB` as is and give the ODE the ABM's step structure,** making the ODE a discrete-time model of the monthly-step ABM.
3. **Use a smaller `dt` in the ABM.** Weekly still leaves a 1.8% transmission bias (16% of prevalence here), at about 4× the run time. Daily leaves 0.5% (5%), at about 30× the run time.

## Caveat for time-since-infection progression

`TBExact` assumes rates that are constant within a step and the same for all agents. With `k_asy`/`k_non` > 0, exits from INFECTION depend on each agent's time since infection. Exact sampling would then need a transition matrix per agent, or agents grouped by time since infection; both are doable. Separately, ABM cohort progression observed once per step, as in tb_LAI_TPT's `cohort_progression`, misses entries to ASY that leave again within the step, so it reads slightly low at monthly `dt`.
