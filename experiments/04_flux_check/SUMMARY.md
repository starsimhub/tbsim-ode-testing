# 04 Flux check: summary

## Background and method

In exp 03, with `beta` calibrated so the ODE equilibrium is 500/100k, the ABM (`tbsim.TB`) drifted below the ODE: prevalence was 5% low at daily `dt`, 16% at weekly and 33% at monthly. This experiment asks which part of the ABM causes that, transmission or the natural-history transitions, by checking each **one-step** transition against what the ODE implies.

The ABM was started at the ODE equilibrium and run for 10 years, with 50k agents and 4 reps at each `dt`. An analyzer recorded every living agent's state at the start and end of each step. Two comparisons follow.

1. **Transmission.** Realised infections among SUSCEPTIBLE and CLEARED agents are compared with the number the ODE force of infection predicts, computed from the ABM's own infectious stocks at the start of each step.
2. **Natural history.** The ODE is a continuous-time model, so it implies an exact probability that an agent starting a step in state A ends it in state B. These probabilities are the entries of the matrix exponential `expm(Q·dt)`, where Q holds the ODE's transition rates. They include **every path** from A to B within the step, not only direct transitions. For example, the ODE has no direct INFECTION → SYMPTOMATIC rate. The ODE has no steps or ordering; all flows run continuously and at once. So someone who reaches ASYMPTOMATIC early in the period can go on to SYMPTOMATIC before it ends, and the probability of INFECTION → SYMPTOMATIC over one step's length is small but non-zero. It is about `inf_asy · asy_sym · dt² / 2`: second order in `dt`, because it needs two moves, and halved because the second move only has the time left after the first. For each pair (A, B), the table reports the ABM's observed one-step frequency divided by this exact probability, **so 1.00 means a match**.

## How tbsim steps through the natural history

Within each step, `tbsim.TB` applies transitions one state at a time, in the order INFECTION → NON_INFECTIOUS → ASYMPTOMATIC → SYMPTOMATIC, and each change takes effect immediately. This has two consequences:
- **Moving forward, an agent can continue in the same step.** An agent who moves INFECTION → ASYMPTOMATIC is then processed again as ASYMPTOMATIC, with the full step's probability of moving on, as if it had been ASYMPTOMATIC since the start of the step. In continuous time it would arrive part-way through the step, so on average it would face only about half of that probability.
- **Moving backward, an agent stops for the step.** An agent who moves ASYMPTOMATIC → NON_INFECTIOUS has already passed the NON_INFECTIOUS stage of the step, so it can't move again until the next step.

The results below show both effects.

## Results

**Transmission is matched at every `dt`.** Realised infections among SUSCEPTIBLE agents divided by the ODE prediction are 1.005, 1.000 and 1.005 at daily, weekly and monthly `dt` (SE ±0.007). So the conversion from ODE `beta` to the ABM's per-contact `beta` (`beta_edge_from_ode`) is correct.

The same check among CLEARED agents gives 0.99, 0.96 and 0.85 (in `transmission.csv`). That apparent shortfall is a counting artefact, not a transmission error:
- An infection is counted when an agent's state changes during the step.
- Transmission runs before the transitions, so a reinfected CLEARED agent can immediately go INFECTION → CLEARED in the same step. It then ends the step where it started and isn't counted.
- The shortfall matches the one-step probability of INFECTION → CLEARED: 0.036 at weekly `dt` and 0.15 at monthly.
- SUSCEPTIBLE agents aren't affected, because a SUSCEPTIBLE agent infected and cleared within the step ends in CLEARED, which still counts as a change.

**Natural-history transitions are biased, in the pattern the step ordering predicts.** Values are ABM ÷ exact one-step probability (1.00 = match), ± SE, from 4 reps. Each column compares the ABM run at that `dt` with the ODE's probability over an interval of the same length. For the two-step flows the ratio stays near 2 as `dt` shrinks, because both probabilities scale with `dt²`. What shrinks is the size of those flows.

| Flow (start of step → end of step) | 1 day | 1 week | 1 month | Why |
|---|---|---|---|---|
| **Forward, then on again in the same step** | | | | |
| INFECTION → SYMPTOMATIC (via ASY) | 1.9 ± 0.6 | 1.8 ± 0.2 | 1.9 ± 0.1 | Entrants to ASY face the full step's onward probability, about twice the continuous-time value |
| NON_INFECTIOUS → SYMPTOMATIC (via ASY) | 2.2 ± 0.3 | 1.9 ± 0.1 | 1.8 ± 0.1 | Same |
| ASYMPTOMATIC → DEAD (via SYM) | 1.1 ± 0.1 | 1.4 ± 0.1 | 1.5 ± 0.1 | Same, for entrants to SYM |
| **Forward into ASY** (some move on to SYM in the same step) | | | | |
| INFECTION → ASYMPTOMATIC | 1.04 ± 0.02 | 1.00 ± 0.02 | 0.95 ± 0.03 | Lost to the onward paths above |
| NON_INFECTIOUS → ASYMPTOMATIC | 1.00 ± 0.01 | 0.99 ± 0.01 | 0.90 ± 0.01 | Same |
| ASYMPTOMATIC → SYMPTOMATIC | 0.97 ± 0.02 | 0.99 ± 0.02 | 0.96 ± 0.02 | Some go SYM → ASY in the same step and end where they started |
| **Backward, so no further move** | | | | |
| SYMPTOMATIC → ASYMPTOMATIC | 1.01 ± 0.02 | 1.02 ± 0.02 | 1.06 ± 0.02 | Entrants to ASY can't move on in the same step, so more end there |
| ASYMPTOMATIC → CLEARED (via NON): observed / expected count | 0 / 2 | 0 / 14 | 0 / 53 | Needs ASY → NON then NON → CLEARED, which can't happen in one step |
| SYMPTOMATIC → NON_INFECTIOUS (via ASY): observed / expected count | 0 / 6 | 0 / 41 | 0 / 159 | Needs SYM → ASY then ASY → NON |
| **Other direct transitions** (to CLEARED, to DEAD, INFECTION → NON_INFECTIOUS) | 0.98–1.05 | 0.96–1.05 | 1.00–1.10 | Little affected. All within about 2 SE of 1, except NON_INFECTIOUS → CLEARED at monthly `dt` (1.03 ± 0.01) |

At monthly `dt` the deviations are clear. At daily `dt` most of them are within noise.

## Interpretation

The net effect of the step ordering is that **agents spend too little time infectious**:
- Agents who move forward into ASYMPTOMATIC or SYMPTOMATIC lose about half a step there.
- Agents who move backward out of them can't move again in that step, so they stay in the earlier, non-infectious states longer.

Less time infectious means less transmission. That's equivalent to a 2–4% lower `beta` at weekly to monthly `dt`, and about 0.5% at daily (exp 03). Because `beta` = 9.917 is close to the epidemic threshold (`beta` ≈ 8.59), prevalence is about 9× as sensitive as `beta`. A 2–4% loss of transmission therefore becomes the 16–33% prevalence shortfall seen in exp 03.

Newly infected agents also go through transitions in the step they are infected, because transmission runs before the transitions. This doesn't shorten their time in INFECTION. It dates each infection to the start of the step, on average half a step earlier than if infections were spread through the step as in the ODE, so end-of-step snapshots of INFECTION read low by about `exp(−r·dt/2)`, with r ≈ 2.1/yr the exit rate from INFECTION (exp 05). It doesn't change where agents go next, and INFECTION isn't infectious.

**Resolution at daily `dt`.** This check can't resolve the gap at daily `dt`. In exp 03, with 40 reps, daily `tbsim.TB` is 5% below the ODE in prevalence, about 0.5% in effective `beta`. Here, each daily flow ratio has an SE of 1–2%, so a bias that small is below their resolution. Exp 05 settles it another way: at daily `dt`, `TBExact`, which removes the step ordering, matches the ODE (1.00 ± 0.02).

## Next

Replace the ordered transitions with sampling each agent's end-of-step state from the exact one-step probabilities `expm(Q·dt)`: [05_exact_transitions](../05_exact_transitions/).
