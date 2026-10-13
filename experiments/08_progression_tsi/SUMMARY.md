# 08 Progression by time since infection: summary

## Background

tbsim's INFECTION-exit hazards can decline with time since infection τ: `inf_non·e^(−kτ)` and `inf_asy·e^(−kτ)`. tb_LAI_TPT's best fit has `k` = 5.846/yr, so progression is strongly front-loaded. The ODE represents this exactly with Poisson-mixture latent classes:
- each new infection joins class m with probability Poisson(m; r₀/k), where r₀ = `inf_non` + `inf_asy`;
- class m progresses at the constant rate m·k.

7–8 classes are needed (see the [README](README.md)). Two ABMs are compared with the ODE:
- **`TBExact`:** samples each agent's end-of-step state from the exact one-step probabilities `expm(Q·dt)`, so it has no step-ordering bias. It uses the same latent classes and dates each new infection within its step (exp 07).
- **Standard `tbsim.TB`:** applies transitions one state at a time and tracks each agent's actual τ.

There are two parameter sets, each with `beta` calibrated to 500/100k:
- **A:** tbsim defaults + `k`. It's near threshold: prevalence moves 8.6× as much as `beta`.
- **B:** tb_LAI_TPT's natural history + `k`, the reference set. Prevalence moves 2.5× as much as `beta`.

Both ABMs use 50k agents, start at the ODE equilibrium with times since infection drawn for latent agents, and run for 100 years; outputs are averaged over the last 50 years.

## Results

### `TBExact` vs the ODE: the ODE matches the feature

Values are ABM ÷ ODE equilibrium (1.00 = match), mean ± SE over 40 reps.

| Output | A, weekly | A, monthly | B, weekly | B, monthly |
|---|---|---|---|---|
| Prevalence (ASY + SYM) | 0.99 ± 0.02 | **0.94 ± 0.02** | **1.00 ± 0.01** | **1.01 ± 0.01** |
| INFECTION | 0.99 ± 0.02 | 0.94 ± 0.02 | 1.00 ± 0.01 | 1.00 ± 0.01 |
| Incidence | 0.98 ± 0.02 | 0.91 ± 0.02 | 1.00 ± 0.01 | 0.97 ± 0.01 |
| TB mortality | 0.99 ± 0.02 | 0.94 ± 0.02 | 1.00 ± 0.01 | 1.01 ± 0.01 |
| % asymptomatic | 1.00 | 1.00 | 1.00 | 1.00 |
| Prevalence, 10–90% across reps (ODE = 500) | 432–552 | 407–531 | 463–536 | 471–536 |
| Deterministic `TBExact` step ÷ ODE, ABM's `1 − e^(−x)` infection probability (exp 10) | 0.986 | 0.941 | 0.999 | 0.994 |

The last row is a deterministic version of the `TBExact` step (`step_map.py` here, corrected in exp 10): the expected start-of-step infections, then exact transitions with latent classes and `'within'` timing, iterated to equilibrium. Here it uses the ABM's per-step infection probability, `1 − e^(−x)` with x = `beta`·w·`dt`. This experiment's `step_map.py` originally used the linear x and predicted 0.999 / 0.995 / 1.000 / 0.999.

- **For set B, the reference set, the ODE matches `TBExact` to within ±1% at both `dt`s.** This includes monthly `dt`, tb_LAI_TPT's setting. Monthly incidence (0.97) is the snapshot artefact seen in exps 05 and 07: entries to ASYMPTOMATIC that leave again within the step are missed. Counted from transitions, it matches ([exp 11](../11_incidence_counting/SUMMARY.md)).
- **For set A at weekly `dt`, the ODE also matches** (0.99 ± 0.02).
- **For set A at monthly `dt`, `TBExact` is 6% low (z = −3.5), and the cause is transmission, not progression.**
  - The ABM's per-step infection probability is `1 − e^(−x)`, with the susceptible pool fixed at the start of the step. That transmits a fraction of about x/2 less than the ODE: 0.7% of effective `beta` here, which near threshold is 6% of prevalence. The deterministic step with that form predicts 0.941 ([exp 10](../10_transmission_step/SUMMARY.md)).
  - It isn't a finite-population effect ([exp 09](../09_population_size/SUMMARY.md)).
  - With a linear per-step probability (mass action), set A gives 0.984 ± 0.010 and set B 1.000 ± 0.006 (exp 10).
- **Starting from the equilibrium works.** Drawing each latent agent's time since infection (and class) from the ODE equilibrium gives flat `TBExact` trajectories from the first year, with no early burst of progression.

**Cohort progression** (`run_cohort.py`; 4 × 100k agents; ABM ÷ exact, where the exact values come from integrating the time-varying hazard directly, without classes):

| Model, `dt` | A: ProgressionY1 | A: 10-year progression | B: ProgressionY1 | B: 10-year progression | B: ProgressionRatio |
|---|---|---|---|---|---|
| `TBExact`, weekly | 0.97 ± 0.01 | 0.98 ± 0.01 | 0.98 ± 0.01 | 0.99 ± 0.01 | 1.00 |
| `TBExact`, monthly | 0.94 ± 0.01 | 0.97 ± 0.01 | 0.94 ± 0.01 | 0.97 ± 0.01 | 0.97 |
| `tbsim.TB`, daily | **1.00 ± 0.01** | **1.02 ± 0.01** | **0.99 ± 0.01** | **0.99 ± 0.01** | **1.00** |
| `tbsim.TB`, weekly | 1.03 ± 0.01 | 1.05 ± 0.01 | 1.03 ± 0.01 | 1.04 ± 0.01 | 0.99 |
| `tbsim.TB`, monthly | **1.15 ± 0.02** | **1.19 ± 0.01** | **1.10 ± 0.01** | **1.16 ± 0.01** | 0.95 |

- **Daily `tbsim.TB` matches the exact curve.** It tracks each agent's actual τ, so this confirms the latent-class representation independently.
- **`TBExact`'s shortfall is the snapshot artefact.** Cohorts are observed only at the end of each step, so an agent who enters ASYMPTOMATIC and returns to NON_INFECTIOUS within the step is missed. The expected loss is about `asy_non·dt/2`: 7% monthly and 1.6% weekly, matching the 0.94 / 0.98.

### `tbsim.TB`: progression by time since infection makes step effects much larger

| Prevalence (40 reps) | A, weekly | A, monthly | B, weekly | B, monthly |
|---|---|---|---|---|
| `tbsim.TB` ÷ ODE | 1.30 ± 0.02 | **2.64 ± 0.02** | 1.18 ± 0.01 | **1.83 ± 0.01** |
| `tbsim.TB` ÷ `TBExact` | 1.32 ± 0.03 | 2.81 ± 0.06 | 1.17 ± 0.02 | 1.82 ± 0.02 |
| Baseline, `k` = 0, tbsim defaults (exps 03, 05) | 0.84 | 0.67 | — | — |

Incidence, mortality and NON_INFECTIOUS move with prevalence. INFECTION moves less (`tbsim.TB` ÷ `TBExact` = 1.46 for B monthly). % asymptomatic is 1% lower.

- **The ABM gap reverses and grows.** At baseline, step ordering made `tbsim.TB` 16% / 33% *low* at weekly / monthly `dt`. With `k` = 5.846 it's 18–30% *high* at weekly `dt`, and 1.8–2.6× at monthly.
- **The cause is how `tbsim.TB` evaluates the declining hazard.** It computes τ at the start of each step and applies that hazard for the whole step. Over a step, the hazard actually falls by a factor of e^(−k·dt): 0.61 at monthly `dt`. So progression in the first months after infection, where most of it happens, is overstated: the integrated hazard is 1.26× the exact value at monthly `dt`, and the cohort runs show 10–19% more progression.
  - This outweighs the step-ordering bias, which pushes the other way.
  - Its effect on prevalence is large because extra progression acts like extra transmission. Even in set B, further from threshold, monthly prevalence is 1.8× the ODE.
- **This is the setting tb_LAI_TPT runs in: monthly `dt` with `k` ≈ 5.85.** tb_LAI_TPT's calibrated parameters therefore produce their fitted epidemic only *together with* this discretization bias. The same parameters in an ODE, or in `TBExact`, give a much smaller epidemic: about 0.55× prevalence for set B's natural history at the same `beta`. Its calibrated rates and `beta` are partly compensating for the time step.

### Progression targets in the exact model: no recalibration of `k` needed

tb_LAI_TPT calibrated `k` (with the other natural-history rates) to two cohort targets taken from data (Ferebee 1970, Sutherland 1968): ProgressionY1 = 0.0523 and ProgressionRatio = 0.663, with log-scale SDs of 0.094 and 0.071. It measured them with a monthly `tbsim.TB` cohort (`tb_LAI_TPT/cohort_progression.py`), which overstates progression as shown above. Measured exactly, set B is closer to the targets:

| Set B's parameters, measured with… | ProgressionY1 | ProgressionRatio | Distance from targets, in target SDs |
|---|---|---|---|
| Exact cohort (ODE, or `TBExact` without the snapshot artefact) | 0.0511 | 0.674 | −0.25, +0.24 |
| Monthly `tbsim.TB` cohort (tb_LAI_TPT's measurement) | 0.0563 | 0.639 | +0.78, −0.52 |

- **No drift away from the targets to correct.** The targets come from data, not from tbsim's dynamics. Without the monthly bias, the same parameters fit the data better.
- **Forcing an exact fit would move the parameters a lot for little gain.** A preview fitted `k` and a common scale on `inf_non`/`inf_asy` so the exact cohort hits both targets. It gave `k` = 4.16 (from 5.85) and rates × 0.77, requiring `beta` = 23.6 (from 25.0) for 500/100k. That's a large move to close gaps of a quarter of an SD each; the two targets constrain `k` only weakly.
- **ODE↔`TBExact` consistency doesn't depend on `k`'s value:** the latent classes are exact for any `k`.

**Decision (user, 2026-10-12):** keep set B as it is. tb_LAI_TPT as a whole will need recalibrating with the corrected ABM if it is to be consistent with the ODE; that's out of scope for now.

## Conclusions

1. **The ODE can match progression by time since infection exactly.** The Poisson-mixture latent classes are exact for `k_asy = k_non`, need only 7–8 compartments, and match `TBExact` to within ±1% for the reference set B at both `dt`s.
2. **Standard `tbsim.TB` at monthly `dt` can't be matched by an ODE with the same parameters.** Evaluating the hazard at the start of each step inflates early progression, giving 1.8–2.6× the ODE prevalence. Weekly `dt` still gives 1.2–1.3×. The bias depends on `dt` and grows with `k·dt`, so an ODE could only imitate it approximately, with a `dt`-specific adjustment to the rates. For tb_LAI_TPT, this is now the largest ABM↔ODE inconsistency found, much larger than the baseline step-ordering bias.
3. **How to handle the bias (user decision, 2026-10-12): no changes to tbsim.** `TBExact` is the corrected ABM, and it accumulates the fixes as features are added. tb_LAI_TPT will need recalibrating with it to be consistent with the ODE, because its parameters currently compensate for the bias; that's deferred. Two upstream fixes were considered and not pursued:
   - the hazard averaged over each step, `e^(−kτ)·(1 − e^(−k·dt))/(k·dt)`, which would leave the smaller step-ordering bias;
   - exact sampling in `tbsim.TB` itself.
4. **No recalibration of `k`.** Measured exactly, set B fits the data-derived progression targets better than tb_LAI_TPT's monthly `tbsim.TB` measurement does (see above).
5. **Set A's monthly shortfall (6%) is a transmission time-step effect** ([exp 10](../10_transmission_step/SUMMARY.md)): the ABM's `1 − e^(−x)` per-step infection probability, with a fixed start-of-step susceptible pool. It isn't related to progression by time since infection.

## Versions

tbsim main @ `b92d7e8` (v0.10.3), starsim @ `3521adb4` (v3.6.1). `../tbsim` had an uncommitted change, but only to `tbsim/compartmental/README.md` (documentation, not code).
