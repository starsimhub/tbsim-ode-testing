# 09 Population size: summary

## Background

`tbode.TBExact` is a version of `tbsim.TB` that samples each agent's end-of-step state from the exact one-step probabilities `expm(Q·dt)`, so it has no step-ordering bias. It uses Poisson-mixture latent classes for progression by time since infection and dates each new infection within its step. In exp 08, set A (tbsim defaults + `k` = 5.846, near threshold) at monthly `dt` gave 0.937 ± 0.018 of the ODE prevalence with 50k agents. A deterministic version of the same step (`step_map.py`) predicts 0.995, with no sampling noise or finite-population effects. Exp 07 had a similar, smaller shortfall (weekly, `k` = 0). The hypothesis tested here was a finite-population effect: near threshold, random fluctuations might hold the stochastic mean below the deterministic equilibrium. Such an effect shrinks roughly like 1/N.

## Result: the shortfall doesn't shrink with N, so it's not a finite-population effect

Set A, monthly `dt`, `TBExact`, 100 years from the ODE equilibrium. Each rep's prevalence is averaged over the last 50 years. Values are ABM ÷ ODE (1.00 = match), mean ± SE across reps.

| N | Reps | Prevalence ÷ ODE | Shortfall vs the deterministic step (0.995) | Rep-to-rep SD | Fade-outs |
|---|---|---|---|---|---|
| 25k | 190 | 0.947 ± 0.011 | 4.8% | 14.5% | 0 |
| 50k | 94 | 0.924 ± 0.010 | 7.1% | 10.1% | 0 |
| 100k | 48 | 0.964 ± 0.011 | 3.1% | 7.6% | 0 |
| 200k | 24 | 0.941 ± 0.007 | 5.4% | 3.3% | 0 |

The 50k row includes exp 08's 40 reps (seeds 0–39). The current code reproduced seed 0 exactly.

- **Fits of the shortfall** (weighted by 1/SE²; χ² on 3 degrees of freedom):
  - constant: χ² = 7.0, about 5.2%;
  - ∝ 1/√N: χ² = 26;
  - ∝ 1/N: χ² = 60.
  
  Only the constant fits acceptably (p ≈ 0.07). The finite-population prediction (12% / 6% / 3% / 1.5%) is ruled out. At 200k agents the shortfall is still 5.4% ± 0.7%.
- **Incidence, mortality and INFECTION** show the same pattern, about 5–8% low at every N. (Monthly incidence also carries its known snapshot artefact.)
- **The time course is the same at every N** (`population_size.png`): the median prevalence declines slowly and smoothly from 500 to about 470 over the century, so the ABM's long-run level is below the ODE's.
- **The rep-to-rep SD shrinks with N** (14.5% at 25k to 3.3% at 200k), roughly as 1/√N or faster. So random fluctuations do get smaller, but the mean shortfall doesn't.

## Interpretation (resolved in exp 10)

`TBExact` in set A at monthly `dt` has a real, systematic shortfall of about 5% in prevalence, and it isn't a finite-population effect. [Exp 10](../10_transmission_step/SUMMARY.md) found the cause: the deterministic step used here assumed a linear per-step infection probability, x = `beta`·w·`dt`. The ABM's network transmission actually gives `1 − e^(−x)`, because the susceptible pool is fixed at the start of the step and each agent can be infected at most once per step. That transmits about x/2 less than the ODE does. With the correct form, the deterministic step predicts 0.941 for this arm, matching the observed 0.943 ± 0.005 (all N, inverse-variance weighted). With a linear infection probability (mass action), the shortfall largely disappears (0.984 ± 0.010; exp 10). It doesn't depend on N because it's a time-step effect.

## Conclusion

Set B, the reference set, still matches (±1%). The near-threshold set A shortfall is not a finite-population effect. It's the ABM's per-step transmission, `1 − e^(−x)` with a susceptible pool fixed at the start of the step, which falls short of the ODE's flow by a fraction of about x/2 (exp 10).

## Versions

tbsim main @ `b92d7e8` (v0.10.3), starsim @ `3521adb4` (v3.6.1). `TBExact` unchanged since exp 08 (seed 0 reproduced exactly).
