# 10 Transmission within a step: summary

## Background

`tbode.TBExact` is a version of `tbsim.TB` that samples each agent's end-of-step state from the exact one-step probabilities `expm(Q·dt)`. It has no step-ordering bias, uses Poisson-mixture latent classes for progression by time since infection, and dates each new infection within its step. Even so, it ran a few percent below the ODE in some near-threshold arms: about 5% in exp 08's set A (tbsim defaults + `k` = 5.846) at monthly `dt`. Exp 09 showed the shortfall doesn't depend on population size.

The comparisons in exps 07–09 used a deterministic version of the `TBExact` step (`step_map.py`), which assumed each susceptible's infection probability per step is linear: x = `beta`·w·`rel_sus`·`dt`, where w is the `trans_asymp`-weighted infectious fraction. The ABM's actual per-step probability is `1 − e^(−x)`: on the random network, an agent with j infectious contacts is infected with probability `1 − (1 − q)^j`, which averages to that.

The two forms differ because of how the susceptible pool is treated within a step:
- **ODE:** at equilibrium, infections over a step are `beta·w·S·dt`, because people who become susceptible during the step, for example by clearing a latent infection, can be infected later in the same step.
- **ABM:** the susceptible pool is fixed at the start of the step, and each agent can be infected at most once per step. So infections are `S·(1 − e^(−x))`, about a fraction x/2 fewer.

This experiment tested whether that explains the shortfall, using a mass-action version of `TBExact` (`tbode.TBExactMassAction`) with either per-step probability.

## Results

Values are prevalence, ABM ÷ ODE (1.00 = match), averaged over each rep's last 50 years, as mean ± SE across reps. The prediction comes from the deterministic step with the arm's infection-probability form. All arms are monthly `dt` with `'within'` infection timing unless noted.

| Arm | N | Reps | ABM ÷ ODE | Prediction: linear | Prediction: `1 − e^(−x)` | z vs the arm's own form |
|---|---|---|---|---|---|---|
| 1. Set A, network | 25k | 190 | 0.947 ± 0.011 | 0.995 | **0.941** | +0.6 |
| 2. Set A, mass action, `1 − e^(−x)` | 25k | 190 | 0.941 ± 0.010 | 0.995 | **0.941** | 0.0 |
| 3. Set A, mass action, linear | 25k | 190 | **0.984 ± 0.010** | **0.995** | 0.941 | −1.1 |
| 4. Set B, network | 25k | 190 | 0.989 ± 0.005 | 0.999 | **0.994** | −1.0 |
| 5. Set B, mass action, linear | 25k | 190 | **1.000 ± 0.006** | **0.999** | 0.994 | +0.3 |
| 6. `k` = 0, `'start'`, network (exp 05's arm) | 50k | 160 | 0.980 ± 0.010 | 1.014 | **0.997** | −1.7 |
| 7. `k` = 0, **weekly**, mass action, `1 − e^(−x)` | 50k | 80 | 0.952 ± 0.017 | 0.999 | **0.995** | −2.5 |

- **Set A is tbsim defaults + `k` = 5.846 (near threshold); set B is tb_LAI_TPT's natural history + `k` (the reference set).** Both are from exp 08. The `k` = 0 arms use tbsim defaults at baseline `beta`.
- **Arm 1 reuses exp 09's 25k runs; arm 6 includes exp 05's 40 reps.** Seed 0 of each was reproduced exactly.
- **No rep faded out.**
- **Incidence, mortality and INFECTION** move with prevalence (`summary.csv`).

### Findings

1. **The near-threshold shortfall is caused by the per-step infection probability.**
   - In set A, mass action with `1 − e^(−x)` gives 0.941 ± 0.010, exactly the deterministic prediction (0.941).
   - Switching to the linear form raises prevalence to 0.984 ± 0.010, consistent with its prediction of 0.995 (z = −1.1).
   - So the ABM's saturating infection probability, together with a susceptible pool fixed at the start of the step, transmits about x/2 too little compared with the ODE. In set A at monthly `dt`, x ≈ 0.014, so that's 0.7% of effective `beta`. Near threshold, prevalence moves 8.6× as much as `beta`, which turns this into 6% of prevalence.
2. **The random network adds nothing.** Network (0.947 ± 0.011) and mass action with the same form (0.941 ± 0.010) agree. Both match the prediction.
3. **Set B matches with the linear form** (1.000 ± 0.006). With the network it is 0.989 ± 0.005, as predicted (0.994). Exp 08's 1.007 ± 0.007 (40 reps), 1.9σ above the prediction, was noise.
4. **Exp 05's monthly shortfall is explained.** With 160 reps, `'start'` timing at monthly `dt` gives 0.980 ± 0.010, against 0.997 predicted (z = −1.7). Exp 05's 0.971 ± 0.022 (40 reps) is consistent with this.
5. **One residual remains: `k` = 0, weekly `dt`, `'within'` timing.** It's predicted at 0.995 for both forms, since the bias is small at weekly `dt`. Results:
   - exp 07, network: 0.959 ± 0.010 (160 reps);
   - here, mass action: 0.952 ± 0.017 (80 reps);
   - exp 07's numpy-timing variant: 0.985 ± 0.013 (80 reps).
   
   Pooled, that's 0.966 ± 0.007 (z ≈ −4). The same setting with `'start'` timing is 0.998 ± 0.008 (exp 07, pooled with its variant), as predicted. The residual doesn't involve the network. It only shows up near threshold at `k` = 0 with weekly `'within'` timing. Monthly `'within'` at `k` = 0 matches (0.981 ± 0.024 vs 0.979), and so does set A at weekly `dt` (0.987 ± 0.016 vs 0.986). The cause is unknown.

## Conclusions

- **The ABM↔ODE transmission gap is a time-step effect, and it's now quantified.**
  - Per-step transmission in any `TBExact` with network transmission, and equally in `tbsim.TB`, matches `1 − e^(−x)`.
  - That's a fraction of about x/2 below the ODE's infection flow: 0.7% of effective `beta` in set A at monthly `dt`, 0.4% in set B monthly, and under 0.1% at weekly `dt`.
  - It mainly matters near threshold, where prevalence moves 6–9× as much as `beta`.
  - Exp 05's statement that "transmission is exact" was right per step, but against the wrong benchmark: compared with the ODE's flow, the ABM transmits slightly too little.
- **A fix is available for mass-action transmission:** a per-step infection probability of x (`TBExactMassAction(infection_prob='linear')`). It still uses the susceptible pool at the start of the step, but its expected infections per step equal the ODE's flow, x·S. It matches the ODE to within about 1–2% in both sets at monthly `dt`.
  - **It isn't exact either, but what's left is small.** Even with the linear form, the deterministic step predicts a small shortfall, which shrinks with `dt`:

    | Linear mass action, `'within'` timing | Deterministic step ÷ ODE | `TBExact` ÷ ODE (mean ± SE, 190 reps) |
    |---|---|---|
    | Set A, monthly | 0.995 | 0.984 ± 0.010 (arm 3) |
    | Set B, monthly | 0.999 | 1.000 ± 0.006 (arm 5) |
    | Set A, weekly | 0.999 | not run |
    | Set B, weekly | 1.000 | 1.001 ± 0.006 ([exp 11](../11_incidence_counting/SUMMARY.md)) |

    Values are prevalence (ASY + SYM). Predictions are from `step_map.py` (`step_map.log`). The observed values agree with the predictions within noise. At ±1%, these reps can't resolve set A's predicted 0.5%.
  - **What causes the remainder hasn't been tested separately.** It is what the deterministic step reproduces once the per-step probability is right. The likely cause is that the force of infection is taken at its start-of-step value and held for the whole step, rather than updating within it. An exact match would need a force of infection that updates within the step, or a smaller `dt`.
  - With network transmission, there's no equally simple per-agent fix, because an agent's chance of infection depends on its number of infectious contacts. A per-step rescaling of `beta` by `−ln(1 − x)/x` would be approximate, because the correction differs for agents with different `rel_sus`.
- **For the corrected ABM:** tb_LAI_TPT's random network (`n_contacts=10`, `dur=0`) behaves like mass action (finding 2). So `TBExactMassAction` with the linear form is a candidate transmission model for `TBExact`. The household network, if it's ever included, would still need network transmission. Whether to adopt it is a decision for the user.
- **The `k` = 0 weekly `'within'` residual is still open.** It's small, confined to one near-threshold setting at `k` = 0, and doesn't affect the reference set B.

## Versions

tbsim main @ `b92d7e8` (v0.10.3), starsim @ `3521adb4` (v3.6.1). `TBExact` unchanged since exp 08 (reused seeds reproduced exactly). `TBExactMassAction` and `make_abm(network=...)` added for this experiment.
