# 07 Infection timing within the step: summary

## Background

`tbode.TBExact` is a version of `tbsim.TB` without tbsim's step-ordering bias. Instead of applying transitions one state at a time, it draws each agent's end-of-step state from the exact one-step probabilities of the continuous-time model, `expm(Q·dt)`. In exp 05 it matched the ODE, apart from INFECTION reading 0.89 of the ODE at monthly `dt`. The explanation was that agents infected during a step went through the whole step's transitions, as if infected at its start (`'start'` timing).

This experiment adds `infection_timing='within'`. Each new infection gets a uniform time within its step and goes through only the rest of the step, `expm(Q·(1 − u)·dt)`. The baseline setup is unchanged: tbsim default natural history, `beta` for 500/100k, 50k agents, 100 years from the ODE equilibrium, outputs averaged over the last 50 years.

## Results

All values are ABM ÷ ODE equilibrium (1.00 = match), mean ± SE across reps.

| Output | Monthly, `'start'` (40) | Monthly, `'within'` (40) | Weekly, `'start'` (160) | Weekly, `'within'` (160) |
|---|---|---|---|---|
| INFECTION | **0.89 ± 0.02** | **0.98 ± 0.02** | 0.98 ± 0.01 | 0.96 ± 0.01 |
| INFECTION ÷ prevalence ratio | 0.91 | 1.00 | 0.98 | 1.00 |
| Prevalence (ASY + SYM) | 0.97 ± 0.02 | 0.98 ± 0.02 | 1.00 ± 0.01 | **0.96 ± 0.01** |
| Incidence | 0.94 ± 0.02 | 0.95 ± 0.02 | 1.00 ± 0.01 | 0.95 ± 0.01 |
| TB mortality | 0.97 ± 0.02 | 0.98 ± 0.02 | 1.00 ± 0.01 | 0.96 ± 0.01 |
| CLEARED | 0.99 ± 0.00 | 0.99 ± 0.00 | 1.00 ± 0.00 | 0.99 ± 0.00 |
| % asymptomatic | 1.00 | 1.00 | 1.00 | 1.00 |

The "INFECTION ÷ prevalence ratio" row divides the INFECTION ratio by the prevalence ratio. It shows the INFECTION artefact apart from the overall level of transmission in each arm.

**The INFECTION artefact is gone, as predicted.** At monthly `dt`, INFECTION rises from 0.89 to 0.98 of the ODE, a paired difference of +0.09 ± 0.03 (z = 3.5). Relative to each arm's level of transmission, INFECTION goes from 0.91 to 1.00 at monthly `dt` and from 0.98 to 1.00 at weekly. This confirms exp 05's explanation: dating infections to the start of the step shifted each cohort half a step further into its time in INFECTION. The new sampling is also checked directly:
- **Unit test (`tests/test_exact.py`):** in one monthly step, new infections end in each state with the probabilities `∫₀¹ expm(Q·s·dt) ds`, computed independently (p > 0.03 in all checks). A control test shows `'start'` fails the same check.
- **Inside running sims (`run_flux_diag.py`, weekly, 12 seeds per timing):**
  - New infections' end-of-step states match the uniform-time average.
  - Infections from SUSCEPTIBLE are 1.003 ± 0.003 (`'within'`) and 0.999 ± 0.003 (`'start'`) of the ODE expectation.
  - Every natural-history flow agrees between the two timings to about 0.1%.

**At monthly `dt`, nothing else changes,** as predicted. Prevalence (0.98 vs 0.97), incidence (0.95 vs 0.94) and mortality all agree within noise. So the monthly incidence shortfall (about 0.95) isn't caused by infection timing; it's brief visits to ASYMPTOMATIC missed between end-of-step snapshots, as [exp 11](../11_incidence_counting/SUMMARY.md) confirmed.

**At weekly `dt`, `'within'` is unexpectedly low, and this isn't fully explained.** Prevalence is 0.959 ± 0.010 for `'within'` and 1.004 ± 0.009 for `'start'`, with 160 reps each. Every output moves together, which is the signature of slightly less transmission: about 0.4% of effective `beta`, given that prevalence moves about 9× as much as `beta` near threshold. The shortfall is consistent across seed batches: 0.948 / 0.966 / 0.960 for seeds 0–39 / 40–79 / 80–159. Several things argue against a real timing effect, or at least one this large:
- **Timing alone shouldn't change the equilibrium much.** A deterministic version of the step (`step_map.py`) applies the start-of-step force of infection and then exact transitions, with no sampling noise. With the ABM's per-step infection probability `1 − e^(−x)` (corrected in [exp 10](../10_transmission_step/SUMMARY.md); originally the linear x was used), it predicts prevalence ÷ ODE of 0.995 / 0.979 for `'within'` and 1.000 / 0.997 for `'start'`, at weekly / monthly `dt`. So weekly `'within'` should be within 0.5% of the ODE.
- **Any timing effect should grow with `dt`.** Monthly would then show about 4× the weekly gap, but shows +0.01 ± 0.03.
- **Two diagnostic variants of `'within'` (`run_variants.py`, 80 weekly reps each):**
  - `u0` draws infection times from the same random stream but ignores them, so its timing is the same as `'start'`.
  - `np` applies the times but draws them from an independent numpy generator.
  - They give 0.981 ± 0.015 (`u0`) and 0.985 ± 0.013 (`np`), the same within noise (difference +0.004 ± 0.022). So `'within'` timing doesn't lower prevalence relative to `'start'` timing here.
- **Pooling runs by timing** gives 0.998 ± 0.008 for `'start'`-like runs and 0.969 ± 0.008 for `'within'`-like runs. That's a difference of −0.029 ± 0.011, against −0.005 predicted, leaving 0.024 ± 0.011 (z ≈ 2.2) unexplained. Exp 10 added 80 mass-action reps of weekly `'within'` (0.952 ± 0.017). Pooled over all its variants, weekly `'within'` is 0.966 ± 0.007, against 0.995 predicted (z ≈ −4).

**Interpretation:** the weekly `'within'` arm is 4% below the deterministic prediction (z ≈ 4), and no mechanism or bug has been found. It isn't specifically a timing effect: the two variants agree with each other, and an effect of timing would grow with `dt`. It doesn't involve the network: mass-action runs show it too (exp 10). It only appears near threshold, at `k` = 0, with weekly `'within'` timing. Monthly `'within'` (0.981 ± 0.024 vs 0.979) and set A weekly (exp 08: 0.987 ± 0.016 vs 0.986) match. The cause is still unknown.

**Pairing by seed doesn't reduce noise.** Runs with the same seed and different timing are uncorrelated (r ≈ 0) after 100 years, so paired differences are as noisy as unpaired ones.

## Conclusions

- **Use `'within'` timing in `TBExact` from now on.** It's exact by construction for a force of infection that is constant over the step, it passes the unit and in-sim checks, and it removes the INFECTION artefact at every `dt`. It matters most with progression by time since infection: with `k` = 5.85/yr, the progression hazard falls by a factor of e about every 2 months, and `'start'` would make the fastest progressors progress about half a step early.
- **Other near-threshold shortfalls are explained, but not this one.** The ABM's per-step infection probability `1 − e^(−x)` transmits a fraction of about x/2 less than the ODE ([exp 10](../10_transmission_step/SUMMARY.md)). That accounts for the monthly shortfalls (exp 05's 0.97, exp 08's set A 0.94), but predicts only 0.5% at weekly `dt`. The 4% weekly `'within'` residual remains unexplained.

## Versions

tbsim main @ `b92d7e8` (v0.10.3), starsim @ `3521adb4` (v3.6.1). `../tbsim` had an uncommitted change, but only to `tbsim/compartmental/README.md` (documentation, not code).
