# 07 Infection timing within the step: summary

## Background

`tbode.TBExact` is a version of `tbsim.TB` without tbsim's step-ordering bias. Instead of applying transitions one state at a time, it draws each agent's end-of-step state from the exact one-step probabilities of the continuous-time model, `expm(Q·dt)`. In exp 05 it matched the ODE, apart from INFECTION reading 0.89 of the ODE at monthly `dt`. The explanation was that agents infected during a step went through the whole step's transitions, as if infected at its start (`'start'` timing).

This experiment adds `infection_timing='within'`. Each new infection gets a uniform time within its step and goes through only the rest of the step, `expm(Q·(1 − u)·dt)`. The baseline setup is unchanged: tbsim default natural history, `beta` for 500/100k, 50k agents, 100 years from the ODE equilibrium, outputs averaged over the last 50 years.

**Two separate step effects in transmission.** Both come from splitting transmission into steps, and both shrink with `dt`, but they are different:

- **Timing (this experiment).** With `'start'` timing, each infection is dated to the start of its step, rather than spread through it as in the ODE. This changes *when* infections happen, not how many. Its symptom is a low INFECTION stock, and `'within'` timing fixes it.
- **Amount ([exp 10](../10_transmission_step/SUMMARY.md), found later).** Each susceptible agent's infection probability per step is `1 − e^(−x)`, where x = `beta`·w·`rel_sus`·`dt` and w is the `trans_asymp`-weighted infectious fraction. The ODE's flow corresponds to x. The ABM's susceptible pool is fixed at the start of the step and each agent can be infected at most once per step, so the ABM has a fraction of about x/2 fewer infections. This changes *how many* infections happen, and lowers prevalence near threshold. It is unaffected by `'start'` vs `'within'` timing.

This experiment didn't test the second effect. It is mentioned below only where it changes what the results here are compared with.

## Results

All values are ABM ÷ ODE equilibrium (1.00 = match), mean ± SE across reps.

| Output | Monthly, `'start'` (40) | Monthly, `'within'` (40) | Weekly, `'start'` (160) | Weekly, `'within'` (160) |
|---|---|---|---|---|
| INFECTION | **0.89 ± 0.02** | **0.98 ± 0.02** | 0.98 ± 0.01 | 0.96 ± 0.01 |
| INFECTION ÷ prevalence ratio | 0.91 | 1.00 | 0.98 | 1.00 |
| Prevalence (ASY + SYM) | 0.97 ± 0.02 | 0.98 ± 0.02 | 1.00 ± 0.01 | **0.96 ± 0.01** |
| Incidence (snapshot count) | 0.94 ± 0.02 | 0.95 ± 0.02 | 1.00 ± 0.01 | 0.95 ± 0.01 |
| TB mortality | 0.97 ± 0.02 | 0.98 ± 0.02 | 1.00 ± 0.01 | 0.96 ± 0.01 |
| CLEARED | 0.99 ± 0.00 | 0.99 ± 0.00 | 1.00 ± 0.00 | 0.99 ± 0.00 |
| % asymptomatic | 1.00 | 1.00 | 1.00 | 1.00 |

The "INFECTION ÷ prevalence ratio" row divides the INFECTION ratio by the prevalence ratio. It shows the INFECTION artefact apart from the overall level of transmission in each arm.

### The INFECTION artefact is gone, as predicted

At monthly `dt`, INFECTION rises from 0.89 to 0.98 of the ODE, a paired difference of +0.09 ± 0.03 (z = 3.5). Relative to each arm's level of transmission, INFECTION goes from 0.91 to 1.00 at monthly `dt` and from 0.98 to 1.00 at weekly. This confirms exp 05's explanation: dating infections to the start of the step shifted each cohort half a step further into its time in INFECTION. The new sampling is also checked directly:
- **Unit test (`tests/test_exact.py`):** in one monthly step, new infections end in each state with the probabilities `∫₀¹ expm(Q·s·dt) ds`, computed independently (p > 0.03 in all checks). A control test shows `'start'` fails the same check.
- **Inside running sims (`run_flux_diag.py`, weekly, 12 seeds per timing):**
  - New infections' end-of-step states match the uniform-time average.
  - Infections from SUSCEPTIBLE are 1.003 ± 0.003 (`'within'`) and 0.999 ± 0.003 (`'start'`) of the expected per-step value.
  - Every natural-history flow agrees between the two timings to about 0.1%.

### At monthly `dt`, nothing else changes

As predicted, prevalence (0.98 vs 0.97), incidence (0.95 vs 0.94) and mortality all agree within noise between the two timings. So infection timing doesn't cause the remaining monthly shortfalls. Both were explained later (see the last section).

### At weekly `dt`, `'within'` runs low, and this is unexplained

Prevalence is 0.959 ± 0.010 for `'within'` and 1.004 ± 0.009 for `'start'`, with 160 reps each. Every output moves together, which is the signature of slightly less transmission: about 0.4% of effective `beta`, given that prevalence moves about 9× as much as `beta` near threshold. The shortfall is consistent across seed batches: 0.948 / 0.966 / 0.960 for seeds 0–39 / 40–79 / 80–159. Several things argue against a real timing effect, or at least one this large:

- **A deterministic version of the step predicts almost no shortfall.** It applies the start-of-step force of infection and then exact transitions, with no sampling noise, and is iterated to equilibrium. There are two versions:
  - This experiment's `step_map.py` (`step_map.log`) uses a linear per-step infection probability, x. It predicts 0.999 for weekly `'within'`.
  - Exp 10's `step_map.py` adds the ABM's actual probability, `1 − e^(−x)`. It predicts 0.995.

  | Prevalence ÷ ODE, deterministic step | Weekly `'start'` | Weekly `'within'` | Monthly `'start'` | Monthly `'within'` |
  |---|---|---|---|---|
  | Linear x (this experiment's map) | 1.004 | 0.999 | 1.014 | 0.995 |
  | `1 − e^(−x)`, as in the ABM (exp 10's map) | 1.000 | **0.995** | 0.997 | 0.979 |

  The comparisons below use the second row. So the x/2 effect accounts for only about 0.5% of the weekly `'within'` arm's 4% shortfall.
- **Any timing effect should grow with `dt`.** Monthly would then show about 4× the weekly gap, but `'within'` − `'start'` at monthly `dt` is +0.01 ± 0.03.
- **Two diagnostic variants of `'within'` (`run_variants.py`, 80 weekly reps each):**
  - `u0` draws infection times from the same random stream but ignores them, so its timing is the same as `'start'`.
  - `np` applies the times but draws them from an independent numpy generator.
  - They give 0.981 ± 0.015 (`u0`) and 0.985 ± 0.013 (`np`), the same within noise (difference +0.004 ± 0.022). So `'within'` timing doesn't lower prevalence relative to `'start'` timing here.
- **Pooling runs by timing** gives 0.998 ± 0.008 for `'start'`-like runs (`'start'` and `u0`) and 0.969 ± 0.008 for `'within'`-like runs (`'within'` and `np`). That's a difference of −0.029 ± 0.011, against −0.005 predicted, leaving 0.024 ± 0.011 (z ≈ 2.2) unexplained.

**Interpretation:** the weekly `'within'` arm is about 4% below the deterministic prediction, and no mechanism or bug has been found. It doesn't look specifically like a timing effect: the two variants agree with each other, and a timing effect would grow with `dt`. Settings that match their predictions are listed in the last section.

**Pairing by seed doesn't reduce noise.** Runs with the same seed and different timing are uncorrelated (r ≈ 0) after 100 years, so paired differences are as noisy as unpaired ones.

## Conclusions

- **Use `'within'` timing in `TBExact` from now on.** It's exact by construction for a force of infection that is constant over the step, it passes the unit and in-sim checks, and it removes the INFECTION artefact at every `dt`. `'start'` timing shifts each infected agent's whole trajectory half a step earlier. In `TBExact`, the latent class (which carries progression by time since infection) shifts with it, so durations and outcomes, including how many progress, are unchanged. What changes is when the infection is dated, and so the generation interval (from one infection to the infections it causes). Because the force of infection is taken from the start-of-step state and held for the whole step, the two timings trade one error for another:
- **`'start'`** dates each infection at the same moment as the infectiousness that caused it, so the mean generation interval is correct (to second order in `dt`). But every infection sits at a step boundary, so end-of-step snapshots of INFECTION read low.
- **`'within'`** spreads infections through the step, so snapshots are correct. But each infection happens on average half a step after the infectiousness that caused it, so the mean generation interval is half a step too long, and the lag adds up over successive generations.

A toy calculation (`dt` = 1 month, exponential latent and infectious periods) confirms this: the mean generation interval is +0.01 months from the continuous-time value with `'start'` and +0.50 months with `'within'`. Neither affects the equilibrium. For runs that change over time, `'within'` lags by about half a step per generation; this hasn't been measured in the ABM. Getting both right would need infections spread through the step *and* a force of infection centred on the step rather than taken at its start. It matters for dynamics that change over time, not for the equilibrium.
- **Weekly `'within'` timing at `k` = 0 runs about 3–4% low, for a reason not yet found.** This doesn't affect the choice of `'within'` timing.

## What later experiments resolved, and what is still open

| Result here | Status | Explanation | Where |
|---|---|---|---|
| Monthly INFECTION 0.89 with `'start'` timing | **Resolved here** | Infections dated half a step early; fixed by `'within'` timing | This experiment |
| Monthly prevalence and mortality 0.97–0.98, with either timing | **Resolved** | The x/2 effect (`1 − e^(−x)` infection probability), predicted at 0.997 (`'start'`) and 0.979 (`'within'`). Rerun with 160 reps (`'start'`): 0.980 ± 0.010. A linear infection probability removes most of it in other near-threshold arms | [Exp 10](../10_transmission_step/SUMMARY.md) |
| Monthly incidence about 0.95, with either timing | **Resolved** | Snapshot artefact: brief visits to ASYMPTOMATIC between end-of-step snapshots are missed (3.3% of episodes at monthly `dt`). Counting from transitions matches the ODE | [Exp 11](../11_incidence_counting/SUMMARY.md) |
| Weekly `'within'` prevalence 0.96 | **Open** | See below | Exps 07, 10 |

**Weekly `'within'` at `k` = 0, where it stands:**
- **More runs:** exp 10 added 80 mass-action reps (0.952 ± 0.017). Pooled over all variants (network, mass action, `np`), weekly `'within'` is 0.966 ± 0.007, against 0.995 predicted (z ≈ −4).
- **Ruled out:**
  - per-step flows and the timing draw itself (this experiment);
  - the source of the random numbers (the `np` variant);
  - the random network: mass action shows the same shortfall (exp 10);
  - the per-step infection probability: it predicts only 0.5% at weekly `dt` (exp 10);
  - `'within'`'s generation-interval lag (see Conclusions): the deterministic step includes it, and a lag doesn't change the equilibrium. It would also be 4× larger at monthly `dt`, where `'within'` matches.
- **Not tested for this arm:** population size. [Exp 09](../09_population_size/SUMMARY.md) found no population-size effect in a different near-threshold arm (set A, monthly).
- **Where it appears:** only near threshold, at `k` = 0, with weekly `'within'` timing. These settings match their predictions:
  - weekly `'start'` (0.998 ± 0.008, pooled with `u0`, vs 1.000);
  - monthly `'within'` (0.981 ± 0.024 vs 0.979);
  - set A weekly ([exp 08](../08_progression_tsi/SUMMARY.md): 0.987 ± 0.016 vs 0.986).

  It doesn't affect set B, the reference set.

## Versions

tbsim main @ `b92d7e8` (v0.10.3), starsim @ `3521adb4` (v3.6.1). `../tbsim` had an uncommitted change, but only to `tbsim/compartmental/README.md` (documentation, not code).
