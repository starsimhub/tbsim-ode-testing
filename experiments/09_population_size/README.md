# 09 Population size: is the near-threshold `TBExact` shortfall a finite-population effect?

## Question

`tbode.TBExact` is a version of `tbsim.TB` that samples each agent's end-of-step state from the exact one-step probabilities of the continuous-time model, `expm(Q·dt)`. It has no step-ordering bias, uses Poisson-mixture latent classes for progression by time since infection, and dates each new infection within its step. Its runs normally match the ODE. But in two near-threshold arms, at 50,000 agents, the ABM's mean prevalence fell a few percent below what a deterministic version of the same step predicts:
- **Exp 08, set A (tbsim defaults + `k` = 5.846), monthly `dt`:** 0.937 ± 0.018 of the ODE, against a deterministic prediction of 0.995 (`step_map.py`). That's 3.2σ below.
- **Exp 07, baseline (`k` = 0), weekly `dt`, `'within'` timing:** 0.959 ± 0.010, against 0.999 predicted.

Further from threshold (exp 08, set B), there's no shortfall (1.00–1.01 ± 0.01).

The deterministic step applies the expected force of infection and expected transitions. It leaves out the randomness of a finite population. Near threshold, random fluctuations in prevalence are large and slow to decay. Because the dynamics are nonlinear (infections ∝ infectious × susceptible), these fluctuations can hold the stochastic mean below the deterministic equilibrium. Such a finite-population effect would shrink as the population grows, roughly like 1/N, and wouldn't depend on `dt`.

**Does the shortfall shrink with more agents?** If it does, roughly in proportion to 1/N, it's a finite-population effect: a property of any 50k-agent ABM near threshold, not an inconsistency between the ODE and the ABM's model. If it doesn't, something else in `TBExact` is responsible and needs finding.

## Design

- **Setup:** exp 08, set A, monthly `dt`, `TBExact(infection_timing='within')`. Parameters and the ODE `beta` (36.753, for 500/100k) come from `../08_progression_tsi/calibration.json`. The ABM starts at the ODE equilibrium, with latent agents' classes and times since infection drawn as in exp 08. It runs for 100 years (1900–2000), using `ss.RandomNet(n_contacts=10, dur=0)` and constant births and deaths. This arm was chosen because it has the largest shortfall and monthly runs are cheap.
- **Arms:** population N = 25k, 50k, 100k and 200k agents.
- **Reps:** chosen so each arm's SE is about 1% of prevalence. The rep-to-rep spread at 50k is about 10% (SD), and should scale as 1/√N.

  | N | Reps | Seeds | Note |
  |---|---|---|---|
  | 25k | 190 | 0–189 | |
  | 50k | 94 | 0–93 | Seeds 0–39 reloaded from exp 08; the script first reruns seed 0 and checks it reproduces exp 08 exactly |
  | 100k | 48 | 0–47 | |
  | 200k | 24 | 0–23 | |

- **Comparison:** each rep's mean prevalence over the last 50 years. For each N, report ABM ÷ ODE as mean ± SE, and the shortfall relative to the deterministic step (0.995). Fit shortfall = c/N and shortfall = c/√N, and check which fits.
- **Fade-out:** count reps in which prevalence reaches zero. There's no importation, so the epidemic can die out in small populations; with 125 active cases expected at 25k, this should be rare. Report means with and without such reps.
- **Outputs:** prevalence, INFECTION, incidence, TB mortality (per-rep time series, `outputs.obj`); table in `summary.csv`; figure `population_size.png` (ABM ÷ ODE against 1/N, and the time course by N).
- **Predictions:**
  - **Finite-population effect:** shortfall of about 12% / 6% / 3% / 1.5% at 25k / 50k / 100k / 200k (∝ 1/N).
  - **Not finite-population:** about 6% at every N.

## Versions

tbsim main @ `b92d7e8` (v0.10.3), starsim @ `3521adb4` (v3.6.1). `TBExact` is unchanged since exp 08.
