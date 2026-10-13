# 10 Transmission within a step: what causes `TBExact`'s near-threshold shortfall?

## Question

`tbode.TBExact` is a version of `tbsim.TB` that samples each agent's end-of-step state from the exact one-step probabilities `expm(Q·dt)`. It has no step-ordering bias, uses Poisson-mixture latent classes for progression by time since infection, and dates each new infection within its step. Even so, in some near-threshold arms it runs a few percent below the ODE: about 5% in exp 08's set A at monthly `dt`. Exp 09 showed that shortfall doesn't depend on population size (25k–200k agents), so it's systematic, not a finite-population effect.

Exps 07–09 compared each arm with a deterministic version of the `TBExact` step (`step_map.py`): the expected update per step, iterated to equilibrium, with no sampling noise. That map assumed each susceptible's infection probability per step is linear, x = `beta`·w·`dt` (w = the `trans_asymp`-weighted infectious fraction). But on the random network, an agent with j infectious contacts is infected with probability `1 − (1 − q)^j`, which averages to `1 − e^(−x)`. Exp 05's flux check confirmed this per-step value.

The two forms differ in what they assume about the susceptible pool:
- **ODE:** at equilibrium, infections over a step are `beta·w·S·dt` exactly, because people who become susceptible during the step, for example by clearing a latent infection, can be infected later in the same step.
- **ABM:** the susceptible pool is fixed at the start of the step, and each agent can be infected at most once per step. So infections are `S·(1 − e^(−x))`, about a fraction x/2 fewer.

With `1 − e^(−x)` (applied per agent, with each agent's `rel_sus`), the deterministic map predicts most of the observed values (exps 05–09). One arm remains clearly off (exp 07 weekly `'within'`, −3.6σ), and set B monthly is 1.9σ high:

| Arm (`TBExact`) | Observed ÷ ODE (reps) | Map, linear | Map, `1 − e^(−x)` |
|---|---|---|---|
| Set A, monthly (exp 08–09) | 0.943 ± 0.005 (356) | 0.995 | 0.941 |
| Set A, weekly (exp 08) | 0.987 ± 0.016 (40) | 0.999 | 0.986 |
| Set B, monthly (exp 08) | 1.007 ± 0.007 (40) | 0.999 | 0.994 |
| Set B, weekly (exp 08) | 1.004 ± 0.009 (40) | 1.000 | 0.999 |
| `k` = 0, monthly, `'start'` (exp 05) | 0.971 ± 0.022 (40) | 1.014 | 0.997 |
| `k` = 0, monthly, `'within'` (exp 07) | 0.981 ± 0.024 (40) | 0.995 | 0.979 |
| `k` = 0, weekly, `'start'` (exp 07) | 1.004 ± 0.009 (160) | 1.004 | 1.000 |
| `k` = 0, weekly, `'within'` (exp 07) | 0.959 ± 0.010 (160) | 0.999 | **0.995** |

Set A is tbsim defaults + `k` = 5.846 (near threshold); set B is tb_LAI_TPT's natural history + `k` (the reference set); both are from exp 08. The `k` = 0 arms use tbsim defaults with the baseline `beta`.

This experiment tests the explanation directly:
1. Is the shortfall entirely due to the `1 − e^(−x)` saturation? If so, a mass-action ABM with per-step probability x should match the ODE (and the linear map).
2. Does the random network add anything? A mass-action ABM with `1 − e^(−x)` should reproduce the network runs if it doesn't.
3. Are the two remaining deviations real?

## Design

**Mass-action `TBExact`** (`tbode.TBExactMassAction`) replaces network transmission. Each step, every susceptible agent j is infected with probability:
- `'exp'`: `1 − exp(−beta·w·rel_sus_j·dt)`;
- `'linear'`: `min(1, beta·w·rel_sus_j·dt)`.

Here `beta` is the ODE value, w = (`trans_asymp`·ASY + SYM) / N over living agents, and `rel_sus` is the agent's susceptibility (the reinfection risk ratio for CLEARED, 0 for latent agents). There is no network (`make_abm(..., network=False)`). Everything else is unchanged: natural history, demography and `'within'` infection timing.

**Arms** (all `TBExact`, monthly `dt` unless noted, 100 years from the ODE equilibrium, outputs averaged over the last 50 years):

| # | Setting | Transmission | N | Reps | Expected if the explanation holds (map) |
|---|---|---|---|---|---|
| 1 | Set A | network | 25k | 190 (reused from exp 09) | 0.941 |
| 2 | Set A | mass action, `'exp'` | 25k | 190 | 0.941 |
| 3 | Set A | mass action, `'linear'` | 25k | 190 | 0.995 |
| 4 | Set B | network | 25k | 190 | 0.994 |
| 5 | Set B | mass action, `'linear'` | 25k | 190 | 0.999 |
| 6 | `k` = 0, `'start'` | network | 50k | 160 (seeds 0–39 from exp 05) | 0.997 |
| 7 | `k` = 0, weekly, `'within'` | mass action, `'exp'` | 50k | 80 | 0.995 |

- **Arms 1–3** are the main test of questions 1 and 2. 25k agents is enough, since exp 09 showed no dependence on N for set A.
- **Arms 4–5** recheck set B's +1.9σ deviation with more reps, and test the linear form far from threshold.
- **Arm 6** repeats exp 05's monthly `'start'` arm with 4× the reps (exp 05: 0.971 ± 0.022). It uses 50k agents, because at baseline `beta` 25k agents would risk the epidemic dying out.
- **Arm 7** asks whether exp 07's weekly `'within'` shortfall (0.959 ± 0.010, against 0.995 predicted) also appears without the network.
- **Reps** give an SE of about 1% in each arm, from the rep-to-rep spreads in exps 05, 07 and 09.
- **Reused runs:** seed 0 of each reused arm is rerun first and checked to reproduce the saved series exactly.

**Deterministic map** (`step_map.py`): exp 08's map extended with a choice of `'linear'` or `'exp'` per-step infection probability, covering both `k` = 0 and sets A and B. It gives the expected values in the last column. Output: `step_map.log`.

**Comparison:** ABM ÷ ODE as mean ± SE across reps; the difference from each map prediction, as a z-score; the number of fade-outs.

**Outputs:**
- per-rep time series (`outputs.obj`) and `summary.csv`;
- a figure of each arm's mean ± SE against its predictions (`transmission_step.png`);
- prevalence and incidence trajectories per arm against the ODE (`trajectories.png`, from `plot_trajectories.py`; plotted series in `trajectories.csv`), added afterwards for context.

## Versions

tbsim main @ `b92d7e8` (v0.10.3), starsim @ `3521adb4` (v3.6.1). `TBExact` is unchanged since exp 08; `TBExactMassAction` and `make_abm(network=...)` are new in `tbode`.
