# 12 Setting up the v4 track: summary

**Track: v4.** tbsim `origin/v4` @ `dcfb522`, starsim `origin/v4` @ `80c81f53`; references from tbsim v0.10.3 @ `b92d7e8` / starsim v3.6.1 @ `3521adb4` (`versions_v4.json`, `versions_v010.json`).

In tbsim v4, `tbsim.TB` is declared as states and transitions, so one model runs as an ODE or as agents. This experiment checks that the v4 ODE is the same model as our v0.10 ODE (`tbode.TBODE`), and how agents should be set up so that "v4 agents ÷ v4 ODE" is comparable with exps 01–11. The setting is the baseline of exps 03–05: tbsim's default natural history, no treatment, births only replacing background deaths, and `beta` calibrated to 500/100k.

## Findings

Finding 1 is about the v4 ODE (part A of the design). Findings 2–6 are about v4's agents (part B). They don't test whether agents match the ODE over a long run, which is exp 13. They test the settings exp 13 needs: the network, the per-contact `beta`, the start, and run time.

### A. The ODE

**1. The v4 ODE is our ODE.** v4 results can be set beside v0.10's directly.

| Check | Result | Criterion |
|---|---|---|
| A1: v4 ODE vs tbsim's `TB_ODE` (defaults, births replacing every death), 520 years | max relative difference 2.8e-5 over all states | 1e-4 |
| A2: v4 ODE vs `TBODE`, baseline, 1% seed, 1900–2100 monthly | 1.6e-6 relative from 1901; 9e-8 absolute throughout | 1e-4 / 1e-5 |
| A2: the same, set B's natural history at `k` = 0 | 3.2e-5 relative from 1901; 1.2e-6 absolute throughout | 1e-4 / 1e-5 |
| A3: started at `TBODE`'s equilibrium, 200 years | prevalence and every share constant to 3e-15 | — |
| A3: `beta` for 500/100k, calibrated on the v4 ODE alone | 9.917338, against `TBODE`'s 9.917338 (5e-11 relative) | — |

In the first month of A2, small compartments filling up from 0 differ by up to 3e-3 relative (2e-9 absolute). This comes from the first integration step: v4 uses RK4 with substeps, and `TBODE` uses `odeint`. Differences fall to ~1e-6 within a year. Prevalence agrees to 2e-7 throughout.

### B. Setting up agents

**2. `ss.RandomNet` is unchanged in the way that matters.** `n_contacts=10` gives 5.0 edges per agent and a mean degree of 10.0 (variance 5.0), as in starsim 3.6.1. No edge survives from one step to the next, so there are fresh contacts every step (default `dur=0`).

**3. As agents on a network, v4's own per-contact `beta` loses about 4% of infections at monthly `dt`; v0.10's matched mapping corrects it.** The ODE has one `beta`; agents on `RandomNet` need a `beta` per contact that gives the same transmission. To test each choice, the infections that actually happened among the agents in each step are compared with the number the ODE's formula predicts from the agents' own susceptible and infectious counts at the start of that step. This isolates transmission from any drift in the dynamics. Monthly `dt`, 50k agents, started at the equilibrium, 4 seeds × 50 years pooled; ± Poisson SE (~100k infections per arm):

| Transmission | ÷ ODE flow, x | ÷ `1 − e^(−x)` |
|---|---|---|
| Homogeneous mixing | 0.998 ± 0.003 | 1.000 ± 0.003 |
| `RandomNet`, `beta_ode / K` per contact (v4's mapping) | **0.956 ± 0.003** | 0.958 ± 0.003 |
| `RandomNet`, `beta_edge = −ln(1 − beta_ode·dt/K)/dt` (v0.10's, `tbode4.beta_edge_from_ode`) | 0.997 ± 0.003 | 0.999 ± 0.003 |

- **Why v4's mapping is low.** Each contact transmits with probability 1 − e^(−b·dt) in a step, so with b = `beta_ode/K` a susceptible's expected infections fall short of the ODE's by (1 − e^(−b·dt))/(b·dt). That is 0.960 at monthly `dt` with b = 0.99 per year. The measured 0.956 ± 0.003 agrees.
- **The ODE approximation assumes v4's mapping.** v4's ODE version of the same network sim (homogeneous mixing at the mean degree, `beta_ode/K` per contact) equals the plain ODE to 1e-15. So `solver='ode'` on a network sim overstates the agents' effective `beta` by about b·dt/2 per step. Near threshold that is large: a 4% change in `beta` at monthly `dt` is roughly a third of prevalence (exp 03).
- **The matched mapping is the one to use,** as in v0.10. It is `make_sim`'s default (`edge_beta='step'`). The x/2 shortfall of `1 − e^(−x)` (exp 10) is ~0.2% here, within the SE.

**4. Agents can start at the ODE equilibrium with `init=`.** Agents start with exactly the given shares: the largest difference is 0.4 of an agent, which is rounding. Exp 03's `init_from_ode` isn't needed while rates have no memory (no time in a state).

**5. Agents are 2.5–5× faster than in v0.10.** At 50k agents, a step takes 2.1 / 1.9 / 1.6 ms under homogeneous mixing and 4.0 / 3.7 / 3.3 ms on the network (monthly / weekly / daily `dt`). v0.10 took ~10.5 ms (exp 06). These are single runs on one core, after compiling the kernels.

**6. v4's `dt` notes and `check_dt()` aren't informative for this model.**
- **The notes are worst-case.** At monthly and weekly `dt`, v4 notes that reinfection after treatment "transmits at up to 31 per year (with everyone infectious)" and suggests a shorter step. That rate assumes everyone is infectious; prevalence here is 0.5%.
- **`check_dt()` is noise.** It compares one stochastic run with one other, using per-step peaks. Its 20–90% changes in each result are run-to-run noise at this population. Step effects need many reps, as in exp 13.

## What this means

- The v4 track is anchored. `v4/tests/test_anchor.py` keeps finding 1 (A1, A2 and the equilibrium) and the `beta_edge` round trip, to rerun whenever the v4 pins are bumped.
- Network runs on v4 must use the matched `beta_edge`. The 4% gap with v4's `beta_ode/K` is worth raising upstream, since v4 says a network run as an ODE is approximated at the network's mean degree.
- Exp 13 (the baseline, agents vs ODE at daily / weekly / monthly `dt`) can start from the equilibrium with `init=`.

## Files

- `export_reference.py` (v0.10 environment) → `reference_v010.json`: `TBODE`'s calibrated `beta`, baseline equilibrium, and trajectories for the baseline and set B at `k` = 0.
- `run.py` (v4 environment) → `outputs.obj` (all results, including per-run transmission sums and prevalence series), `summary.json`, `v4_setup.png`, `run.log`.
