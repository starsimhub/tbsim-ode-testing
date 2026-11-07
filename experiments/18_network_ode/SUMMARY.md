# 18: Summary

**Track:** v4, in this experiment's own environment: starsim `origin/v4` @ `bc2704b2` (adds `5ba36ebb`, "add beta correction for ODEs on networks") and tbsim `origin/v4` @ `efeca56` (`versions.json`). `v4/tests/test_anchor.py` passes in it: the mass-action v4 ODE still reproduces `TB_ODE` and `TBODE`, and calibrates to the same `beta` as before (set A 9.926, set B two-state 21.143). starsim's newer commits (`308a27a5` on) don't yet work with tbsim `efeca56`: they reserve the name `transmission`, which tbsim uses.

**Setting:** `ss.RandomNet(n_contacts=10, dur=0)`, fresh contacts each step. Set A is tbsim's defaults (near threshold; reinfection of cleared people at the susceptibles' rate). Set B is exp 17's two-state model at φ = 0.01, with every reinfection risk at π = 0.24 and latent reinfection on. In set B most infections are reinfections.

**The fix:** Starsim's ODE of a network now multiplies each contact's transmission by (1 − e^(−b·D))/(b·D), the share of mass action that agents achieve when each contact transmits at most once in a step. Here b is the per-contact rate and D the step.

## Finding 1: the network ODE is the same model as our convention

Calibrated to 500/100k on the network ODE, the per-contact b equals our matched value, −ln(1 − `beta`·dt/K)/dt from the mass-action `beta`, to five significant figures. That holds for both sets at monthly and weekly steps (`run_q3.log`). For set B, b is 2.326 monthly and 2.158 weekly, against `beta`/K = 2.114. So the network ODE behaves exactly like the mass-action ODE with an effective `beta` of K·(1 − e^(−b·D))/D. Agents get the same per-contact `beta` either way.

**So Starsim now does internally what our convention did by hand.** Calibrating and running the ODE on the network, as Cliff suggests, gives the same model and the same agents. No mapping is needed, except to carry over a `beta` calibrated with mass action.

## Finding 2: neither handles reinfection at a reduced risk

The ODE applies the susceptibles' factor to every transition by contact. Each reinfection transition (cleared, recovered, and remote → recent latent) runs at rr·b on the same contacts. As agents, each contact then transmits with probability 1 − e^(−rr·b·D), so the share of mass action is closer to 1 for rr < 1.

| Agents ÷ ODE, per person over 24 steps (1M agents, 4 seeds; mean ± SE) | Set B, monthly | Predicted | Set B, weekly | Predicted |
|---|---|---|---|---|
| Transmission (susceptibles, rr = 1) | 0.994 ± 0.005 | 1.000 | 0.995 ± 0.007 | 1.000 |
| Reinfection of cleared (rr = π) | **1.078 ± 0.003** | 1.075 | 1.011 ± 0.007 | 1.016 |
| Reinfection of recovered (rr = π) | **1.076 ± 0.008** | 1.075 | 1.021 ± 0.008 | 1.016 |
| Remote → recent (rr = π) | **1.071 ± 0.010** | 1.075 | 0.989 ± 0.016 | 1.016 |

The prediction is agents' factor for each transition's rr over the ODE's (the susceptibles'). Set A's reinfection of recovered people (rr = 0.21) is 1.020 ± 0.010 monthly (predicted 1.034). Its reinfection of cleared people (rr = 1) matches. The two workflows give identical numbers (`q1_summary.csv`).

## Finding 3: at equilibrium, set B's agents run 20% above the ODE at monthly steps

| Prevalence, agents ÷ ODE, 1950–2000 (50k agents, 40 reps; mean ± SE) | Observed | Predicted from finding 2 |
|---|---|---|
| Set A, monthly | 1.008 ± 0.018 | 1.023 |
| Set A, weekly | 1.017 ± 0.018 | 1.005 |
| Set B, monthly | **1.204 ± 0.010** | 1.201 |
| Set B, weekly | **1.052 ± 0.009** | 1.037 |

- **How the prediction is made:** the mass-action ODE, with each reinfection risk scaled by agents' factor over the ODE's (`explain.json`). It accounts for the whole gap, so there's nothing else wrong.
- **Other outputs:** latent infection and TB deaths move with prevalence (set B monthly 1.22 and 1.20). The cleared share is within 0.2% (`q2_summary.csv`).
- **Why set B is so sensitive:** most of its infections are reinfections at π, and raising every reinfection risk by 7.5% raises prevalence by 20%.
- **The same holds under our convention:** agents on the network with the matched `beta` would run 20% above the mass-action ODE in set B at monthly steps. This wasn't seen before because exps 12–13's network runs used set A, where reinfection of cleared people is at the susceptibles' rate. Exps 14–17 used homogeneous mixing, which this doesn't affect.

## What this means

- **The override itself is no longer needed.** Starsim's network ODE gives the same model as our matched `beta`. The convention can become "run and calibrate the ODE on the network; use `beta_edge_from_ode` only to carry over a mass-action `beta`".
- **But neither version makes agents on the network match the ODE for a model like set B,** where reinfection runs at a reduced risk. The gap is 20% of prevalence at monthly steps and 5% at weekly. Ways to close it:
  - **Fix in Starsim:** apply each transition's own factor, worked out from its own rate (rr·b), not the first transition's ([proposal 03](../../proposals/03_network_reinfection/proposal.md)). The factor already varies by group with susceptibility, so tbsim expressing reinfection risk as susceptibility might also work; not tested.
  - **Use homogeneous mixing for v4 agents,** as in exps 13–17, where agents match the ODE within about 1%. A `RandomNet` with fresh contacts every step is mass action plus this per-step saturation, so it adds no epidemiology. tb_LAI_TPT's household network, if ever included, would still need a network.
  - **Correct the ODE's reinfection risks ourselves,** as the prediction does: it matched agents (1.201 against 1.204).
  - **A shorter step** shrinks the gap (5% weekly) but doesn't remove it.

**Decided (user, 2026-11-06):** the ODE recalibration uses homogeneous mixing, so none of this affects it.

## Not covered

- **Networks whose edges last several steps** (`dur` > 0), such as households.
- **Moving the v4 track's pins** to these commits: the anchor test passes, but every v4 result would need checking against the new indexing of results by step.
- **Starsim's newest commits,** which need a tbsim update first.
