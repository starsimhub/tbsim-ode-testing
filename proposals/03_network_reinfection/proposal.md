# Proposal 03: a saturation factor for each transition in Starsim's ODE of a network

**For:** Starsim (Cliff). A follow-up to `5ba36ebb`, "add beta correction for ODEs on networks".
**Status:** draft for the user to review before sending.
**Evidence:** [exp 18](../../experiments/18_network_ode/SUMMARY.md), starsim `bc2704b2` with tbsim `efeca56`.

## Summary

- **What the fix gets right:** for transmission at the disease's own per-contact rate, the ODE of a network now matches agents. Calibrated on the network ODE, the per-contact `beta` comes out exactly where agents need it.
- **What it misses:** the saturation factor is worked out once per exposure, from the first transition by contact that reads it. It's then applied to every transition on those contacts. A transition at a different per-contact rate, such as reinfection at a reduced risk rr·b, saturates differently as agents. The ODE gets it wrong by the ratio of the two factors.
- **The size of it:** for `tbsim.TB` with tb_LAI_TPT's calibrated natural history, on `ss.RandomNet(n_contacts=10, dur=0)` at monthly steps, agents' prevalence is 1.20× the network ODE's (1.05× at weekly steps).
- **Who else it affects:** any declared model with more than one transition by contact at different rates, e.g. reinfection with partial immunity, or vaccinated people infected at a reduced rate.

## The mechanism

As agents, a contact transmits at most once in a step. Over a step D, a contact at per-contact rate x/D transmits with probability 1 − e^(−x), a share (1 − e^(−x))/x of mass action. That share depends on the transition's own rate.

| Set B, monthly (b·D = 0.194) | Susceptible (rate b) | Cleared, reinfected at π·b (π = 0.24) |
|---|---|---|
| Agents: chance per infectious contact per step | 1 − e^(−0.194) = 0.176 | 1 − e^(−0.0465) = 0.0454 |
| ODE: the susceptibles' factor applied to both | 0.176 | 0.24 × 0.176 = 0.0423 |
| Agents ÷ ODE | 1.00 | **1.075** |

In `tbsim.TB`, the transmission transitions read the same exposure, β·(trans_asymp·asymptomatic + symptomatic):

- `transmission`: susceptible → latent, at rate 1;
- `reinfect_cle`, `reinfect_rec`, `reinfect_trd`: cleared, recovered and treated → latent, at `rr_reinfection_*`;
- `reinfect_remote`: remote latent → recent, at `rr_reinfection_inf`.

`Compartmental.exposure_flow()` picks the first of these, `transmission`. `saturation()` works out its factor, and `set_contacts()` multiplies the exposure's contacts by it. So every reinfection transition gets the susceptibles' factor.

## Evidence (exp 18)

**Per transition:** 1,000,000 agents × 4 seeds, 24 steps from equilibrium. Agents ÷ network ODE, flow per person, mean ± SE; the prediction is (1 − e^(−rr·b·D))/(rr·b·D) ÷ (1 − e^(−b·D))/(b·D).

| tb_LAI_TPT's natural history, rr = 0.24 | Monthly | Predicted | Weekly | Predicted |
|---|---|---|---|---|
| `transmission` (rr = 1) | 0.994 ± 0.005 | 1.000 | 0.995 ± 0.007 | 1.000 |
| `reinfect_cle` | 1.078 ± 0.003 | 1.075 | 1.011 ± 0.007 | 1.016 |
| `reinfect_rec` | 1.076 ± 0.008 | 1.075 | 1.021 ± 0.008 | 1.016 |
| `reinfect_remote` | 1.071 ± 0.010 | 1.075 | 0.989 ± 0.016 | 1.016 |

**At equilibrium:** 50,000 agents, 40 reps, 1950–2000, prevalence agents ÷ network ODE: 1.204 ± 0.010 monthly and 1.052 ± 0.009 weekly. With tbsim's default rates, where reinfection of cleared people runs at rr = 1, the results are 1.008 ± 0.018 and 1.017 ± 0.018.

**The cause accounts for all of it:** the mass-action ODE with each `rr_reinfection_*` scaled by the ratio of the agents' factor to the ODE's predicts 1.201 (monthly) and 1.037 (weekly).

## Suggested change

Work out the factor for each transition by contact, from that transition's own per-contact rate:

- **The rate:** its rate parameters (here `rr_reinfection_*`), its modifiers (`flow_mult`), the exposure's transmissibility and β. That gives x_f = b_f·D and the factor (1 − e^(−x_f))/x_f.
- **Where to apply it:** contacts are stored per exposure and term, so either apply each transition's factor relative to the exposure's as a multiplier on that transition's rate, or store contacts per transition.
- **Smaller, optional:** the factor uses the mean transmissibility of the source states (`exposure_trans`). As agents, a contact with someone asymptomatic saturates at trans_asymp·b and one with someone symptomatic at b. A per-source-state factor would be exact; with trans_asymp = 0.94 the difference is negligible.

**Acceptance test:** exp 18's `run.py` (`q1`, `q2`), against a starsim with the change. Agents ÷ network ODE should come out about 1 for every transition by contact, and set B's equilibrium within noise of 1. We can run it on request.

## For our repos, until then

- **Our ODE calibration uses homogeneous mixing,** where agents match the ODE (exps 13–17), so it doesn't depend on this.
- **Our matched-`beta` convention for `RandomNet`** ([proposal 01](../01_network_beta/proposal.md), CONVENTIONS.md 1) is the same model as Starsim's network ODE, and has the same gap. Proposal 01's "not submitted upstream" decision is superseded by `5ba36ebb`.
