# Proposal 03: a saturation factor for each transition in Starsim's ODE of a network

**For:** Starsim (Cliff). A follow-up to `5ba36ebb`, "add beta correction for ODEs on networks".
**Status:** resolved upstream in starsim `6dcad9fd` (2026-11-09), and confirmed by [exp 25](../../experiments/25_starsim_network_fix/SUMMARY.md) at starsim `d83f4f87` with tbsim `619c296`.

- **In Starsim's ODE of a network:** each transition by contact, and each source state, has its own saturation factor, as suggested below.
- **For agents, done differently from below:** rather than leaving source infectiousness linear, every factor of a contact's hazard is inside the exponent (κ, `rel_trans`, `rel_sus`, modifiers), and the ODE matches that. This meets the acceptance criterion that agents and the ODE agree.
- **Exp 25's acceptance test:** set B's monthly prevalence is 1.016 ± 0.011 of the network ODE, against 1.204 before. Every reinfection transition runs with first infection.

The rest of this document is kept as the record of the request.
**Evidence:** [exp 18](../../experiments/18_network_ode/SUMMARY.md), starsim `bc2704b2` with tbsim `efeca56`; [exp 22](../../experiments/22_network_conversion/README.md)'s per-step diagnostic (`diagnose.py`), starsim `80c81f53` with tbsim `dcfb522`. Agents' network transmission is the same code at both commits.

## Summary

- **What the fix gets right:** for transmission at the disease's own per-contact rate, the ODE of a network now matches agents. Calibrated on the network ODE, the per-contact `beta` comes out exactly where agents need it.
- **What it misses:** the saturation factor is worked out once per exposure, from the first transition by contact that reads it. It's then applied to every transition on those contacts. A transition at a different per-contact rate, such as reinfection at a reduced risk rr·b, saturates differently as agents. The ODE gets it wrong by the ratio of the two factors.
- **The size of it:** for `tbsim.TB` with tb_LAI_TPT's calibrated natural history, on `ss.RandomNet(n_contacts=10, dur=0)` at monthly steps, agents' prevalence is 1.20× the network ODE's (1.05× at weekly steps).
- **Source infectiousness doesn't saturate as agents:** an asymptomatic contact transmits with κ times a symptomatic contact's chance, where κ is `trans_asymp`. The ODE's factor should leave it out.
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

### Which factors saturate

As agents (`transmit_from` and `transmit_together` in `solvers/agents_branching.py`), an edge from u to v transmits in a step with chance

  rel_trans[u] · s[v] · (1 − e^(−r·D)),

where:

- **r, which saturates,** is the transition's rate coefficient for its most infectious source state: β times the parameters written into the rate (here `rr_reinfection_*`).
- **rel_trans[u], linear,** is the source's coefficient relative to that largest one (κ for an asymptomatic source when κ < 1), times its per-agent `rel_trans`.
- **s[v], linear,** is the target's susceptibility: its `rel_sus`, its state's `rel_sus`, and the modifiers on the transmission (`susceptibility()` in `solvers/agents_stepping.py`).

So whether a factor saturates depends on where the model writes it, not on whether it is about transmissibility or susceptibility. A reduced reinfection risk written as `rr` in the rate saturates; the same risk written as `rel_sus` doesn't.

For tbsim, the chance per contact per step is then:

| Source ↓ / target → | Susceptible | Reinfected at rr |
|---|---|---|
| Symptomatic | 1 − e^(−b·D) | 1 − e^(−rr·b·D) |
| Asymptomatic | κ·(1 − e^(−b·D)) | κ·(1 − e^(−rr·b·D)) |

## Evidence

**Per transition (exp 18):** 1,000,000 agents × 4 seeds, 24 steps from equilibrium. Agents ÷ network ODE, flow per person, mean ± SE; the prediction is (1 − e^(−rr·b·D))/(rr·b·D) ÷ (1 − e^(−b·D))/(b·D).

| tb_LAI_TPT's natural history, rr = 0.24 | Monthly | Predicted | Weekly | Predicted |
|---|---|---|---|---|
| `transmission` (rr = 1) | 0.994 ± 0.005 | 1.000 | 0.995 ± 0.007 | 1.000 |
| `reinfect_cle` | 1.078 ± 0.003 | 1.075 | 1.011 ± 0.007 | 1.016 |
| `reinfect_rec` | 1.076 ± 0.008 | 1.075 | 1.021 ± 0.008 | 1.016 |
| `reinfect_remote` | 1.071 ± 0.010 | 1.075 | 0.989 ± 0.016 | 1.016 |

**At equilibrium (exp 18):** 50,000 agents, 40 reps, 1950–2000, prevalence agents ÷ network ODE: 1.204 ± 0.010 monthly and 1.052 ± 0.009 weekly. With tbsim's default rates, where reinfection of cleared people runs at rr = 1, the results are 1.008 ± 0.018 and 1.017 ± 0.018.

**The cause accounts for all of it:** the mass-action ODE with each `rr_reinfection_*` scaled by the ratio of the agents' factor to the ODE's predicts 1.201 (monthly) and 1.037 (weekly). With exp 18's κ = 0.94, how κ enters makes no visible difference.

**Source infectiousness (exp 22):** 50,000 agents on `RandomNet(10, dur=0)` at monthly steps, two parameter sets with κ = 0.63 and 0.73, 2 reps each, 1900–2000. Each step, each target's chance was worked out from its contacts at the start of the step. Infections observed ÷ expected, with SE about 0.003:

| Expected with asymptomatic contacts at | First infections | Reinfection of cleared people |
|---|---|---|
| 1 − e^(−κ·rr·b·D) (saturating) | 0.975, 0.978 | 0.993, 0.997 |
| κ·(1 − e^(−rr·b·D)) (linear) | **0.995, 0.998** | **1.002, 1.003** |

## Suggested change

Work out the factor for each transition by contact, from that transition's own per-contact rate, as agents do:

- **The rate:** β times the transition's rate parameters (here `rr_reinfection_*`) and the coefficient of its most infectious source state. That gives x_f = r_f·D and the factor (1 − e^(−x_f))/x_f.
- **Left out of the factor, as linear:** the other source states' coefficients relative to the largest, and `rel_trans`, `rel_sus` and the modifiers on the transmission. `saturation()` currently includes the mean infectiousness of the source states (`exposure_trans`) and `flow_mult`. The modifiers point is from reading the code and hasn't been tested.
- **Where to apply it:** contacts are stored per exposure and term, so either apply each transition's factor relative to the exposure's as a multiplier on that transition's rate, or store contacts per transition.

**An alternative on the agents' side:** apply rate parameters such as `rr` linearly, as `rel_sus` already is. Only β would then saturate, one factor per exposure, and `5ba36ebb`'s correction would be exact on its own. Either works for us, as long as the ODE and agents agree.

**Acceptance test:** exp 18's `run.py` (`q1`, `q2`) and exp 22's `diagnose.py`, against a starsim with the change. Agents ÷ network ODE should come out about 1 for every transition by contact, and set B's equilibrium within noise of 1. We can run it on request.

## For our repos, until then

- **Our ODE calibration uses homogeneous mixing,** where agents match the ODE (exps 13–17, 21), so it doesn't depend on this.
- **Our matched-`beta` convention for `RandomNet`** ([proposal 01](../01_network_beta/proposal.md), CONVENTIONS.md 1) is the same model as Starsim's network ODE, and has the same gap. Proposal 01's "not submitted upstream" decision is superseded by `5ba36ebb`.
- **A homogeneous calibration should convert to `RandomNet(dur=0)` in closed form** ([exp 22](../../experiments/22_network_conversion/README.md), being run): the per-contact `beta` b and each `rr_reinfection_*` are converted to match per contact and per step, and `trans_asymp` is left unchanged, since agents apply it linearly.
