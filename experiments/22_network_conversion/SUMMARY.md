# 22: Summary

**Track:** v4 (tbsim `origin/v4` @ `dcfb522`, starsim `origin/v4` @ `80c81f53`; `versions.json`).

**Setting:**

- **Parameter sets:** exp 21's 12 draws from exp 20's calibrated posterior (homogeneous mixing; `beta` 10.9–29.0, π 0.20–0.66, `trans_asymp` κ 0.57–0.88).
- **Agents:** 50,000 on `ss.RandomNet(n_contacts=10, dur=0)`, started at each set's equilibrium, 1900–2000, 20 reps per set and arm.
- **The reference:** the calibrated homogeneous ODE. Incidence is the first entry to asymptomatic TB per episode, as in exp 21.

## Finding 1: the closed-form conversion reproduces the calibrated ODE on the network

Convert the per-contact `beta` and each reinfection risk, and leave `trans_asymp` unchanged. With K = 10 contacts, step D and q = 1 − e^(−b·D):

- b = −ln(1 − `beta`·D/K)/D;
- π′ = −ln(1 − π·q)/(b·D), for each `rr_reinfection_*`;
- κ unchanged.

Across the 12 sets, b·D is 0.10–0.28 monthly (0.02–0.06 weekly), and π′/π is 0.90–0.98 monthly (0.98–1.00 weekly).

| Agents on the network ÷ homogeneous ODE, 1950–2000, averaged over 12 sets (mean ± SE between sets) | Converted, monthly | Converted, weekly | `beta` only, monthly |
|---|---|---|---|
| Active prevalence | **1.000 ± 0.006** (0.98–1.04) | **0.999 ± 0.006** (0.98–1.03) | 1.160 ± 0.018 (1.07–1.26) |
| TB mortality | 0.998 ± 0.005 | 0.999 ± 0.006 | 1.157 ± 0.018 |
| Share asymptomatic | 0.999 ± 0.001 | 1.000 ± 0.001 | 0.999 ± 0.001 |
| Incidence (first entry) per prevalent case | 0.999 ± 0.001 | 1.001 ± 0.001 | 0.999 ± 0.001 |
| Prevalence ÷ incidence | 1.002 ± 0.001 | 0.999 ± 0.001 | 1.002 ± 0.001 |
| Mortality ÷ incidence | 1.000 ± 0.002 | 0.999 ± 0.003 | 0.999 ± 0.002 |
| Remote latent | 0.997 ± 0.005 | 0.997 ± 0.004 | 1.151 ± 0.017 |
| Cleared | 0.997 ± 0.001 | 0.998 ± 0.001 | 1.017 ± 0.001 |
| Recent latent (end-of-step snapshot) | 1.42 ± 0.04 | 1.09 ± 0.01 | 1.65 ± 0.05 |

Ranges in brackets are of the 12 set means. Sources: `pooled.csv`; by set, `by_set.csv` and `agents_vs_ode.png`.

- **Every set is within about 2 SE of 1** for the converted arms. The largest |z| for prevalence is 2.3 monthly and 1.5 weekly, over 12 sets × 2 steps.
- **That's as close as exp 21's homogeneous agents** (prevalence 0.997 ± 0.007 monthly, 0.987 ± 0.007 weekly).
- **Recent latent infection reads high in monthly snapshots,** by the same amount as with homogeneous mixing (exp 21: 1.42). v4 dates infections at the end of their step. The cleared share is 0.2–0.3% low, the small finite-population effect seen in exps 17 and 21.
- **Converting `beta` alone,** our convention (CONVENTIONS.md 1), overstates prevalence by 16% on average and by up to 26%. The excess comes from reinfection: on the network, a reduced risk saturates less than the susceptibles' risk, so reinfection runs above π times a susceptible's rate (exp 18). Ratios per prevalent case are unaffected, so the error is in the level of the epidemic, not in its natural history.

## Finding 2: v4's agents scale a contact's chance linearly by the source's infectiousness

How agents transmit over an edge in a step (`starsim/solvers/agents_branching.py`):

- **Saturates:** the transition's rate for its most infectious source state, r = rr·b, giving 1 − e^(−r·D).
- **Linear:** everything else multiplies that chance. That's the other source states' coefficients relative to the largest (κ for asymptomatic TB), per-agent `rel_trans` and `rel_sus`, and modifiers on the transmission.

Whether a factor saturates therefore depends on where the model writes it, not on whether it's about transmissibility or susceptibility.

The per-step diagnostic (`diagnose.py`; sets 3 and 7, 2 reps each, every step 1900–2000) measured this directly:

| Infections observed ÷ expected from each target's contacts at the start of the step (SE about 0.003) | First infections | Reinfection of cleared people |
|---|---|---|
| Asymptomatic contacts at 1 − e^(−κ·rr·b·D) (saturating) | 0.975, 0.978 | 0.993, 0.997 |
| Asymptomatic contacts at κ·(1 − e^(−rr·b·D)) (linear) | **0.995, 0.998** | **1.002, 1.003** |

The network itself is as intended: mean degree 10.000 (variance 5), 5 self-loops per step out of 250,000 edges.

**What converting κ costs:** three arms also converted κ, to κ′ = −ln(1 − κ·q)/(b·D) < κ, on the saturating assumption. Transmission from asymptomatic sources then runs κ′/κ too low. Each of these arms follows the network ODE's prediction under the linear rule:

| Prevalence, agents on the network ÷ homogeneous ODE (mean ± SE over 12 sets) | Agents | Predicted (linear κ) | Agents ÷ predicted |
|---|---|---|---|
| b, π′ and κ′, monthly | 0.951 ± 0.008 | 0.958 (0.92–0.98) | 0.992 ± 0.005 |
| b, κ′ and π′ matched at the equilibrium, monthly | 0.974 ± 0.010 | 0.976 (0.95–0.99) | 0.997 ± 0.008 |
| The same, weekly | 0.982 ± 0.006 | 0.995 (0.99–1.00) | 0.987 ± 0.006 |
| `beta` only, monthly | 1.160 ± 0.018 | 1.155 (1.08–1.24) | 1.005 ± 0.004 |

The predictions use the standalone ODE's network option (`tbode4.fastode`, `net_contacts`/`net_step`): each transition by contact at (K/D)·(1 − e^(−rr·b·D))·(κ·A + Y), where A and Y are the asymptomatic and symptomatic shares. The weekly κ′ arm's 0.987 is 2 SE from its prediction, and is likely noise, as exp 21's weekly figure was.

## What this means

- **A homogeneous calibration carries over to agents on `RandomNet(dur=0)` without recalibrating:** convert b and each `rr_reinfection_*` as above, and leave `trans_asymp` and every non-contact parameter unchanged. This needs no change to Starsim.
- **Converting `beta` alone isn't enough when reinfection runs at a reduced risk:** prevalence is 7–26% too high at monthly steps, depending on the set.
- **For Starsim's ODE of a network** ([proposal 03](../../proposals/03_network_reinfection/proposal.md)): the per-transition factor should use the transition's own rate for its most infectious source, and leave source infectiousness and the other per-agent factors linear, as the proposal sets out. When the fix lands, its ODE should reproduce the predictions above.
- **Not covered:** networks whose edges persist (households), which have no conversion.
