# 25: Summary

**Track:** v4, in this experiment's own environment: starsim `origin/v4` @ `d83f4f87` (4.0.0rc4, which includes `6dcad9fd`) and tbsim `origin/v4` @ `619c296` (`versions.json`). The homogeneous v4 ODE still reproduces `TB_ODE` and `TBODE` (`v4/tests/test_anchor.py`, `test_fastode.py` and `test_latent.py` pass). `tbode4.TBTreat` doesn't build at these pins (`test_treat.py`), because `ss.State` changed its arguments.

**What `6dcad9fd` changed:**

- **Agents on a network:** every factor of a contact's hazard is inside the exponent. An edge transmits with chance 1 − e^(−rr·κ·b·D·`rel_trans`·`rel_sus`·modifiers). Here b is the per-contact `beta`, D the step, rr the transition's relative risk (1 for susceptibles, `rr_reinfection_*` for reinfection) and κ the source state's infectiousness (`trans_asymp` for asymptomatic TB, 1 for symptomatic).
- **Starsim's ODE of a network:** each transition by contact, and each source state, has its own saturation factor (1 − e^(−x))/x, with x = rr·κ·b·D. This was [proposal 03](../../proposals/03_network_reinfection/proposal.md)'s request.

**Setting:** `ss.RandomNet(n_contacts=10, dur=0)`, fresh contacts every step.

- **Questions 1–2:** [exp 18](../18_network_ode/SUMMARY.md)'s sets. Set A is tbsim's defaults. Set B is exp 17's two-state model at φ = 0.01, with every reinfection risk at π = 0.24.
- **Questions 3–5:** exp 21's 12 draws from exp 20's homogeneous calibration.

## Finding 1: agents on the network now match Starsim's network ODE, reinfection included

**Per transition:** 1,000,000 agents × 4 seeds, 24 steps from the network ODE's equilibrium. The values are agents ÷ network ODE, flow per person, mean ± SE.

| Transition | Set B, monthly | Exp 18 | Set B, weekly | Exp 18 |
|---|---|---|---|---|
| First infection (rr = 1) | 0.982 ± 0.004 | 0.994 ± 0.005 | 0.998 ± 0.005 | 0.995 ± 0.007 |
| Reinfection of cleared (rr = π) | 0.986 ± 0.004 | **1.078 ± 0.003** | 0.993 ± 0.002 | 1.011 ± 0.007 |
| Reinfection of recovered (rr = π) | 0.973 ± 0.008 | **1.076 ± 0.008** | 0.996 ± 0.035 | 1.021 ± 0.008 |
| Remote → recent (rr = π) | 0.969 ± 0.003 | **1.071 ± 0.010** | 0.982 ± 0.039 | 0.989 ± 0.016 |

- **Set A:** 0.993–1.009 for every transition at both steps (`q1_summary.csv`).
- **The reinfection excess is gone:** reinfection now runs at the same ratio as first infection.
- **A small shortfall remains in set B at monthly steps:** 1.5–3%, in first infections as well, so it isn't specific to reinfection. Exp 18 had about 1% in set A (0.991–0.992). It doesn't show at equilibrium (below), and isn't explained here.

**At equilibrium:** 50,000 agents, 40 reps, 1900–2000, started at the network ODE's equilibrium. The values are prevalence agents ÷ network ODE over 1950–2000, mean ± SE (`q2_summary.csv`).

| Set | Monthly | Exp 18 | Weekly | Exp 18 |
|---|---|---|---|---|
| A | 1.005 ± 0.019 | 1.008 ± 0.018 | 0.983 ± 0.017 | 1.017 ± 0.018 |
| B | **1.016 ± 0.011** | **1.204 ± 0.010** | 0.999 ± 0.008 | 1.052 ± 0.009 |

- **TB deaths and the cleared share move with prevalence,** all within 1.5 SE of 1.
- **Latent infection reads high at monthly steps** (A 1.09, B 1.04), as with homogeneous mixing (exp 21: 1.04). v4 dates infections at the end of their step, so end-of-step snapshots of recent latent infection read high.

**What the network ODE's calibration now gives:** the per-contact b that gives 500/100k on Starsim's network ODE (`run_q3.log`), against the matched value from the mass-action `beta`, −ln(1 − `beta`·D/K)/D. At exp 18's pins the two were equal.

| Set | Monthly | Weekly |
|---|---|---|
| A | 0.994 × matched | 0.999 × matched |
| B | **0.939 × matched** | 0.987 × matched |

The fixed ODE transmits more per unit of b. Reinfection at a reduced risk saturates less than first infection does, and so do asymptomatic sources (κ < 1). Set B, where most infections are reinfections, needs 6% less b at monthly steps.

## Finding 2: agents put κ inside the exponent

This is exp 22's per-step diagnostic: sets 3 and 7, 2 reps each, monthly, 1900–2000, on exp 22's κ′ parameters. Each step, each target class's infections are compared with the chance worked out from its contacts at the start of the step (`diagnose_summary.csv`). SE is about 0.003 for the first two rows and about 0.01 for the others.

| Infections observed ÷ expected, set 3, set 7 | κ inside the exponent | κ linear (exp 22's finding) |
|---|---|---|
| First infection | **0.996, 1.003** | 1.017, 1.023 |
| Reinfection of cleared | **0.993, 0.999** | 1.002, 1.005 |
| Reinfection of recovered | 1.020, 1.006 | 1.029, 1.012 |
| Remote → recent | 0.964, 0.948 | 0.972, 0.954 |

- **The rule has reversed:** exp 22 found 0.975 and 0.978 for first infections with κ inside the exponent, and 0.995 and 0.998 with κ linear.
- **Remote → recent reads 4–5% low under both rules.** It did in exp 22 too (0.93–0.94), so it isn't new. It probably comes from how the diagnostic counts this class: people in remote latent infection also leave it by other routes within the step. Not pursued.
- **The network is as intended:** mean degree 10.000 (variance 5).

## Finding 3: Starsim's network ODE matches the standalone network ODE

`tbode4.fastode` has a new network option, `net_source_exp=1`, that puts κ inside the exponent as agents now do. For each of the 12 sets × 4 converted parameter sets, Starsim's network ODE gives 0.997–1.004 of its prevalence (`conversions.csv`). So the standalone ODE is a faithful copy of the fixed network ODE, and serves as the prediction for finding 4.

## Finding 4: converting a homogeneous calibration to the network is now approximate

These are exp 22's conversions, now on 50,000 agents, 20 reps per set, monthly, 1900–2000. With K = 10 and q = 1 − e^(−b·D):

- the per-contact b is −ln(1 − `beta`·D/K)/D;
- each `rr_reinfection_*` becomes −ln(1 − rr·q)/(b·D);
- κ is either left unchanged (exp 22's rule) or converted in the same way, κ′ = −ln(1 − κ·q)/(b·D).

| Agents on the network ÷ homogeneous ODE, 1950–2000 (mean ± SE between 12 sets; range of set means) | κ unchanged | Predicted | κ converted | Predicted |
|---|---|---|---|---|
| Active prevalence | 1.018 ± 0.009 (0.97–1.09) | 1.023 (1.01–1.05) | **0.978 ± 0.006** (0.94–1.02) | 0.981 (0.97–0.99) |
| TB mortality | 1.017 ± 0.009 | 1.023 | 0.979 ± 0.006 | 0.981 |
| Share asymptomatic | 0.999 ± 0.001 | 1.000 | 0.999 ± 0.001 | 1.000 |
| Incidence per prevalent case | 0.999 ± 0.002 | 1.000 | 0.998 ± 0.001 | 1.000 |
| Remote latent | 1.016 ± 0.008 | 1.022 | 0.975 ± 0.005 | 0.981 |
| Cleared | 1.004 ± 0.002 | 1.007 | 0.993 ± 0.001 | 0.997 |

Sources: `conv_pooled.csv`, and by set, `conv_by_set.csv`.

- **Both arms follow their predictions** (largest |z| against the prediction for prevalence: 1.4 with κ converted, 4.6 for one set with κ unchanged).
- **Neither matches the homogeneous ODE exactly.**
  - With κ unchanged, asymptomatic contacts transmit too much, so prevalence runs about 2% high on average and up to 5% (predicted).
  - With κ converted, reinfection from asymptomatic contacts runs too low, because κ′·rr′ falls short of the converted product of κ and rr. So prevalence runs about 2% low, up to 3.4%.
  - No choice of b, rr′ and κ′ is exact, since every factor multiplies inside the exponent.
- **Ratios per case are unaffected:** the share asymptomatic and incidence per prevalent case are within 0.2%. The error is in the level of the epidemic, not in its natural history.
- **At weekly steps (predictions only; the agents weren't run):** prevalence 1.005 with κ unchanged and 0.996 with κ converted.

## Finding 5: homogeneous agents are unchanged

These are exp 21's monthly arm at the new pins (12 sets × 20 reps, no network): agents ÷ homogeneous ODE, mean ± SE between sets (`homog_pooled.csv`).

| Output | Exp 25 | Exp 21 |
|---|---|---|
| Active prevalence | 1.009 ± 0.005 | 0.997 ± 0.007 |
| Latent infection | 1.051 ± 0.006 | 1.041 ± 0.007 |
| Recent latent (end-of-step snapshot) | 1.43 ± 0.04 | 1.42 ± 0.04 |

The share asymptomatic, incidence per prevalent case, mortality and the cleared share are all within 0.2–0.6% of 1.

## What this means

- **Proposal 03 is resolved upstream.** Agents and Starsim's network ODE agree for every transition by contact, and set B's equilibrium gap (1.20 at monthly steps) is gone. Cliff fixed it on both sides: each transition's own factor in the ODE, and every factor inside the exponent for agents. That differs from the proposal, which had κ linear, but meets its criterion that the two agree.
- **On a network, calibrate the ODE on the network.** Starsim's network ODE is now the same model as agents on `RandomNet(dur=0)`, so a calibration on it carries over to agents with no conversion. Its per-contact b differs from the matched b by up to 6% (set B, monthly), so `beta_edge_from_ode` alone is no longer the right carry-over.
- **Our closed-form conversion (CONVENTIONS.md 1) is now approximate:** about ±2% of prevalence at monthly steps, under 0.5% at weekly steps, whether or not κ is converted. Exp 22's rule (κ unchanged) overstates prevalence, so `tbode4.rr_edge_from_ode`'s note that κ needs no conversion holds only at the old pins.
- **Homogeneous mixing, the default, is unaffected.**
- **Before moving the v4 pins:** fix `TBTreat` for the new `ss.State`, and rerun the treatment checks (exp 23).

## Not covered

- Networks whose edges last several steps (households).
- The weekly conversion arms as agents (user, 2026-11-10): predicted within 0.5%.
- Set B's 1.5–3% monthly shortfall per transition over 24 steps, and the diagnostic's remote → recent reading.
