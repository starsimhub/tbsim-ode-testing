# 25: Does Starsim's network fix make agents on a random network match the network ODE, and what does it change for our conversion?

**Status:** done; see [SUMMARY.md](SUMMARY.md). Design confirmed by the user (2026-11-10). Question 4's weekly agent arms were not run (user, 2026-11-10): their predicted errors are 0.5% or less. Their predictions are in `conversions.csv`.
**Track:** v4, in this experiment's own pinned environment (`experiments/25_starsim_network_fix/.venv`), so the v4 track's pins (`v4/`) and exps 12–24 stay reproducible.

- **starsim** `origin/v4` @ `d83f4f87` (4.0.0rc4). It includes `6dcad9fd`, Cliff's response to [proposal 03](../../proposals/03_network_reinfection/proposal.md).
- **tbsim** `origin/v4` @ `619c296` ("ss.fitting update").
- **Checked before designing:** tbsim runs at these pins as agents and as an ODE, on `RandomNet` and without it. `v4/tests/test_anchor.py`, `test_fastode.py` and `test_latent.py` pass, so the homogeneous v4 ODE still reproduces `TB_ODE` and `TBODE`. `test_treat.py` fails, because `ss.State` changed its arguments and `tbode4.TBTreat` passes `label` twice; treatment isn't used here.

Run from the repo root with `experiments/25_starsim_network_fix/.venv/bin/python`, with `PYTHONPATH=v4`.

## Background

On `ss.RandomNet(n_contacts=10, dur=0)`, with fresh contacts every step, a contact transmits at most once in a step of length D. So its chance saturates: a contact at hazard h transmits with probability 1 − e^(−h·D), less than h·D.

**What `6dcad9fd` changes:**

- **Agents:** every factor of a contact's hazard is now inside the exponent. An edge from u to v transmits with chance 1 − e^(−b·D·rr·κ_u·`rel_trans`[u]·`rel_sus`[v]·modifiers), where b is the per-contact `beta`, rr the transition's relative risk (1 for susceptibles, `rr_reinfection_*` for reinfection) and κ the source state's infectiousness (`trans_asymp` for asymptomatic TB, 1 for symptomatic). At the v4 track's pins, only rr·b saturated, and κ and the other factors scaled the chance linearly ([exp 22](../22_network_conversion/SUMMARY.md), finding 2).
- **Starsim's ODE of a network:** each transition by contact has its own saturation factor, (1 − e^(−x))/x with x = rr·κ·b·D, for each source state. Before, the first transition's factor (the susceptibles', rr = 1) was applied to every transition, so the ODE under-transmitted reinfection at a reduced risk. In set B at monthly steps, agents ran at 1.20× the network ODE's prevalence ([exp 18](../18_network_ode/SUMMARY.md)).

**What it means for converting a homogeneous calibration to the network.** Exp 22's closed-form conversion (b, and each rr′ = −ln(1 − rr·q)/(b·D) with q = 1 − e^(−b·D); κ unchanged) was exact while κ scaled the chance linearly. With κ in the exponent:

- **κ unchanged:** an asymptomatic contact transmits with 1 − e^(−κ·b·D), more than the κ·q the homogeneous model implies. Expected: prevalence too high, by a few percent at monthly steps.
- **κ converted too,** κ′ = −ln(1 − κ·q)/(b·D): susceptibles and reinfection from symptomatic sources match exactly, and so do susceptibles from asymptomatic sources. Reinfection from asymptomatic sources gets rr′·κ′, which is less than what's needed, −ln(1 − rr·κ·q)/(b·D). Expected: prevalence slightly too low, about 1% at monthly steps.

No choice of b, rr′ and κ′ matches all four source–target pairs exactly, because every factor multiplies inside the exponent.

## Questions

1. **Per transition:** for each transition by contact, do agents' flows on `RandomNet(10, dur=0)` now match Starsim's network ODE (exp 18 found 1.07–1.08 for reinfection at monthly steps)?
2. **Equilibrium:** do agents' prevalence, latent infection, cleared share and TB deaths match the network ODE in sets A and B (exp 18: set B 1.20 monthly, 1.05 weekly)?
3. **Per step:** do agents now transmit with κ inside the exponent, as the changelog says? (Exp 22's per-step diagnostic found κ linear at the old pins.)
4. **Conversion:** how far from the calibrated homogeneous ODE are agents on the network with κ unchanged (exp 22's rule) and with κ converted? Do they follow the predictions?
5. **Homogeneous mixing:** do agents without a network still match the homogeneous ODE, as in exp 21?

## Design

**Network:** `ss.RandomNet(n_contacts=ss.constant(10), dur=0)`. **Steps:** monthly and weekly. **Seeds:** new, from 25000 (v4's draws are keyed by seed).

**Questions 1–2 (exp 18's `q3`, `q1`, `q2`, network workflow only):**

- **Parameter sets:** exp 18's. Set A is tbsim's defaults (near threshold; reinfection of cleared people at rr = 1). Set B is exp 17's two-state model at φ = 0.01, every reinfection risk at π = 0.24, latent reinfection on.
- **Calibration first (`q3`):** the per-contact b that gives 500/100k active prevalence on Starsim's network ODE, at each step. Reported against the matched b from the mass-action `beta`, −ln(1 − `beta`·D/K)/D, which it equalled at exp 18's pins. The fixed network ODE transmits more per unit of b than before, in two ways: reinfection at a reduced risk gets its own factor, nearer 1 than the susceptibles', and asymptomatic sources (κ < 1) do too. So the calibrated b comes out below the matched b, by more in set B, where most infections are reinfections. This difference follows from the fix and isn't part of the test.
- **Per transition (`q1`):** 1,000,000 agents × 4 seeds, 24 steps from the network ODE's equilibrium. Each transition by contact's flow per person ÷ the network ODE's from the same start.
- **Equilibrium (`q2`):** 50,000 agents, 40 reps per arm, 1900–2000, started at the network ODE's equilibrium; agents ÷ network ODE averaged over 1950–2000.
- **Predicted:** about 1 everywhere (within about 2 SE), against exp 18's 1.07–1.08 per reinfection transition and 1.20 / 1.05 for set B's equilibrium.

**Question 3 (exp 22's `diagnose.py`):** sets 3 and 7 of exp 21's posterior draws, 2 reps each, monthly, 1900–2000, exp 22's κ′ equilibrium-matched arm (so the two rules differ visibly). Each step, each target class's infections ÷ the chance worked out from its contacts at the start of the step, with κ inside the exponent (1 − e^(−κ·rr·b·D) per asymptomatic contact) or linear (κ·(1 − e^(−rr·b·D))). Predicted: the in-exponent expectation now gives about 1.00, and the linear one about 1.02 (the mirror of exp 22's 0.975 against 0.995).

**Question 4 (exp 22's conversion):** exp 21's 12 draws from exp 20's posterior (homogeneous calibration), 50,000 agents on the network, 20 reps per set and arm, 1900–2000, started at each set's homogeneous equilibrium. Compared with the calibrated homogeneous ODE, as in exp 22.

| Arm | Per-contact `beta` | Reinfection risk | `trans_asymp` | Step |
|---|---|---|---|---|
| κ unchanged (`convk_month`, exp 22's rule) | b | rr′ | κ | Monthly |
| κ converted (`conv_month`) | b | rr′ | κ′ | Monthly |
| κ unchanged, weekly (`convk_week`; prediction only) | b, rr′ for weekly steps | | κ | Weekly |
| κ converted, weekly (`conv_week`; prediction only) | b, rr′, κ′ for weekly steps | | | Weekly |

**Predictions** come from the standalone ODE's network option (`tbode4.fastode`, `net_contacts`, `net_step`), with a new switch, `net_source_exp=1`, that puts κ inside the exponent as agents now do: each transition by contact runs at (K/D)·[(1 − e^(−rr·κ·b·D))·A + (1 − e^(−rr·b·D))·Y], where A and Y are the asymptomatic and symptomatic shares. They are worked out before the agents run (`conversions.csv`). As a cross-check, Starsim's own network ODE is run at each converted set and compared with them.

**Question 5 (homogeneous mixing):** exp 21's monthly arm (12 sets × 20 reps, 50,000 agents, no network) at the new pins, against the homogeneous ODE. Predicted: as exp 21 (prevalence 0.997 ± 0.007), since the change is about network edges.

**Measured, as in exps 18, 21 and 22:** agents ÷ ODE, mean ± SE over reps (and between sets for question 4), with a 10–90% range: active prevalence, the share asymptomatic, first-entry incidence per prevalent case (exp 19's `Incidence` analyzer), TB mortality, latent infection (recent and remote) and the cleared share.

## What would answer the questions

- **Proposal 03 is resolved** if questions 1–2 give about 1 for every transition and both sets, and question 3 confirms κ inside the exponent. Then the network convention can become "calibrate and run the ODE on the network directly" (CONVENTIONS.md 1, CLAUDE.md's network decisions), and we report back to Cliff.
- **For the conversion:** if the κ-converted arms are within about 1% of the homogeneous ODE at monthly steps, the conversion remains usable as an approximation (with κ converted). If not, a homogeneous calibration needs recalibrating on the network ODE before running network agents.
- **If question 5 moves,** the change affects homogeneous agents too, which would need its own experiment before any pin bump.

## Not covered

- Networks whose edges last several steps (households).
- Moving the v4 track's pins, and fixing `tbode4.TBTreat` for the new `ss.State`: decided after this experiment.

## Files

- `pyproject.toml`, `uv.lock`: this experiment's environment.
- `run.py`: `q3`, `q1`, `q2` (questions 1–2), `diagnose` (question 3), `conversions` and `convert` (question 4), `homog` (question 5), `summary`.
- `v4/tbode4/fastode.py`: the `net_source_exp` switch (off by default, so exp 22's predictions are unchanged).
- `SUMMARY.md`: findings, written after running.
