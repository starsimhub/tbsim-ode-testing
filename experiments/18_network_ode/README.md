# 18: Does Starsim v4's ODE of a network make our per-contact `beta` convention unnecessary?

**Track:** v4, in this experiment's own pinned environment (`experiments/18_network_ode/.venv`), so the v4 track's pins (`v4/`) and exps 12–17 stay reproducible.

- **starsim** `origin/v4` @ `bc2704b2`: the merge that adds `5ba36ebb`, "add beta correction for ODEs on networks". Later commits (`308a27a5` on) rename transitions that tbsim `efeca56` uses, so it can't go newer yet.
- **tbsim** `origin/v4` @ `efeca56` ("updates for v4-rc3").
- **Also new since the v4 track's pins** (starsim `80c81f53`): about 40 other commits, including a change to how results are indexed by step and a chain-binomial fix. `v4/tests/test_anchor.py` passes in this environment, so the mass-action v4 ODE still reproduces `TB_ODE` and `TBODE`.

Run from the repo root with `experiments/18_network_ode/.venv/bin/python`, with `PYTHONPATH=v4` for `tbode4` and Cliff's helpers.

## Background

- **The gap:** on `ss.RandomNet(n_contacts=K, dur=0)`, with fresh contacts every step, each contact transmits at most once in a step, with probability 1 − e^(−b·D) for a per-contact rate b over a step D. So K contacts give a force of infection of K·(1 − e^(−b·D))/D per infectious share, not the K·b of mass action. At monthly steps with `beta` ≈ 10 and K = 10, that's about 4% less transmission. Near threshold it's about 25% less prevalence (exps 12–13).
- **Our convention** ([CONVENTIONS.md](../../CONVENTIONS.md) 1, user 2026-11-05): calibrate the mass-action ODE, and give agents on the network the matched per-contact `beta`, −ln(1 − `beta`·dt/K)/dt (`tbode4.beta_edge_from_ode`). The mass-action ODE is then the model.
- **The Starsim fix** (`5ba36ebb`): Starsim's ODE of a network now multiplies each contact's transmission by (1 − e^(−b·D))/(b·D), what agents achieve per step. D is the longer of the network's step and the disease's. So the same model, TB on `RandomNet` with per-contact b = `beta`/K, should give the same answer as an ODE and as agents. Cliff's suggestion is to calibrate on the network ODE and drop the mapping, keeping it only to carry over a mass-action `beta`.
- **A possible gap in the fix:** the factor is worked out once per exposure, from the first transition by contact that reads it. In `tbsim.TB` that is `transmission` (susceptible → latent), at rate b. The reinfection transitions (cleared, recovered, and remote latent → recent) run at π·b on the same contacts. As agents, each of their contacts transmits with probability 1 − e^(−π·b·D), which is closer to π·b·D. If the ODE applies the susceptibles' factor to them, it under-transmits reinfection. In set B at monthly steps (b·D ≈ 0.18 with b = `beta`/K), the factor is 0.91 for susceptibles against about 0.98 for reinfection at π = 0.24, and about 60% of set B's infections are reinfections. Our convention has the mirror image: the matched b makes agents' susceptible infections equal the ODE's, but agents' reinfections then run about (1 − π)·b·D/2 above it.

## Questions

1. **One step:** on `RandomNet(n_contacts=10)`, do the network ODE's flows match agents' for each transition by contact (transmission and each reinfection)? Against our convention's flows: agents with the matched b, against the mass-action ODE.
2. **Equilibrium:** do agents on the network match the network ODE (new workflow), and how does that compare with our convention (mass-action ODE, agents with the matched b)?
3. **Calibration:** how far apart are `beta` calibrated on the network ODE and on the mass-action ODE, and does it matter which is calibrated?

## Design

**Network:** `ss.RandomNet(n_contacts=10, dur=0)`, as in tb_LAI_TPT and exps 12–13. **Steps:** monthly and weekly.

**Parameter sets:**

- **Set A:** tbsim's default rates (near threshold). Its reinfection risks are tbsim's defaults: cleared 1.0, recovered 0.21, so most reinfection is at the susceptibles' rate.
- **Set B, two-state:** exp 17's model, set B with φ = 0.01 and every reinfection risk at π = 0.24, latent reinfection on (remote → recent). This is the current model and the stronger test of reinfection.

**Two workflows are compared:**

| Workflow | ODE | Agents' per-contact `beta` | `beta` calibrated on |
|---|---|---|---|
| Network ODE (new) | ODE on the same `RandomNet` | b, calibrated | the network ODE, at this step |
| Our convention | Mass-action ODE | matched b, −ln(1 − `beta`·dt/K)/dt | the mass-action ODE |

**Question 3 (calibration), run first:** `beta` for 500/100k active prevalence on the mass-action ODE, and the per-contact b on the network ODE at monthly and weekly steps (by root-finding near the matched b), for both sets. If the calibrated b equals the matched b, the two workflows give the same agents.

**Question 1 (per transition):** 1,000,000 agents × 4 seeds, started at each set's equilibrium, 24 steps. For each transition by contact, the agents' flow per person, summed over the steps, ÷ the ODE's from the same start, in both workflows. The prediction is the agents' per-contact factor for the transition's risk ratio rr, (1 − e^(−rr·b·D))/(rr·b·D), over the factor the ODE applies.

**Question 2 (equilibrium):** 50,000 agents, started at the ODE's equilibrium, 1900–2000, 40 reps per arm, for each set at monthly and weekly steps (4 arms, 160 runs). Measured: agents ÷ ODE over 1950–2000 (mean ± SE, 10–90% range), as in exps 13 and 17, for active prevalence, latent infection, the cleared share and TB deaths. Only the network workflow is run, because question 3 showed the convention gives the same agents and the same ODE.

**Explaining the gap (`explain`):** the mass-action ODE at the calibrated `beta`, with each reinfection risk scaled by the agents' factor for it over the ODE's. It predicts the agents' equilibrium if the per-transition factor is the whole gap.

**Seeds:** new (from 18000), since v4's draws are keyed by seed.

## What would answer the questions

- **If the network ODE matches agents in both sets** (within about 2 SE, as exps 13 and 17), the convention can change to "calibrate and run the model on the network as an ODE; use `beta_edge_from_ode` only to carry over a mass-action `beta`" (Cliff's note).
- **If it matches in set A but not set B,** and question 1 traces the gap to reinfection, the fix needs each transition's own factor; that would go to Cliff, and the convention stays until it's fixed.
- **Our convention's own reinfection mismatch,** if it shows in question 1, applies today, and its size at equilibrium matters for exps 12–13's network results (set A only).

## Not covered

- Networks whose edges last several steps (`dur` > 0), such as tb_LAI_TPT's household network: the fix takes them a step at a time, which Cliff notes undercorrects.
- Moving the v4 track's pins: decided after this experiment.

## Files

- `pyproject.toml`, `uv.lock`: this experiment's environment.
- `run.py`: the questions (`q3`, `q1`, `q2`, `explain`, in that order) → `q3.obj`, `q1.obj`, `q1_summary.csv`, `q2.obj` (per-rep yearly means and the ODE runs), `q2_summary.csv`, `explain.json`, `run_*.log`, `explain.log`, `versions.json`.
- `SUMMARY.md`: findings, written after running.
