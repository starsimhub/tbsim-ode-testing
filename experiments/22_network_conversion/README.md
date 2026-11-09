# 22: Converting the calibrated homogeneous-mixing parameters for a random network

**Status:** done; see [SUMMARY.md](SUMMARY.md). Design confirmed by the user (2026-11-07). The converted arms, with κ unchanged, were added after a per-step diagnostic (`diagnose.py`) showed how v4's agents apply source infectiousness.
**Track:** v4 (`v4/.venv`; tbsim `origin/v4` @ `dcfb522`, starsim `origin/v4` @ `80c81f53`).

## Background

- **The calibration** (exps 20–21) uses homogeneous mixing, where v4's agents match its ODE. tb_LAI_TPT instead runs its agents on `ss.RandomNet(n_contacts=10, dur=0)`, with fresh contacts every step.
- **How agents transmit on such a network** (`starsim/solvers/agents_branching.py`; measured in `diagnose.py`): a contact transmits at most once in a step D.
  - **The target's relative risk saturates.** From a symptomatic contact, someone at relative risk rr (1 for susceptibles, π for reinfection) is infected with chance 1 − e^(−rr·b·D), where b is the per-contact `beta`. So reinfection does not scale with π as in mass action ([exp 18](../18_network_ode/SUMMARY.md)).
  - **The source's relative infectiousness doesn't.** An asymptomatic contact transmits with κ times a symptomatic contact's chance, where κ is `trans_asymp`.

  Starsim's ODE of a network applies the susceptibles' factor to every transition, so it under-transmits reinfection relative to agents ([proposal 03](../../proposals/03_network_reinfection/proposal.md) asks Cliff to fix this). Our convention (CONVENTIONS.md 1) converts only `beta`, so it has the same gap.
- **The conversion:** a homogeneous parameter set can be re-expressed for the network in closed form, with K = 10 and q = 1 − e^(−b·D):
  - b = −ln(1 − `beta`·D/K)/D (`tbode4.beta_edge_from_ode`), so that K·q = `beta`·D;
  - π′ = −ln(1 − π·q)/(b·D), for each `rr_reinfection_*`, so that a reinfection's chance per contact is π·q;
  - κ unchanged.

  Then agents on the network transmit, per contact and per step, exactly what the homogeneous model's flows give over a step for every combination of target and source: q, κ·q, π·q and κ·π·q. The other parameters carry over unchanged, because they don't involve contacts.
- **This needs no change to Starsim:** only Starsim's ODE of a network needs proposal 03's fix.

## Question

Do v4 agents on `RandomNet(10, dur=0)`, with the converted parameters, reproduce the calibrated homogeneous ODE? If they do, a calibration with homogeneous mixing can be reused for agents on a random network, without recalibrating. And how far off is converting `beta` alone, as our convention does?

## Design

**Parameter sets:** exp 21's 12 draws from exp 20's posterior (`../21_calibrated_comparison/sets.obj`), with the same equilibria.

**Arms:** agents, 50,000, on `ss.RandomNet(n_contacts=10, dur=0)`, started at each set's equilibrium, 1900–2000, 20 reps per set and arm (seeds 22000 + 1000·set + 100·arm + rep).

| Arm | Per-contact `beta` | Reinfection risk | `trans_asymp` | Step | Network standalone ODE's prediction, prevalence ÷ calibrated (range over sets) |
|---|---|---|---|---|---|
| **Converted** (`convk_month`) | b | π′ | κ | Monthly | 1 |
| **Converted, weekly** (`convk_week`) | b, π′ for weekly steps | | κ | Weekly | 1 |
| **`beta` only** (`beta_month`, our convention) | b | π | κ | Monthly | 1.08–1.24 |
| κ′, per combination (`conv_month`) | b | π′ | κ′ | Monthly | 0.92–0.98 |
| κ′, equilibrium-matched (`conveq_month`) | b | π′ matched at equilibrium | κ′ | Monthly | 0.95–0.99 |
| κ′, equilibrium-matched, weekly (`conveq_week`) | b, κ′, π′ for weekly steps | | | Weekly | 0.99–1.00 |

The three κ′ arms show what converting `trans_asymp` costs. They convert it as κ′ = −ln(1 − κ·q)/(b·D), which assumes the source's infectiousness saturates as a target's risk does. Since agents apply it linearly, transmission from asymptomatic sources then runs κ′/κ too low. The equilibrium-matched variant chooses π′ so that the total force of reinfection at the calibrated equilibrium matches under that same assumption.

**References:**

- **The calibrated homogeneous ODE** at each set: v4's ODE for state shares and mortality, and the standalone ODE for first-entry incidence, as in exp 21.
- **Exp 21's homogeneous agents,** reused, not rerun, as the benchmark for what "matching" looks like.
- **Predictions:** the standalone ODE with the network option (`net_contacts`, `net_step` in `tbode4.fastode`, off by default). Each transition by contact runs at (K/D)·(1 − e^(−rr·b·D))·(κ·A + Y), as agents do, where A and Y are the asymptomatic and symptomatic shares. It's a reference version of proposal 03's fix.

**Diagnostic** (`diagnose.py`): two sets (3 and 7) × 2 reps, κ′ equilibrium-matched arm, 1900–2000. Each step, an intervention records the network's edges and each agent's state before transmission. Each class's infections in the step are then compared with the chance worked out from its contacts, with asymptomatic contacts saturating (1 − e^(−κ·rr·b·D)) or linear (κ·(1 − e^(−rr·b·D))). It also records the network's realised degree.

**Measured:** agents ÷ homogeneous ODE, the mean over 1950–2000, for each set and arm (mean ± SE over reps), and averaged over the 12 sets. The outputs are as in exp 21:

- active prevalence and the share asymptomatic;
- first-entry incidence per prevalent case, and prevalence ÷ incidence;
- TB mortality and mortality ÷ incidence;
- latent infection (recent and remote), and the cleared share.

Plus each set's converted values (b, π′, κ′) against the originals.

**Figures:** agents ÷ ODE by set for each arm, with each arm's prediction. The raw outputs behind it are saved with `sc.save`.

## What would answer the question

- **The converted arms match the ODE** if every set is within about 2 SE of 1 and the average over sets is within about 1%, as exp 21's homogeneous agents did.
- **The other arms follow their predictions,** which tests the per-contact mechanism, including how κ enters.
- **Then:** a homogeneous calibration converts to network parameters in closed form. Once Starsim's network ODE has proposal 03's fix, it should reproduce these predictions (exp 18's acceptance test).

## Not covered

- Networks whose edges persist (households): no conversion exists.
- **Starsim's own network ODE, which still lacks the per-transition fix.** It is left out on purpose rather than patched locally:
  - the v4 track's pinned starsim predates even Starsim's first network correction (`5ba36ebb`);
  - a local patch of solver internals would break the rule of pinning exact upstream commits;
  - Cliff's implementation may differ, and will need its own test.

  The standalone ODE's network option stands in for the fixed ODE. When the fix lands, rerun exp 18's `q1`/`q2`, this experiment's converted arms and `diagnose.py` with Starsim's network ODE, in an environment pinned to that commit.

## Files

- `run.py`: `conversions` → `conversions.obj`, `conversions.csv`, `conversions.log` (each set's converted parameters and the predictions); with no argument it also runs the agents (resuming) → `runs.obj`, `by_set.csv`, `pooled.csv`, `agents_vs_ode.png`, `versions.json`; `summary` redraws. Uses exp 21's sets, ODE runs and outputs, and exp 19's `Incidence` analyzer. Logs: `run.log` (the κ′ and `beta`-only arms) and `run2.log` (the converted arms, and the summary of all six).
- `diagnose.py` → `diagnose.obj` (per step, by class: infections observed and expected); `diagnose_predictions.csv` (the κ′ and `beta`-only arms against predictions with κ linear, worked out before the network option was corrected).
- `v4/tbode4/fastode.py`: the network option (`net_contacts`, `net_step`; `force()`). Off by default; `v4/tests/test_fastode.py` passes.
- `SUMMARY.md`: findings, written after running.
