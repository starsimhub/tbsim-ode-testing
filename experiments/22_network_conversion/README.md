# 22: Converting the calibrated homogeneous-mixing parameters for a random network

**Status:** design confirmed by the user (2026-11-07), with the equilibrium-matched arm added after the predictions showed the cross term. Conversions and predictions done (`conversions.csv`); agent runs not yet started (`run.py`, resumes).
**Track:** v4 (`v4/.venv`; tbsim `origin/v4` @ `dcfb522`, starsim `origin/v4` @ `80c81f53`).

## Background

- **The calibration** (exps 20–21) uses homogeneous mixing, where v4's agents match its ODE. tb_LAI_TPT instead runs its agents on `ss.RandomNet(n_contacts=10, dur=0)`, with fresh contacts every step.
- **Why the network differs** ([exp 18](../18_network_ode/SUMMARY.md)): on such a network, a contact transmits at most once in a step, with probability 1 − e^(−x) for a per-contact rate x/D over a step D. A transition at relative risk rr (reinfection at π, say), from a source of relative infectiousness κ (`trans_asymp` for asymptomatic TB), therefore saturates at its own rate, rr·κ·b·D. Starsim's ODE of a network applies the susceptibles' factor to every transition, so it under-transmits reinfection relative to agents ([proposal 03](../../proposals/03_network_reinfection/proposal.md) asks Cliff to fix this). Our convention (CONVENTIONS.md 1) matches only transmission to susceptibles, so it has the same gap.
- **The conversion:** because each transition saturates at its own rate, a homogeneous parameter set can be re-expressed for the network with three parameters converted. The rest carry over unchanged, because they don't involve contacts. With K = 10 and step D:
  - b = −ln(1 − `beta`·D/K)/D (`tbode4.beta_edge_from_ode`)
  - π′ = −ln(1 − π·(1 − e^(−b·D)))/(b·D)
  - κ′ = −ln(1 − κ·(1 − e^(−b·D)))/(b·D)

  Agents on the network then transmit, per contact and per step, exactly what the homogeneous ODE does for three of the four combinations of target (susceptible or reinfected) and source (symptomatic or asymptomatic). The fourth, reinfection from asymptomatic sources (rate π′·κ′·b), can't also be matched by two adjusted parameters.
  - **Its size, from the network standalone ODE:** at monthly steps (b·D ≈ 0.12–0.28) it puts that flow about 2% low, and prevalence 0.7–3.4% low (0.966–0.993 by set); at weekly steps, 0.994–0.999.
  - **An equilibrium-matched conversion removes it at equilibrium.** Keep b and κ′, and choose π′ so that the total force of reinfection at the calibrated equilibrium equals the homogeneous model's. Every rate of change is then zero at the same state, so the network model has exactly the calibrated equilibrium. Only dynamics away from it differ, slightly.
- **This needs no change to Starsim:** agents on the network already saturate per transition. Only Starsim's ODE of a network needs the fix.

## Question

Do v4 agents on `RandomNet(10, dur=0)`, with the converted parameters, reproduce the calibrated homogeneous ODE? If they do, a calibration with homogeneous mixing can be reused for agents on a random network, without recalibrating. And how far off is converting `beta` alone, as our convention does?

## Design

**Parameter sets:** exp 21's 12 draws from exp 20's posterior (`../21_calibrated_comparison/sets.obj`), with the same equilibria.

**Arms:** agents, 50,000, on `ss.RandomNet(n_contacts=10, dur=0)`, started at each set's equilibrium, 1900–2000, 20 reps per set and arm (new seeds, from 22000):

| Arm | Per-contact `beta` | Reinfection risk | `trans_asymp` | Step | Network standalone ODE's prediction, prevalence ÷ calibrated (range over sets) |
|---|---|---|---|---|---|
| 1. Converted per combination | b | π′ (per combination) | κ′ | Monthly | 0.966–0.993 |
| 2. Converted, equilibrium-matched | b | π′ (equilibrium-matched) | κ′ | Monthly | 1 |
| 3. `beta` only (our convention) | b | π | κ | Monthly | 1.10–1.28 |
| 4. Converted, equilibrium-matched, weekly | b, κ′, π′ for weekly steps | | | Weekly | 1 |

Arm 4 checks that the conversion works at another step, with its own values. That's 960 runs, about 50 minutes on 4 cores (network runs are slower than homogeneous ones). Arm 1 tests the per-combination mechanism, including the cross term; arm 2 tests the conversion to use.

**References:**

- **The calibrated homogeneous ODE** at each set: v4's ODE for state shares and mortality, and the standalone ODE for first-entry incidence, as in exp 21.
- **Exp 21's homogeneous agents,** reused, not rerun, as the benchmark for what "matching" looks like.
- **A prediction for arm 2:** the standalone ODE with each transition by contact scaled by its per-step saturation, for each combination of target (susceptible or reinfected) and source (asymptomatic or symptomatic). This is a reference version of proposal 03's fix, added to `tbode4.fastode` as an option that is off by default. With converted parameters it should reproduce the homogeneous ODE (arm 1, apart from the cross term). With unconverted ones it predicts arm 2's gap.

**Measured:** agents ÷ homogeneous ODE, the mean over 1950–2000, for each set and arm (mean ± SE over reps), and averaged over the 12 sets. The outputs are as in exp 21:

- active prevalence and the share asymptomatic;
- first-entry incidence per prevalent case, and prevalence ÷ incidence;
- TB mortality and mortality ÷ incidence;
- latent infection (recent and remote), and the cleared share.

Plus each set's converted values (b, π′, κ′) against the originals.

**Figures:** agents ÷ ODE by set for each arm, with arm 2's prediction; and the agents' target outputs against the target intervals, for arms 1 and 2. The raw outputs behind each figure are saved with `sc.save`.

## What would answer the question

- **Arms 2 and 4 match the ODE** if every set is within about 2 SE of 1 and the average over sets is within about 1%, as exp 21's homogeneous agents did.
- **Arm 1 follows its prediction** (0.966–0.993 by set), which tests the per-combination mechanism, cross term included.
- **Arm 3 departs from the ODE** by roughly the predicted amount (1.10–1.28 by set). If it doesn't, the per-transition explanation (exp 18) is incomplete.
- **Then:** a homogeneous calibration converts to network parameters in closed form. Once Starsim's network ODE has proposal 03's fix, it should reproduce arm 1's prediction (exp 18's acceptance test).

## Not covered

- Networks whose edges persist (households): no conversion exists.
- **Starsim's own network ODE, which still lacks the per-transition fix.** It is left out on purpose rather than patched locally:
  - the v4 track's pinned starsim predates even Starsim's first network correction (`5ba36ebb`);
  - a local patch of solver internals would break the rule of pinning exact upstream commits;
  - Cliff's implementation may differ, and will need its own test.

  The standalone ODE's network option stands in for the fixed ODE. When the fix lands, rerun exp 18's `q1`/`q2` and this experiment's arms 1 and 2 with Starsim's network ODE, in an environment pinned to that commit.

## Files

- `run.py`: `conversions` → `conversions.obj`, `conversions.csv` (each set's converted parameters and the network standalone ODE's predictions); with no argument it also runs the agents → `runs.obj`, `by_set.csv`, `pooled.csv`, `agents_vs_ode.png`, `versions.json`; `summary` redraws. Uses exp 21's sets, ODE runs and outputs, and exp 19's `Incidence` analyzer.
- `v4/tbode4/fastode.py`: the network option (`net_contacts`, `net_step`; `force()`), a reference version of proposal 03's fix. Off by default; `v4/tests/test_fastode.py` passes.
- `SUMMARY.md`: findings, written after running.
