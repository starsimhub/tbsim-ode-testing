# 11 Counting incidence from transitions: summary

## Background

`tbode.TBExact` is a version of `tbsim.TB` that samples each agent's end-of-step state from the exact one-step probabilities `expm(Q·dt)`, so it has no step-ordering bias. Incidence follows tb_LAI_TPT's IncidenceTracker: the first entry into ASYMPTOMATIC or SYMPTOMATIC per episode. Until now the ABM counted it from end-of-step snapshots (`StateCounts`).

At monthly `dt`, that count was about 3% further below the ODE than prevalence was. The hypothesis was a snapshot artefact: an agent who enters ASYMPTOMATIC and leaves again within the step, and is never observed active later in the episode, is never counted.

`TBExact` now also counts incidence from its sampled transitions (`incidence_event_py`). For each agent starting a step in INFECTION or NON_INFECTIOUS, it draws whether the path visited ASYMPTOMATIC or SYMPTOMATIC during the step, given its start and end states. The probability comes from an augmented chain in which entering either state marks the path. The draw uses its own random stream after the end state is drawn, so the dynamics don't change.

## Results

All arms are `TBExact` with `'within'` timing and 25k agents, run for 100 years from the ODE equilibrium; each rep is averaged over its last 50 years. Values are ABM ÷ ODE (1.00 = match), mean ± SE over 190 reps. The last two columns are paired: they compare the two counts within the same runs, which makes them very precise.

| Arm | Prevalence | Incidence, snapshot | Incidence, from transitions | Transitions ÷ snapshot | Incidence per prevalent case (transitions), ABM ÷ ODE |
|---|---|---|---|---|---|
| Set B, monthly, mass action, linear | 1.000 ± 0.006 | 0.970 ± 0.005 | **1.002 ± 0.005** | 1.0332 ± 0.0003 | 1.003 ± 0.001 |
| Set B, monthly, network | 0.989 ± 0.005 | 0.960 ± 0.005 | **0.991 ± 0.005** | 1.0328 ± 0.0002 | 1.003 ± 0.001 |
| Set A, monthly, mass action, linear | 0.984 ± 0.010 | 0.954 ± 0.010 | **0.984 ± 0.010** | 1.0318 ± 0.0003 | 1.002 ± 0.001 |
| Set B, weekly, mass action, linear | 1.001 ± 0.006 | 0.993 ± 0.005 | **1.000 ± 0.005** | 1.0075 ± 0.0001 | 1.000 ± 0.002 |

Set A is tbsim defaults + `k` = 5.846 (near threshold); set B is tb_LAI_TPT's natural history + `k` (the reference set). Both are from exp 08; transmission options as in exp 10.

- **The hypothesis is confirmed: the shortfall was a snapshot artefact.** Counted from transitions, incidence moves with prevalence in every arm. Incidence per prevalent case, which is the quantity the counting method affects, matches the ODE to within 0.3%.
- **The snapshot misses 3.2–3.3% of incident episodes at monthly `dt` and 0.75% at weekly**, the same in every arm and setting. That's a fixed property of the step length and the rates out of ASYMPTOMATIC, not of transmission. The arms come out within 0.1% of each other, with paired SEs of 0.01–0.03%.
- **The counting doesn't change the dynamics.** In arms 1–3, all 190 reps' prevalence series are bit-identical to exp 10's.
- **A small residual:** incidence per prevalent case at monthly `dt` is 0.2–0.3% above the ODE (z ≈ 2). It's too small to matter for calibration, and probably reflects the force of infection being fixed at its start-of-step value. It hasn't been investigated.
- **Unit test (`tests/test_exact.py`):** in one monthly step, the fraction of INFECTION and NON_INFECTIOUS starters counted as new episodes matches an independently computed probability of visiting ASYMPTOMATIC or SYMPTOMATIC within the step, within 4 SE. All agents ending the step active are counted, and some agents ending inactive are too: exactly the episodes the snapshot misses.

`incidence_counting.png` shows the two counts against the ODE over time. The snapshot count sits a steady 3% below the event count from the start.

## Conclusions

- **Use `incidence_event_py` (counted from transitions) wherever incidence is compared with the ODE or used as a calibration target.** The snapshot count (`incidence_py`) understates incidence by 3.3% at monthly `dt` and 0.75% at weekly, whatever the parameters.
- **The same approach will be needed for notifications once treatment is added.** Treatment starts are triggered by interventions between natural-history steps, so they may be countable directly. But any target defined by passing through a state within a step needs either a transition-based count or the augmented-chain draw used here.
- **tb_LAI_TPT's own IncidenceTracker counts from snapshots with `tbsim.TB`.** There, the step-ordering bias means agents moving backward stop for the step, so ASYMPTOMATIC → NON_INFECTIOUS within a step can't happen, and the artefact should be smaller. That hasn't been measured. If tb_LAI_TPT is recalibrated with `TBExact`, its incidence targets should be compared with the transition-based count.

## Note: which incidence definition for calibration (open, user, 2026-10-13)

**Decision (user):** use the transition-based count (`incidence_event_py`) for every ODE↔ABM comparison. It's the ABM's equivalent of the ODE's incidence flow, and it doesn't depend on `dt`. The snapshot count is a measurement that changes with the step length, not a model quantity.

**Deferred:** whether calibration targets should exclude brief subclinical episodes, for example an agent who enters ASYMPTOMATIC and leaves again within a month.
- **The snapshot count isn't a way to exclude them.** It drops an episode only if the episode falls entirely between two step boundaries, so whether a short episode counts depends on its timing and on `dt`, not on its length.
- **If short episodes should be excluded, define the rule by duration** (e.g. active for at least one month in total per episode) and implement it in both models. The ABM would track time spent active; the ODE would need extra compartments for it, e.g. a chain or phase-type structure.
- **Whether to exclude them depends on what the target data measure.** Subclinical incidence isn't observed directly. For notifications or prevalence-survey targets, it may be better to model those observations directly (treatment starts; point prevalence) than to adjust incidence.
- **tb_LAI_TPT's current incidence already effectively counts every entry.** It uses snapshots with `tbsim.TB`, where an agent can't move backward within the step it entered ASYMPTOMATIC, so every visit lasts at least until the end of the step and is seen. This follows from exps 04–05 but hasn't been measured.

## Versions

tbsim main @ `b92d7e8` (v0.10.3), starsim @ `3521adb4` (v3.6.1). `TBExact` gained event-based incidence counting (results `new_episodes`; states `been_active`, `new_episode`); its dynamics are unchanged (bit-identical runs).
