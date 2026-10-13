# 11 Counting incidence from transitions: is `TBExact`'s monthly incidence shortfall a snapshot artefact?

## Question

`tbode.TBExact` is a version of `tbsim.TB` that samples each agent's end-of-step state from the exact one-step probabilities `expm(Q·dt)`, so it has no step-ordering bias. It uses Poisson-mixture latent classes for progression by time since infection, and dates each new infection within its step.

Incidence is defined as in tb_LAI_TPT's IncidenceTracker: the first entry into ASYMPTOMATIC or SYMPTOMATIC per episode. An episode starts when an agent leaves SUSCEPTIBLE, INFECTION or CLEARED, so cycling between ASYMPTOMATIC and NON_INFECTIOUS within an episode is counted once. The ODE counts this exactly as a flow, INFECTION → ASYMPTOMATIC plus NON_INFECTIOUS → ASYMPTOMATIC for agents not yet active in the episode. The ABM's `StateCounts` analyzer counts it from end-of-step snapshots.

At monthly `dt`, `TBExact`'s incidence is consistently about 3% further below the ODE than its prevalence is. For example, set B with linear mass-action transmission has prevalence 1.000 ± 0.006 and incidence 0.970 ± 0.005 (exp 10). At weekly `dt` the gap mostly disappears. The hypothesis, never tested, is a snapshot artefact. An agent who enters ASYMPTOMATIC and leaves again for NON_INFECTIOUS within the same step is never seen as active. If that agent then recovers without being observed active again, the episode is never counted.

This matters because incidence and notifications will be calibration targets. **Does counting incidence from the sampled transitions, rather than from snapshots, remove the shortfall?**

## Change tested: event-based incidence in `TBExact`

For each agent that starts a step in INFECTION or NON_INFECTIOUS, `TBExact` draws whether its path passed through ASYMPTOMATIC or SYMPTOMATIC during the step, given its start and end states:

- **Augmented chain.** The probability comes from an augmented version of the step's continuous-time chain, in which every state has an "unmarked" and a "marked" copy. Entering ASYMPTOMATIC or SYMPTOMATIC moves the path to the marked copy, and it stays marked. Then `P(visited, end = j | start = i)` is the marked entry of `expm(Q_aug·dt)`, and the probability of having visited given the end state is that divided by `P(end = j | start = i)`. For new infections, the step is shortened to the time remaining after infection, as for the end state. Each latent class has its own augmented chain.
- **The dynamics don't change.** The end state is sampled exactly as before; the visit is drawn afterwards from a separate random stream (`tb_rng_visit`). So runs are bit-identical to exp 10's.
- **Episode bookkeeping.** `TBExact` keeps a per-agent flag for "active in this episode" and counts a new episode when the path visited ASYMPTOMATIC or SYMPTOMATIC, the agent started in INFECTION or NON_INFECTIOUS, and the flag was off. The flag resets when the agent ends the step in CLEARED. At most one episode can start per step, because infection only happens at the start of a step.
- **Outputs.** The count is stored in a new result, `new_episodes`. `abm_outputs` reports it as `incidence_event_py`, alongside the snapshot-based `incidence_py`.

## Design

- **Arms:** `TBExact` with `'within'` timing, 25k agents, started at the ODE equilibrium, 100 years (1900–2000), each rep averaged over its last 50 years. Settings and calibrated `beta` as in exp 10:

  | # | Setting | `dt` | Transmission | Reps | Seeds also run in |
  |---|---|---|---|---|---|
  | 1 | Set B | monthly | mass action, linear | 190 | exp 10 arm 5 |
  | 2 | Set B | monthly | network | 190 | exp 10 arm 4 |
  | 3 | Set A | monthly | mass action, linear | 190 | exp 10 arm 3 |
  | 4 | Set B | weekly | mass action, linear | 190 | (new) |

  Set A is tbsim defaults + `k` = 5.846, near threshold; set B is tb_LAI_TPT's natural history + `k`, the reference set (exp 08).
- **Bit-identity check:** for arms 1–3, every rep's prevalence series must equal exp 10's exactly, confirming the counting doesn't perturb the dynamics.
- **Comparison:** for each arm, event-based and snapshot-based incidence ÷ ODE (mean ± SE across reps), and their paired ratio, which is very precise because it comes from the same runs. These are compared with prevalence ÷ ODE.
- **Predictions, if the hypothesis is right:**
  - **Event-based incidence ÷ ODE tracks prevalence ÷ ODE:** about 1.00 for arms 1 and 4, about 0.99 for arm 2, about 0.98 for arm 3.
  - **Event ÷ snapshot is about 1.03 at monthly `dt`** and about 1.01 or less at weekly.
- **Unit test (`tests/test_exact.py`):** in one step, the fraction of agents starting in INFECTION or NON_INFECTIOUS (with the episode flag off) that are counted matches `Σ_j expm(Q_aug·dt)[start, marked j]`, computed independently.
- **Outputs:** per-rep time series (`outputs.obj`), `summary.csv`, and a figure of the snapshot and event-based incidence trajectories against the ODE (`incidence_counting.png`).

## Versions

tbsim main @ `b92d7e8` (v0.10.3), starsim @ `3521adb4` (v3.6.1). `TBExact` gains event-based incidence counting, which leaves its dynamics unchanged.
