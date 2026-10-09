# 02 Time step: summary

## Result

**Confirmed: the baseline gap is a time-step effect in the ABM.** It almost disappears at a weekly `dt`. Ratios are ABM median / ODE in 1959:

| Output | `beta` | 1 month | 1 week | 1 day |
|---|---|---|---|---|
| prev_active | 12 | 0.78 | 0.97 | 1.02 |
| prev_active | 20 | 0.91 | 0.99 | 1.03 |
| INFECTION | 12 | 0.74 | 0.96 | 1.00 |
| INFECTION | 20 | 0.87 | 0.98 | 1.02 |
| CLEARED (1901) | both | 1.20 | 1.03–1.05 | 0.93 |

At weekly and daily steps, the ABM's 10–90% band contains the ODE for all states (see `timestep_beta12.png`).

## Interpretation

Each `tbsim.TB` step runs transmission first, then the transitions, which are applied immediately in the order INFECTION → NON_INFECTIOUS → ASYMPTOMATIC → SYMPTOMATIC. Newly infected agents can leave INFECTION in the step they enter it, and agents can pass through several states in one step, but only in the forward direction. The mean stays in INFECTION (≈ 5.7 months) and ASYMPTOMATIC (≈ 4.7 months) are only a few monthly steps long, so the error is large at a monthly `dt`. The direction of the bias depends on the balance of flows, which explains why it varies with `beta`.

## Implications

- **Monthly `dt` (the user's tb_LAI_TPT setting) is not equivalent to the continuous-time model with the same rates.** Calibration will have absorbed this, so calibrated rates are effectively "rates for a monthly-step ABM".
- **Time-since-infection progression will make this worse.** With the calibrated `k_asy ≈ 5.85`/yr, the progression rate falls by a factor of e roughly every 2 months, which is shorter than a few monthly steps. The ordering effect on new infections (exiting INFECTION in the infection step) falls where the hazard is highest.
- **A decision for the user.** For ODE ↔ ABM consistency, either compare against the ABM at a weekly `dt`, or build the ABM's monthly-step behaviour into the ODE (e.g. discrete-time transition probabilities), or propose an upstream change to tbsim's step ordering.

## Caveats

- `tb_deaths_ppy` here is per-step deaths per person-year, which is very noisy at short steps. Future experiments should compare annual totals.
- In 1901 at daily `dt`, several states are 5–10% low with only 5 seeds. This is probably early-transient noise; not investigated.
