# 01 Baseline: summary

## Result

At a monthly `dt`, the ABM runs consistently below the ODE. Ratios are ABM median / ODE:

| ODE `beta` | prev_active (1950 / 2099) | INFECTION (1950 / 2099) |
|---|---|---|
| 9 (20k agents) | 0.64 / 0.65 | 0.59 / 0.56 |
| 9 (100k agents) | 0.74 / 0.44 | 0.67 / 0.41 |
| 12 | 0.77 / 0.66 | 0.72 / 0.68 |
| 20 | 0.92 / 0.88 | 0.86 / 0.84 |

- **The first few years agree.** prev_active is within about 1% in 1901, and the gap builds up as transmission does.
- **CLEARED starts about 20% high.** That holds at every `beta` in 1901, which points to the time step rather than transmission.
- **The ABM keeps drifting down after the ODE has flattened** (`beta` = 12, 1950–2100). This is consistent with the ABM behaving like an ODE with lower effective transmission, closer to threshold, where the approach to equilibrium is slow (compare the `beta` = 9 ODE).
- **Not a small-population effect.** No runs went extinct, and 100k agents at `beta` = 9 shows the same or a larger gap than 20k.
- **The ODE's demography choice matters little.** Replacing TB deaths with births (the original `TB_ODE`) instead of births = μN changes prevalence by 2–7%.

## Also established

- `beta_edge_from_ode` was first written with an assumed mean degree of 2K. Measured, `ss.RandomNet(n_contacts=K)` gives a mean degree of K, so the formula was corrected before these runs. The two-strain tbsim test's `beta_edge = BETA_ODE / K` is consistent with this.
- `TBODE` with default options reproduces `TB_ODE` (`tests/test_ode.py`).

## Next

The pattern pointed to the ABM's discrete time step; tested in [02_timestep](../02_timestep/).
