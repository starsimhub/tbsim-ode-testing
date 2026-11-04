# Per-contact `beta` on a network: docs, and `ode_compare.py`'s shortfall is a step effect, not the network

*Draft issue for starsimhub/tbsim. Versions: tbsim `origin/v4` @ `dcfb522`, starsim `origin/v4` @ `80c81f53`. The underlying cause is in Starsim: see starsimhub/starsim#NNN.*

## Summary

On a network, agents transmit with probability 1 − e^(−b·dt) per contact per step. So a per-contact `beta` of `beta_ode / K` gives about b·dt/2 less transmission than the ODE: 4% at monthly `dt` on `RandomNet(n_contacts=10)`. Starsim's ODE approximation of a network assumes `beta_ode / K` (starsimhub/starsim#NNN).

This has two consequences for tbsim:
- **The docs** say `beta` is per contact on a network, but not which per-contact value matches an ODE `beta`.
- **`tests/devtests/ode_compare.py`** passes `beta/K`. Most or all of the gap it shows, which its docstring attributes to a sparse network, is this step effect.

## `ode_compare.py`

The script passes `beta/K_CONTACTS` per contact, and its docstring says agents "fall below [the ODE] on a sparse network". We reran its setting:
- `tbsim.TB` defaults (`beta` = 9 per year, births replacing every death);
- 1% symptomatic at the start, 10k agents, `RandomNet(n_contacts=10)`, 1950–2100;
- prevalence averaged over 2050–2100.

The matched per-contact value is b = −ln(1 − `beta`·dt/K)/dt, which makes K·(1 − e^(−b·dt)) = `beta`·dt.

Agents ÷ ODE prevalence, mean ± SE across runs:

| `dt` | Agents | ÷ ODE | Runs |
|---|---|---|---|
| Monthly | Homogeneous mixing | 1.07 ± 0.07 | 40 |
| Monthly | Network, `beta/K` (as the script) | **0.59 ± 0.05** | 40 |
| Monthly | Network, matched | 0.99 ± 0.06 | 40 |
| Weekly | Network, `beta/K` | 0.68 ± 0.09 | 20 |
| Weekly | Network, matched | 0.95 ± 0.10 | 20 |

At monthly `dt`, the shortfall is the per-contact step effect, not the sparsity of the network. With the matched value, agents on a network with fresh contacts each step behave like homogeneous mixing, as they should when contacts are redrawn every step. `beta` = 9 is just above threshold, which is why a 4% difference in effective `beta` becomes 40% of prevalence. The weekly runs are too noisy to say more than that the gap is smaller.

We suggest:
- using the matched value in the script;
- correcting its docstring;
- considering turning it into a regression test at monthly `dt`, comparing agents with the ODE using the matched value.

## Docs

`tb.py` (the `beta` parameter's `doc`), `sim.py`, `docs/tutorials/tuberculosis_sim.qmd` and `docs/user_guide/examples.qmd` say that with a network, `beta` is per contact. They don't say which per-contact value reproduces an ODE `beta`, and dividing by K, which users porting a calibration will naturally do, is about b·dt/2 low. We suggest giving the matched value, or pointing to a Starsim helper if one is added (option (c) in the Starsim issue).

The script behind the table is [`proposals/01_network_beta/check_ode_compare.py`](https://github.com/starsimhub/tbsim-ode-testing/blob/main/proposals/01_network_beta/check_ode_compare.py) in [starsimhub/tbsim-ode-testing](https://github.com/starsimhub/tbsim-ode-testing), with its results in `check_ode_compare.json` in the same folder. It uses `tbode4.make_sim` (`v4/tbode4/sims.py`), where `edge_beta='rate'` is `beta/K` and `'step'` is the matched value. The step-by-step comparison in the Starsim issue is in [`experiments/12_v4_setup/`](https://github.com/starsimhub/tbsim-ode-testing/tree/main/experiments/12_v4_setup).
