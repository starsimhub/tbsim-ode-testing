# ODE approximation of a network assumes more transmission than agents have (by about b·dt/2)

*Draft issue for starsimhub/starsim. Versions: starsim `origin/v4` @ `80c81f53` (4.0.0rc1 + 21), tbsim `origin/v4` @ `dcfb522`. Related tbsim issue: starsimhub/tbsim#NNN (docs and `ode_compare.py`).*

## Summary

On a network, each contact transmits in a step with probability 1 − e^(−b·dt), where b is the per-contact `beta`. When the same model runs as an ODE, the network is replaced by homogeneous mixing at `n_contacts × b`, which treats b as a rate. The two differ by a factor of (1 − e^(−b·dt))/(b·dt) ≈ 1 − b·dt/2. For TB on `RandomNet(n_contacts=10)` at monthly `dt`, agents get 4% fewer infections than the ODE assumes.

Near the epidemic threshold this becomes a large difference in prevalence. In the setting of tbsim's `tests/devtests/ode_compare.py`, agents reach 0.59 ± 0.05 of the ODE's prevalence with `beta/K` per contact, against 0.99 ± 0.06 with the matched value (details in the tbsim issue).

## Background

We are checking how closely tbsim's agents reproduce its ODE, feature by feature, for a TB model calibrated to tb_LAI_TPT (India, monthly `dt`). In v4 the comparison is between the two solvers of one declared `tbsim.TB`. The v4 ODE reproduces our independent ODE and `TB_ODE` to 1e-5 or better, so the ODE side is not in question here.

## What happens now

- **As agents,** a declared disease on a network gives each edge a probability of 1 − e^(−b·dt) per step (`solvers/agents.py`, `Transmission`). A susceptible person with c infectious contacts is infected with probability 1 − e^(−c·b·dt).
- **As an ODE,** an individual-level network becomes homogeneous mixing at its mean number of contacts times `beta` (`solvers/mixing.py:153`, and `:538` for pools): a contact rate of K·b.

The force of infection the agents actually have is K·(1 − e^(−b·dt))/dt per unit of infectious share, which is below K·b by about b·dt/2.

## Evidence (agents, measured step by step)

We used `tbsim.TB` at the equilibrium of an ODE calibrated to 500/100k prevalence (ODE `beta` = 9.92 per year): 50k agents, monthly `dt`, 4 seeds × 50 years, about 100k infections per arm. In each step we compared the infections that happened with the number the ODE's formula predicts from the agents' own counts at the start of that step. This measures transmission alone, independent of any drift in the dynamics.

| Agents | Infections ÷ ODE's expected (± Poisson SE) |
|---|---|
| Homogeneous mixing | 0.998 ± 0.003 |
| `RandomNet(n_contacts=10)`, b = `beta`/K | **0.956 ± 0.003** (predicted: (1 − e^(−b·dt))/(b·dt) = 0.960) |
| `RandomNet(n_contacts=10)`, b = −ln(1 − `beta`·dt/K)/dt | 0.997 ± 0.003 |

We also checked the network itself: `n_contacts=10` gives a mean degree of 10.0, and no edge is kept from one step to the next (`dur=0`).

**Why homogeneous mixing doesn't show this.** Its per-step quantity is each person's chance of infection, x = `beta`·w·dt, which is small (here x ≈ 0.004, so x/2 ≈ 0.2%). On a network, the per-step quantity is each contact's chance, b·dt, which is much larger (here 0.08). It doesn't shrink with prevalence.

## Why it matters

The bias is about b·dt/2 in effective `beta`: 4% at monthly `dt`, 1% at weekly, 0.1% at daily, at b ≈ 1 per year. Models near threshold amplify it, as TB models calibrated to low prevalence often are. At tbsim's default `beta` (9, just above threshold), the 4% becomes a 40% shortfall in prevalence. The agreement tests (`test_agents_match_ode`) use quarter-day steps, so they wouldn't catch it.

## Options

- **(a) The ODE approximation uses the agents' per-step transmission:** a contact rate of K·(1 − e^(−b·dt))/dt rather than K·b. This makes the approximation depend on `dt`, but what it approximates, agents taking whole steps, does too. The existing note could say so. This is our preferred option.
- **(b) Keep the approximation, and have its note give the size of the bias** (about b·dt/2 in effective `beta`, worst near threshold) and the matched per-contact value.
- **(c) A helper** that gives the per-contact `beta` for an ODE `beta`, a network and a `dt`, so users porting a calibrated ODE don't divide by K. For a network with fresh contacts each step:

  ```python
  def beta_edge_from_ode(beta_ode, K, dt):
      q = beta_ode*dt/K                # Each contact's chance of transmission in a step, to match the ODE
      return -np.log1p(-q)/dt          # The per-contact rate that gives it: K·(1 − e^(−b·dt)) = beta_ode·dt
  ```

Options (b) and (c) could go with (a).

**Scope.** This is about the per-step edge probability, so it applies to any disease on any network as agents. For networks whose edges last longer than a step (`dur > 0` or `dynamic=False`), the ODE approximation is coarser anyway, because of repeated contacts and saturation of each person's neighbours. The b·dt/2 term is separate from that and still applies.

The scripts and outputs behind the tables are in [starsimhub/tbsim-ode-testing](https://github.com/starsimhub/tbsim-ode-testing):
- **The step-by-step comparison:** [`experiments/12_v4_setup/`](https://github.com/starsimhub/tbsim-ode-testing/tree/main/experiments/12_v4_setup). See `run.py` (`check_transmission`), `SUMMARY.md` (finding 3) and `summary.json` (`B2_transmission`).
- **The `ode_compare.py` rerun:** [`proposals/01_network_beta/check_ode_compare.py`](https://github.com/starsimhub/tbsim-ode-testing/blob/main/proposals/01_network_beta/check_ode_compare.py), with its results in `check_ode_compare.json`.
