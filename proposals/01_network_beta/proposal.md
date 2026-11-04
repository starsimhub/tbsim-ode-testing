# Per-contact `beta` on networks: fixing it in Starsim and TBsim

*For Cliff, to discuss. Versions: starsim `origin/v4` @ `80c81f53` (4.0.0rc1 + 21), tbsim `origin/v4` @ `dcfb522`.*

## Summary

Your v4 redo found that on `RandomNet` at monthly `dt`, a per-contact rate of `beta/K` gives 22% less prevalence near threshold, and that the matched per-contact rate removes it (`v4/v4_findings.qmd`, arms `A0_month_net` and `A0_month_net_matched`). We found the same thing independently, measured step by step.

The cause is that each contact transmits with probability 1 − e^(−b·dt) per step, which is about b·dt/2 below the rate b. Since `beta/K` is what Starsim's ODE approximation of a network assumes, and what tbsim's docs and `ode_compare.py` lead users to, we suggest fixing it in both places:
1. **Starsim:** the ODE approximation of a network should use the per-step transmission the agents actually have.
2. **TBsim:** the docs should give the matched per-contact value, and `ode_compare.py` should use it.

## The mechanism

- **As agents,** a declared disease on a network gives each edge a probability of 1 − e^(−b·dt) per step (`solvers/agents.py`, `Transmission`).
- **As an ODE,** an individual-level network becomes homogeneous mixing at its mean number of contacts times `beta` (`solvers/mixing.py:153`, and `:538` for pools): a contact rate of K·b.

The agents' force of infection is K·(1 − e^(−b·dt))/dt per unit of infectious share. That is below K·b by a factor of (1 − e^(−b·dt))/(b·dt) ≈ 1 − b·dt/2: 4% at monthly `dt`, 1% at weekly and 0.1% at daily, at b ≈ 1 per year.

Homogeneous mixing as agents doesn't show this. Its per-step quantity is each person's chance of infection, which is small (≈ 0.004 here). On a network, the per-step quantity is each contact's chance, which is much larger (≈ 0.08) and doesn't shrink with prevalence.

## Evidence

**Step by step** ([exp 12](../../experiments/12_v4_setup/SUMMARY.md), finding 3). We used `tbsim.TB` at the equilibrium of an ODE calibrated to 500/100k (ODE `beta` = 9.92 per year): 50k agents, monthly `dt`, 4 seeds × 50 years, about 100k infections per arm. In each step we compared the infections that happened with the number the ODE's formula predicts from the agents' own counts at the start of that step:

| Agents | Infections ÷ ODE's expected (± Poisson SE) |
|---|---|
| Homogeneous mixing | 0.998 ± 0.003 |
| `RandomNet(n_contacts=10)`, b = `beta`/K | **0.956 ± 0.003** (predicted 0.960) |
| `RandomNet(n_contacts=10)`, b = −ln(1 − `beta`·dt/K)/dt | 0.997 ± 0.003 |

`n_contacts=10` gives a mean degree of 10.0, and no edge is kept from one step to the next.

**In `ode_compare.py`'s setting** ([`check_ode_compare.py`](check_ode_compare.py), results in `check_ode_compare.json`):
- `tbsim.TB` defaults (`beta` = 9, births replacing every death);
- 1% symptomatic at the start, 10k agents, `RandomNet(n_contacts=10)`, 1950–2100;
- prevalence averaged over 2050–2100.

Agents ÷ ODE prevalence, mean ± SE across runs:

| `dt` | Agents | ÷ ODE | Runs |
|---|---|---|---|
| Monthly | Homogeneous mixing | 1.07 ± 0.07 | 40 |
| Monthly | Network, `beta/K` (as the script) | **0.59 ± 0.05** | 40 |
| Monthly | Network, matched | 0.99 ± 0.06 | 40 |
| Weekly | Network, `beta/K` | 0.68 ± 0.09 | 20 |
| Weekly | Network, matched | 0.95 ± 0.10 | 20 |

So the shortfall the script's docstring attributes to a sparse network is this per-contact effect. `beta` = 9 is just above threshold, which is why a 4% difference in effective `beta` becomes 40% of prevalence. The weekly runs are too noisy to say more than that the gap is smaller.

## Suggested changes

### Starsim

- **(a) The ODE approximation uses the agents' per-step transmission:** a contact rate of K·(1 − e^(−b·dt))/dt rather than K·b. This makes the approximation depend on `dt`, but what it approximates, agents taking whole steps, does too. The existing note could say so. This is our preferred option.
- **(b) Or keep it, and have the note give the size of the bias** (about b·dt/2 in effective `beta`, worst near threshold) and the matched per-contact value.
- **(c) A helper** that gives the per-contact `beta` for an ODE `beta`, a network and a `dt`, so users porting a calibrated ODE don't divide by K. For a network with fresh contacts each step:

  ```python
  def beta_edge_from_ode(beta_ode, K, dt):
      q = beta_ode*dt/K                # Each contact's chance of transmission in a step, to match the ODE
      return -np.log1p(-q)/dt          # The per-contact rate that gives it: K·(1 − e^(−b·dt)) = beta_ode·dt
  ```

Options (b) and (c) could go with (a).

The agreement tests (`test_agents_match_ode`) use quarter-day steps, so they wouldn't catch this. For networks whose edges last longer than a step, the ODE approximation is coarser anyway, but the b·dt/2 term is separate from that and still applies.

### TBsim

- **Docs.** `tb.py` (the `beta` parameter's `doc`), `sim.py`, `docs/tutorials/tuberculosis_sim.qmd` and `docs/user_guide/examples.qmd` say that with a network, `beta` is per contact, but not which per-contact value reproduces an ODE `beta`. We suggest giving the matched value, or pointing to the Starsim helper if one is added.
- **`tests/devtests/ode_compare.py`.** We suggest:
  - using the matched value;
  - correcting the docstring's "fall below it on a sparse network";
  - considering turning it into a regression test at monthly `dt`.
