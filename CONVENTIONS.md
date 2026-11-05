# Conventions to carry over

The modelling conventions decided in this repo, written so they can be adopted in a new repo or brought into tb_LAI_TPT. Each says what to do, the code, when it's valid, the pitfalls, and where it was decided. The evidence is in the linked experiments and proposals.

**They assume tbsim v4 on starsim v4** (tbsim `origin/v4` @ `dcfb522`, starsim `origin/v4` @ `80c81f53`, the versions tested here), unless a section says otherwise. tb_LAI_TPT is on v0.10. Bringing these conventions into it means moving it to v4 (or to `tbode.TBExact` for v0.10 agents), and recalibrating (see [Progression targets](#4-progression-targets)).

| # | Convention | Status |
|---|---|---|
| 1 | [Per-contact `beta` on a random network](#1-per-contact-beta-on-a-random-network) | Decided (2026-11-05) |
| 2 | [Latent progression: two states with remote progression](#2-latent-progression-two-states-with-remote-progression) | Tentative, for now (2026-11-05) |
| 3 | [Latent reinfection and superinfection](#3-latent-reinfection-and-superinfection) | Decided for one strain; multistrain to revisit |
| 4 | [Progression targets](#4-progression-targets) | Definitions fixed by tb_LAI_TPT |
| 5 | [Incidence: first entry to asymptomatic TB, counted from transitions](#5-incidence-first-entry-to-asymptomatic-tb-counted-from-transitions) | Decided (2026-10-13, 2026-11-05); not yet in v4 |
| 6 | [Starting and comparing runs](#6-starting-and-comparing-runs) | Practice in this repo |
| 7 | [Versions](#7-versions) | Practice in this repo |

## 1. Per-contact `beta` on a random network

**Decided (user, 2026-11-05):** match the ODE's infections per step, in our own code, not in base Starsim or TBsim, for now.

On `ss.RandomNet(n_contacts=K, dur=0)`, with fresh contacts each step, each contact transmits in a step with probability 1 − e^(−b·dt), where b is its per-contact `beta`. The ODE `beta` and b correspond when:

> `beta_ode` = K · (1 − e^(−b·dt)) / dt

```python
import numpy as np

def beta_edge_from_ode(beta_ode, n_contacts, dt):
    """ Per-contact beta (per year) for agents on ss.RandomNet that matches an ODE beta; dt in years """
    q = beta_ode*dt/n_contacts  # Each contact's chance of transmission in a step
    return -np.log1p(-q)/dt

def beta_ode_from_edge(beta_edge, n_contacts, dt):
    """ The inverse: the ODE beta that matches a per-contact beta, e.g. one calibrated in an ABM """
    return n_contacts*(-np.expm1(-beta_edge*dt))/dt
```

- **Code here:**
  - `tbode4.beta_edge_from_ode` and `tbode4.beta_ode_from_edge` (`v4/tbode4/sims.py`); `tbode4.make_sim` uses the first by default (`edge_beta='step'`).
  - On v0.10: `tbode.beta_edge_from_ode`.
  - Tested in `v4/tests/test_anchor.py::test_beta_edge_round_trip`.
- **Which direction:**
  - **Starting from an ODE `beta`** (porting a calibrated ODE to agents): `beta_edge_from_ode`.
  - **Starting from an ABM's per-contact `beta`** (as tb_LAI_TPT calibrates): `beta_ode_from_edge` gives the ODE equivalent. The ABM itself doesn't change.
- **Valid only:**
  - for networks with fresh contacts every step (`dur=0`), which behave like mass action;
  - with K the **mean degree**. That's what `n_contacts` is on `ss.RandomNet`: about K/2 edges per agent, with transmission both ways. Checked on starsim 3.6.1 and v4.
- **Not covered: networks whose edges persist,** such as tb_LAI_TPT's `HouseholdNetFixed`. They have no ODE counterpart, and their `beta` (tb_LAI_TPT's `beta * rel_beta_hh`) isn't converted this way.
- **Pitfalls:**
  - **v4's own mapping, `beta/K`, is low.** It's what Starsim v4's ODE approximation of a network assumes, and what tbsim's docs suggest. It gives about b·dt/2 fewer infections: 4% at monthly `dt`, and 22–26% less prevalence near threshold.
  - **`dt` matters.** A per-contact `beta` calibrated at one `dt` doesn't mean the same at another; convert through the ODE `beta`.
  - **Don't reuse tb_LAI_TPT's old `beta_ode_from_tb`.** It used −ln(1 − p) per contact, which is K·b, v4's low mapping. It was in `ode_testing/run_ode_vs_abm.py`, deleted in tb_LAI_TPT commit `bb19486` (2026-10-28).
- **Evidence:** [exp 12](experiments/12_v4_setup/SUMMARY.md) (0.956 ± 0.003 of the ODE's infections per step with `beta/K`, 0.997 ± 0.003 matched); [exp 13](experiments/13_v4_pinned_rerun/SUMMARY.md); [proposal 01](proposals/01_network_beta/proposal.md).

## 2. Latent progression: two states with remote progression

**Tentative (user, 2026-11-05): option B of the [options note](proposals/02_latent_progression/options.md), for now.** v4's two latent states, where remote infections progress at a fraction φ of the recent rates:

| From | To NON | To ASY | To remote | To cleared |
|---|---|---|---|---|
| Recent (`infection`) | `inf_non` | `inf_asy` | `k_latent` | `inf_cle` |
| Remote (`infection_remote`) | `inf_non_remote` = φ · `inf_non` | `inf_asy_remote` = φ · `inf_asy` | — | `inf_cle` |

```python
phi = 0.01  # Remote progression ÷ recent; swept for now
tb = tbsim.TB(inf_non=inf_non, inf_asy=inf_asy, k_latent=k_latent,
              inf_non_remote=phi*inf_non, inf_asy_remote=phi*inf_asy, ...)
```

- **φ > 0, always.** With φ = 0 a remote infection never progresses, so the model is effectively a single latent state with raised clearance, outside tb_LAI_TPT's prior for clearance. Set aside (user, 2026-11-04).
- **φ for now:** sweep a couple of values, about 0–0.03 (up to 0.1 if clearance may rise within its prior; [exp 16](experiments/16_two_state_vs_lai_tpt/SUMMARY.md)). Calibrate it eventually; it may only be identifiable with all the targets, burden data included.
- **Fitting:** fit the recent rates and `k_latent` to the progression targets (convention 4), not by converting tb_LAI_TPT's `k`. The two models aren't the same, and tb_LAI_TPT's `k` was fitted with v0.10's monthly bias.
- **Pitfalls:**
  - **Set both remote rates explicitly.** v4's defaults are `inf_non_remote=None` (the same as `inf_non`, so no front-loading to NON) and `inf_asy_remote=0`.
  - **Our φ is not v4's `phi`,** which is the treatment-failure rate.
  - **`k_latent=0`** (v4's default) means no remote state at all.
- **What would reopen it:** see the [options note](proposals/02_latent_progression/options.md#decision-tentative). Examples: the full recalibration can't fit the data within the priors, or the underlying behaviour lacks face validity against tb_LAI_TPT's.

## 3. Latent reinfection and superinfection

**Mirror tb_LAI_TPT: every relative risk of reinfection is π** (≈ 0.24 in tb_LAI_TPT's best fit). In v4:

| Parameter | Who | v4 default | Set to |
|---|---|---|---|
| `rr_reinfection_cleared`, `rr_reinfection_rec`, `rr_reinfection_treat` | Cleared, recovered, treated | 1.00, 0.21, 3.15 | π |
| `rr_reinfection_inf` | Remote latent (→ recent) | 0.21 | π |
| `rr_superinfection_latent` (with strains) | Recent latent, new strain | 1.0 | π |
| `rr_superinfection_non` (with strains) | Non-infectious TB, new strain | 1.0 | π |

- **The rule in v4:** reinfecting a remote infection makes it recent again. Reinfecting a recent one leaves it recent, which is what a clock reset means for a state without memory. Whether the rule matters with φ > 0 is exp 17.
- **Multistrain: to revisit later** (user, 2026-11-05). The aim is v0.10 `TBResistant`'s behaviour where possible. The table above gives its single σ; what still differs is in the [options note](proposals/02_latent_progression/options.md#with-multiple-strains-drug-resistance), for example that v4 has no strain counts.

## 4. Progression targets

tb_LAI_TPT's targets, from household-contact studies (Ferebee 1970, Sutherland 1968):

- **ProgressionY1** = B(1) = 0.0523, log-space σ 0.094;
- **ProgressionRatio** = B(1) / B(10) = 0.663, log-space σ 0.071.

B(t) is the share of a cohort, all latent at t = 0, that has ever reached asymptomatic TB, with no transmission, reinfection or deaths. z = log(model ÷ target) ÷ σ.

- **Compute them exactly in v4:**
  - `tbode4.first_passage` or `tbode4.progression_targets` (`v4/tbode4/latent.py`), a deterministic chain;
  - or v4's own ODE cohort, `tbode4.make_sim(beta=0, solver='ode', init=dict(infection=1.0), mu=0, replace_tb_deaths=True, ...)`, which matches the chain to 2.3e-6.
- **tb_LAI_TPT computed them with v0.10's monthly ABM** (`cohort_progression.py`). That overstates early progression (exp 08), so its calibrated values don't carry over as is. At the same parameters, the exact model gives 8% less lifetime progression ([exp 16](experiments/16_two_state_vs_lai_tpt/SUMMARY.md)).
- **Weighting:** fit to the data matters much more than consistency with tb_LAI_TPT as run (user, 2026-11-05).

## 5. Incidence: first entry to asymptomatic TB, counted from transitions

**Decided (user):**

- **Count from transitions,** not from end-of-step snapshots (2026-10-13).
- **Count only the first entry to asymptomatic TB in each episode,** per person-year, as tb_LAI_TPT does (2026-11-05). Re-entries after ASY → NON → ASY don't count.

**Why:**

- **Transitions:** at monthly `dt`, snapshots miss 3.3% of episodes, because they're brief ([exp 11](experiments/11_incidence_counting/SUMMARY.md)).
- **First entry:** this is tb_LAI_TPT's definition (`IncidenceTracker`). Re-entries could be a sizeable share of all entries, because people leave ASY for NON quickly (`asy_non` = 1.66 per year by default); not measured.

**Code:**

- **v0.10, done:**
  - **ODE:** `tbode.TBODE`'s `incidence_py` counts entries to ASY from latent infection and from NON0, the non-infectious people not yet asymptomatic in this episode (`tbode/ode.py`). NON0 is a split of NON with the same rates, so it doesn't change the dynamics.
  - **Agents:** `tbode.TBExact`'s `new_episodes` (`incidence_event_py` in `tbode.compare`).
- **v4, not yet:**
  - **What's there:** `v4/common.py`'s `outputs()` counts **every** entry to ASY from latent or non-infectious TB, on both solvers. So it suits agent ↔ ODE checks, but it isn't this definition. Don't use it for incidence.
  - **The likely fix:** split `non_infectious` the same way, into "not yet asymptomatic this episode" and "has been", with identical rates. Both solvers then run the same declared model, so both can count first entries. An agents-only flag wouldn't do, because the ODE couldn't count the same thing.
- **Open:** whether calibration targets should leave out brief subclinical episodes. If so, define the rule by duration in both models, not by step boundaries.
- **Notifications,** once treatment is on, need a transition-based count too.

## 6. Starting and comparing runs

- **Start agents at the ODE's equilibrium** (`init=` in v4, exact) rather than burning in. Near threshold the approach to equilibrium is slow.
- **Compare proportions, not counts.** Use the same outputs on both sides: state shares, prevalence (ASY + SYM), TB mortality, and incidence as in convention 5 (not yet available in v4). ABM infection counts include reinfections.
- **Reps:** show the mean ± SE across seeds and a 10–90% interval, with enough reps that the SE is small against the gap being tested. Near threshold, the spread across reps is about 14% of prevalence at 50k agents.
- **v4 seeds are deterministic:** reps with the same seed and code repeat a run, so extra reps need new seeds.
- **Step:** v4 agents match the ODE at monthly `dt`. v0.10's `tbsim.TB` doesn't, because of its step ordering and start-of-step hazards (exps 03–08); use `tbode.TBExact` there.

## 7. Versions

- **Pin tbsim and starsim to exact git commits** in `pyproject.toml` / `uv.lock`, not editable checkouts shared with other projects. tbsim v4 still reports `__version__ = '0.10.3'`, so only the commit identifies it.
- **Record the versions with every run's outputs** (`tbode4.save_versions`).
- **After bumping the pins,** rerun an anchor test like `v4/tests/test_anchor.py` (the v4 ODE reproduces the reference ODE).
