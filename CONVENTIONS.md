# Conventions to carry over

The modelling conventions decided in this repo, written so they can be adopted in a new repo or brought into tb_LAI_TPT. Each says what to do, the code, when it's valid, the pitfalls, and where it was decided. The evidence is in the linked experiments and proposals.

**They assume tbsim v4 on starsim v4** (tbsim `origin/v4` @ `dcfb522`, starsim `origin/v4` @ `80c81f53`, the versions tested here), unless a section says otherwise. tb_LAI_TPT is on v0.10. Bringing these conventions into it means moving it to v4 (or to `tbode.TBExact` for v0.10 agents), and recalibrating (see [Progression targets](#4-progression-targets)).

| # | Convention | Status |
|---|---|---|
| 1 | [Per-contact `beta` and reinfection risks on a random network](#1-per-contact-beta-and-reinfection-risks-on-a-random-network) | Decided (2026-11-05; reinfection risks 2026-11-08) |
| 2 | [Latent progression: two states with remote progression](#2-latent-progression-two-states-with-remote-progression) | Tentative, for now (2026-11-05) |
| 3 | [Latent reinfection and superinfection](#3-latent-reinfection-and-superinfection) | Decided for one strain; multistrain to revisit |
| 4 | [Progression targets](#4-progression-targets) | Definitions fixed by tb_LAI_TPT |
| 5 | [Incidence: first entry to asymptomatic TB, counted from transitions](#5-incidence-first-entry-to-asymptomatic-tb-counted-from-transitions) | Decided (2026-10-13, 2026-11-05); in v4 since exp 19 |
| 6 | [Starting and comparing runs](#6-starting-and-comparing-runs) | Practice in this repo |
| 7 | [Versions](#7-versions) | Practice in this repo |
| 8 | [Calibrating the ODE: a standalone copy at equilibrium, with IMIS](#8-calibrating-the-ode-a-standalone-copy-at-equilibrium-with-imis) | Practice in this repo (2026-11-07) |
| 9 | [Treatment](#9-treatment) | Decided and implemented (2026-11-09; exp 23): fixed 6-month course; loss to follow-up on treatment to add |

## 1. Per-contact `beta` and reinfection risks on a random network

**Decided (user, 2026-11-05):** match the ODE's infections per step, in our own code, not in base Starsim or TBsim, for now.

**Decided (user, 2026-11-08): convert each reinfection risk too, and leave `trans_asymp` unchanged**, until Starsim's ODE of a network handles reinfection ([proposal 03](proposals/03_network_reinfection/proposal.md)). See "Reinfection risks" below.

**Decided (user, 2026-11-06): v4 agents mix homogeneously for the ODE recalibration.** With homogeneous mixing, agents match the ODE without any conversion (exps 13–17), so this section applies only to work on `RandomNet`: running a homogeneous calibration as agents on the network. Starsim's ODE of a network (starsim `5ba36ebb`) converts `beta` as we do but not the reinfection risks, so it under-transmits reinfection relative to agents ([exp 18](experiments/18_network_ode/SUMMARY.md)). To combine homogeneous mixing with a household network later, add the homogeneous part as an explicit `ss.MixingPool` of everyone: in v4, homogeneous mixing is the default only when a sim has no networks (not yet tested).

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
  - **Converting `beta` alone matches only transmission at the susceptibles' rate.** Reinfection at a reduced risk rr runs on the same contacts at rr·b, and agents then transmit more than the ODE's rr·`beta`, because a smaller per-contact rate saturates less. In set B (rr = π = 0.24) at monthly `dt`, reinfection is 7–8% higher in agents than in the ODE, and prevalence is 1.20× the ODE's (1.05× weekly; exp 18). Over exp 20's posterior it's 1.07–1.26× (exp 22). Starsim's ODE of a network (starsim `5ba36ebb`) has the same gap, and is otherwise the same model as the `beta` conversion. Hence the reinfection risks below.
  - **Don't reuse tb_LAI_TPT's old `beta_ode_from_tb`.** It used −ln(1 − p) per contact, which is K·b, v4's low mapping. It was in `ode_testing/run_ode_vs_abm.py`, deleted in tb_LAI_TPT commit `bb19486` (2026-10-28).
- **Evidence:** [exp 12](experiments/12_v4_setup/SUMMARY.md) (0.956 ± 0.003 of the ODE's infections per step with `beta/K`, 0.997 ± 0.003 matched); [exp 13](experiments/13_v4_pinned_rerun/SUMMARY.md); [proposal 01](proposals/01_network_beta/proposal.md).

**Reinfection risks.** As agents, a factor written into a transmission's rate saturates per contact, and every other factor scales a contact's chance linearly. On a contact in a step D:

- **The rate saturates:** for the transition's most infectious source state, a contact transmits with chance 1 − e^(−rr·b·D).
- **Linear factors multiply that chance:** another source state's coefficient (κ, `trans_asymp`, for asymptomatic TB), per-agent `rel_trans` and `rel_sus`, and modifiers on the transmission.

So convert what saturates and leave what's linear. With q = 1 − e^(−b·dt), set each `rr_reinfection_*` to

> rr′ = −ln(1 − rr·q)/(b·dt)

so that a reinfection's chance per contact is rr·q, and leave `trans_asymp` unchanged. Every combination of target and source then transmits what the homogeneous model's flows give over a step: q, κ·q, rr·q and κ·rr·q. The other parameters don't involve contacts and carry over unchanged.

```python
def rr_edge_from_ode(rr, beta_edge, dt):
    """ The relative risk that gives agents on ss.RandomNet the ODE's reinfections per step; beta_edge from beta_edge_from_ode, dt in years """
    q = -np.expm1(-beta_edge*dt)
    return -np.log1p(-rr*q)/(beta_edge*dt)

def rr_ode_from_edge(rr_edge, beta_edge, dt):
    """ The inverse: the ODE's relative risk that matches agents with a given rr' """
    return np.expm1(-rr_edge*beta_edge*dt)/np.expm1(-beta_edge*dt)
```

- **Code here:** `tbode4.rr_edge_from_ode` and `tbode4.rr_ode_from_edge` (`v4/tbode4/sims.py`), tested in `v4/tests/test_anchor.py::test_rr_edge_round_trip`. For the standalone ODE, the network option of `tbode4.fastode` (`net_contacts`, `net_step`) gives the network model with unconverted parameters.
- **Result:** a homogeneous calibration converted this way reproduces the calibrated ODE on `RandomNet(10, dur=0)`: prevalence 1.000 ± 0.006 monthly and 0.999 ± 0.006 weekly over 12 posterior sets, with every ratio per prevalent case within 0.2% ([exp 22](experiments/22_network_conversion/SUMMARY.md)).
- **Pitfalls:**
  - **Don't convert `trans_asymp`.** Agents apply it linearly, so converting it as if it saturated cuts transmission from asymptomatic TB, and prevalence comes out 2–5% low at monthly `dt` (exp 22).
  - **Factors applied linearly need no conversion.** That includes products that act on `rel_sus` or `rel_trans` (e.g. TPT or vaccines) and `ss.modify` of a transmission. Only parameters written into a transmission's rate, here `rr_reinfection_*`, are converted.
  - **rr′ depends on b and `dt`.** Recompute it whenever `beta` or `dt` changes. Starting from an agent model calibrated per contact (as tb_LAI_TPT is), `rr_ode_from_edge` gives the ODE's relative risk.
  - **It depends on how Starsim's agents transmit** (checked at starsim `80c81f53`; the same code at `bc2704b2`). After bumping Starsim, rerun exp 22's `diagnose.py`, which measures each class's infections per step against both rules.
- **When to revisit:** once Starsim's ODE of a network applies a factor per transition (proposal 03), calibrating on the network ODE with unconverted risks is the same model, parametrised per contact. Then choose between that and calibrating with homogeneous mixing and converting.

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
- **φ for now:** sweep a couple of values, about 0–0.03 (up to 0.1 if clearance may rise within its prior; [exp 16](experiments/16_two_state_vs_lai_tpt/SUMMARY.md)). Calibrate it eventually; it may only be identifiable with all the targets, burden data included. Exp 20's calibration without treatment (prior 0–0.1, five targets: burden and progression) didn't identify it. The likelihood is flat for φ from 0.01 to 0.1, so its posterior (0.008–0.095) reflects the prior.
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

- **The rule in v4:** reinfecting a remote infection makes it recent again. Reinfecting a recent one leaves it recent, which is what a clock reset means for a state without memory. With φ = 0.003–0.03 it adds 83–93% as much progression as tb_LAI_TPT's model solved exactly (continuous time). Against tb_LAI_TPT as run (v0.10, monthly `dt`), which overstates progression after each infection, it's 75–84%. Either way that's about 0.5% of cumulative progression at set B's force of infection (exp 17).
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
- **Count relapses after treatment as new episodes** (2026-11-09), as tb_LAI_TPT and WHO do. A relapse goes straight to symptomatic TB, so the first-entry rule alone would miss it. A failed course continues the same episode. Implemented in [exp 23](experiments/23_treatment/README.md).

**Why:**

- **Transitions:** at monthly `dt`, snapshots miss 3.3% of episodes, because they're brief ([exp 11](experiments/11_incidence_counting/SUMMARY.md)).
- **First entry:** this is tb_LAI_TPT's definition (`IncidenceTracker`). Re-entries are about half of all entries (46% in set A, 50% in set B at 500/100k), because people leave ASY for NON quickly (`asy_non` ≈ 1.7 per year). So counting every entry roughly doubles incidence and halves prevalence ÷ incidence ([exp 19](experiments/19_v4_incidence/SUMMARY.md)).

**Code:**

- **v0.10, done:**
  - **ODE:** `tbode.TBODE`'s `incidence_py` counts entries to ASY from latent infection and from NON0, the non-infectious people not yet asymptomatic in this episode (`tbode/ode.py`). NON0 is a split of NON with the same rates, so it doesn't change the dynamics.
  - **Agents:** `tbode.TBExact`'s `new_episodes` (`incidence_event_py` in `tbode.compare`).
- **v4, done ([exp 19](experiments/19_v4_incidence/SUMMARY.md)):**
  - **ODE:** `tbode4.fastode`, a standalone copy of v4's ODE with NON split into NON0 and NON1 (identical rates). `flows(x, p).incidence` counts first entries; `summary()` gives the calibration outputs. It reproduces v4's ODE (`v4/tests/test_fastode.py`; rerun after bumping the v4 pins). v4's own declared model still has a single NON state.
  - **Agents:** the `Incidence` analyzer in `experiments/19_v4_incidence/run.py`, with a per-agent flag for "has been asymptomatic this episode". It matches the ODE per prevalent case within 0.4%.
  - **Don't use** `v4/common.py`'s `outputs()` for incidence. It counts every entry to ASY, and with the two-state model it also misses remote latent → ASY.
- **Open:** whether calibration targets should leave out brief subclinical episodes. If so, define the rule by duration in both models, not by step boundaries.
- **Notifications** are counted from transitions too ([convention 9](#9-treatment)).

## 6. Starting and comparing runs

- **Start agents at the ODE's equilibrium** (`init=` in v4, exact) rather than burning in. Near threshold the approach to equilibrium is slow.
- **Compare proportions, not counts.** Use the same outputs on both sides: state shares, prevalence (ASY + SYM), TB mortality, and incidence as in convention 5. ABM infection counts include reinfections.
- **Reps:** show the mean ± SE across seeds and a 10–90% interval, with enough reps that the SE is small against the gap being tested. Near threshold, the spread across reps is about 14% of prevalence at 50k agents.
- **v4 seeds are deterministic:** reps with the same seed and code repeat a run, so extra reps need new seeds.
- **Step:** v4 agents match the ODE at monthly `dt`. v0.10's `tbsim.TB` doesn't, because of its step ordering and start-of-step hazards (exps 03–08); use `tbode.TBExact` there.

## 7. Versions

- **Pin tbsim and starsim to exact git commits** in `pyproject.toml` / `uv.lock`, not editable checkouts shared with other projects. tbsim v4 still reports `__version__ = '0.10.3'`, so only the commit identifies it.
- **Record the versions with every run's outputs** (`tbode4.save_versions`).
- **After bumping the pins,** rerun an anchor test like `v4/tests/test_anchor.py` (the v4 ODE reproduces the reference ODE).

## 8. Calibrating the ODE: a standalone copy at equilibrium, with IMIS

**What:** calibrate the ODE with `tbode4.fastode`, a standalone copy of v4's TB ODE (our own code, outside Starsim's solver). Its outputs come from its endemic equilibrium:

- **A short burn-in** (300 years with an adaptive solver), to land near the endemic equilibrium rather than the one where TB has died out.
- **Then root-finding** for the state where every rate of change is zero.

Sample the posterior with IMIS.

**Why:** one equilibrium takes about 10 ms, against 0.1–0.3 s through Starsim, which steps a full sim for thousands of years. [Exp 20](experiments/20_ode_calibration/SUMMARY.md)'s calibration (228,000 model runs) took about 10 minutes on 4 cores, against roughly 2–5 hours through Starsim. The equilibria are the same: they agree to 1e-14 ([exp 19](experiments/19_v4_incidence/SUMMARY.md)).

**Workflow** (`experiments/20_ode_calibration/calib.py` and `run.py`, reusable):

1. **Coverage check:** about 5,000 prior draws. Does each target's interval fall inside the spread of the model's outputs, and do some draws satisfy every target at once? A target out of reach is a model or prior problem, not something calibration can fix (exp 20's mortality ÷ incidence without treatment).
2. **Re-identification:** IMIS on synthetic targets from a known parameter set. It should recover the truth's best fit, and cover the truth in most marginals. A truth outside its interval can mean the parameter isn't identified; check whether the likelihood is flat in it (φ in exp 20).
3. **Real targets.**

**Settings and likelihoods:**

- **IMIS:** 10,000 initial prior draws; 1,000 draws per added component; stop at (1 − 1/e)·3,000 expected unique points in a resample of 3,000. Allow 300 iterations: broad, weakly identified posteriors can need more than 150.
- **Likelihoods:** Binomial for survey proportions, log-Normal (data σ) for ratios and the progression targets. The ODE is deterministic, so there's no simulation-noise term.
- **Progression targets:** from `tbode4.latent`'s exact chain (convention 4).
- **No endemic equilibrium:** a draw without one gets zero likelihood (`fastode.equilibrium` returns None).

**Keeping it valid:**

- **It is a second implementation of the model.** v4's own ODE stays the reference for ODE↔agents comparisons.
  - **Every feature added to the model** (treatment next) goes into `fastode` too, and is tested against v4's ODE before it's used for calibration.
  - **After bumping the v4 pins,** rerun `v4/tests/test_fastode.py`.
- **The equilibrium shortcut needs constant rates.** With time-varying inputs (treatment scale-up, demography, trends over time), integrate trajectories with `fastode.run`, which is slower than root-finding but still far faster than Starsim. Or start from an equilibrium at an early year and integrate from there.
- **Check identifiability:** compare each marginal with its prior. With few targets, many parameters stay near their priors (exp 20: `beta` and π constrained; φ and the non-infectious cascade not).

## 9. Treatment

**Status:** design agreed (user, 2026-11-09), to implement in our own code: a `tbsim.TB` subclass in `tbode4`, and `fastode`. Implemented as `tbode4.TBTreat` (v4's ODE and agents) and in `tbode4.fastode`; [exp 23](experiments/23_treatment/SUMMARY.md) checks it, and compares the design with tb_LAI_TPT and with tbsim v4 as it stands. The course is fixed (user, 2026-11-09): it gives about 1% less prevalence than a rate, through TB death on treatment, and v4's ODE runs it exactly. The parts that may belong in tbsim v4 itself are in [proposal 04](proposals/04_treatment/proposal.md).

**The design:**

- **Starting treatment:** a declared transition from symptomatic TB at θ(t) = −ln(1 − c(t)) per year, where c is the annual coverage (the share of people with symptomatic TB treated within a year). Make it time-varying with `ss.timeseries` and `ss.modify`, in both solvers. Using c itself as the rate understates treatment by 31–40% at c = 0.55–0.68.
- **The outcome, drawn at the start:** a person starting treatment enters `treatment_success` with probability s = adherence × efficacy, or `treatment_fail` otherwise. A time-varying efficacy then applies at its value when treatment starts.
- **Failure, back to the starting state:** one pair of treatment states for each state treatment can start from. That's exact in both solvers, with no per-person memory. Starting from symptomatic TB only, there's one pair.
- **The course:** fixed, D = 0.5 years (user, 2026-11-09): `TBTreat(course='fixed')`, the default. In `fastode`, pass `tx_stages=12`: Erlang stages approximate the fixed course, within 0.1% of prevalence (exp 23); the default, 1, is a rate.
- **After cure:** v4's `treated` state, reinfected at `rr_reinfection_treat` = π like every other reinfection (convention 3), not v4's default of 3.15 (user, 2026-11-09). Keeping a separate state costs nothing, keeps the number previously treated visible, and leaves room for a different risk later.
- **Relapse:** a share `p_relapse` of successes goes to `relapsing` instead of `treated`, then to symptomatic TB at rate 1/1.5 per year. People waiting to relapse can't be reinfected.
- **On treatment:** nobody transmits or is reinfected. People whose treatment will fail keep the TB death rate of the state they started from (`sym_dead`) for the course; people who will be cured don't die of TB (user, 2026-11-09). Count these deaths in TB mortality, and include them in any births that replace TB deaths (`tbsim.TB`'s replacement covers symptomatic TB only).
- **Notifications:** every start of treatment, counted from the transitions into treatment, so restarts after failure or relapse count again. Report notifications ÷ incidence too.

**Declaring it in Starsim v4:**

- **A plain duration is a rate:** `dur=ss.years(0.5)` is the mean of an exponential dwell time, in both solvers. A fixed time needs `dur=ss.constant(ss.years(0.5))`. Agents then leave after exactly D, and the ODE runs it as slots of whole steps, exact to the step (or as k Erlang stages, with `stages=k`).
- **A fixed exit can branch:** `ss.Transition(treatment_success, [treated, relapsing], dur=ss.constant(...), p=[None, p_relapse])` runs in both solvers.
- **No flags:** keep every part of treatment as states and transitions, not per-person flags such as `diagnosed`, `pending_relapse` or `prior_state`. Stale flags caused most of v0.10's treatment pitfalls (re-initiation, test results carried into later courses).
- **Restarting in the step of failure is correct:** v4 agents can fail and restart within one step, for a share 1 − e^(−θ·dt/2) of failures, as in continuous time.
- **Dates of results:** results at index i hold the transitions of the step ending at i (index 0 is the starting state). An `ss.timeseries` change applies from the first step starting on or after its date, and monthly steps (30.44 days) don't line up with calendar months.

**Starting runs:** before treatment starts, use the untreated equilibrium with nobody on treatment. A start at a treated equilibrium needs people in the treatment and `relapsing` states and, with a fixed course, a time already on treatment (uniform over the course). Otherwise start with nobody on treatment and burn in for a few years.

**To add later (user, 2026-11-09):**

- **Loss to follow-up during treatment:** leaving treatment early, back to the starting state.
- **Lower coverage for children,** with age structure.
- **Not planned:** partial benefit from incomplete treatment, and transmission on treatment.
