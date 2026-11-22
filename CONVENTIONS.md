# Conventions to carry over

The modelling conventions decided in this repo, written so they can be adopted in a new repo or brought into tb_LAI_TPT. Each says what to do, the code, when it's valid, the pitfalls, and where it was decided. The evidence is in the linked experiments and proposals.

**They assume tbsim v4 on starsim v4** (tbsim `origin/v4` @ `619c296`, starsim `origin/v4` @ `d83f4f87`, the versions tested here; exps 12–25 ran at tbsim `dcfb522` and starsim `80c81f53`), unless a section says otherwise. tb_LAI_TPT is on v0.10. Bringing these conventions into it means moving it to v4 (or to `tbode.TBExact` for v0.10 agents), and recalibrating (see [Progression targets](#4-progression-targets)).

| # | Convention | Status |
|---|---|---|
| 1 | [Contact networks: homogeneous mixing by default; on a random network, calibrate on Starsim's ODE of the network](#1-contact-networks-homogeneous-mixing-by-default-on-a-random-network-calibrate-on-starsims-ode-of-the-network) | Decided (2026-11-08; network ODE 2026-11-10) |
| 2 | [Latent progression: two states with remote progression](#2-latent-progression-two-states-with-remote-progression) | Tentative, for now (2026-11-05) |
| 3 | [Latent reinfection and superinfection](#3-latent-reinfection-and-superinfection) | Decided for one strain; multistrain to revisit |
| 4 | [Progression targets](#4-progression-targets) | Definitions fixed by tb_LAI_TPT |
| 5 | [Incidence: first entry to asymptomatic TB, counted from transitions](#5-incidence-first-entry-to-asymptomatic-tb-counted-from-transitions) | Decided (2026-10-13, 2026-11-05); in v4 since exp 19 |
| 6 | [Starting and comparing runs](#6-starting-and-comparing-runs) | Practice in this repo |
| 7 | [Versions](#7-versions) | Practice in this repo |
| 8 | [Calibrating the ODE: a standalone copy at equilibrium, with IMIS](#8-calibrating-the-ode-a-standalone-copy-at-equilibrium-with-imis) | Practice in this repo (2026-11-07; with ages and demography, exp 31) |
| 9 | [Treatment](#9-treatment) | Decided and implemented (2026-11-09; exp 23): fixed 6-month course; loss to follow-up on treatment to add |
| 10 | [Demography: births, ages and the starting population](#10-demography-births-ages-and-the-starting-population) | Practice in this repo (2026-11-19; exps 28–30); tested with TB and treatment (exp 30), and across the calibrated posterior (exp 34) |

## 1. Contact networks: homogeneous mixing by default; on a random network, calibrate on Starsim's ODE of the network

**Decided (user, 2026-11-08): agents mix homogeneously by default,** with no network, and use one only when a feature needs real contacts (a household network, contact-targeted interventions). Homogeneous agents match the ODE with no conversion (exps 13–17, 21, 25). In v4, homogeneous mixing is the default only when a sim has no networks. To combine it with a household network, add the homogeneous part as an explicit `ss.MixingPool` of everyone (not yet tested).

**Decided (user, 2026-11-10): on a random network, calibrate on Starsim's ODE of the network,** which is the same model as agents on it. Its parameters are per contact: b, the per-contact `beta`, and the reinfection risks as they act on a contact.

**How transmission works on `ss.RandomNet(n_contacts=K, dur=0)`,** with fresh contacts every step and K the mean degree (about K/2 edges per agent, transmitting both ways):

- **Agents:** an edge from u to v transmits in a step of length D with chance 1 − e^(−rr·κ·b·D·`rel_trans`[u]·`rel_sus`[v]·modifiers). Here rr is the transition's relative risk (1 for susceptibles, `rr_reinfection_*` for reinfection) and κ the source state's infectiousness (`trans_asymp` for asymptomatic TB, 1 for symptomatic).
- **Starsim's ODE of the network:** each transition and each source state has its own saturation factor (1 − e^(−x))/x, with x = rr·κ·b·D.
- **They agree:** per transition, agents give 0.97–1.01 of the network ODE's flows, reinfection included; at equilibrium, prevalence is 1.016 ± 0.011 of the ODE in set B and 1.005 ± 0.019 in set A at monthly steps (exp 25).

**Code here:** `tbode4.make_sim(b, n_contacts=K, edge_beta='given', solver='ode')` runs the network ODE with a per-contact b, and `solver='abm'` the matching agents. `tbode4.fastode`'s network option with `net_source_exp=1` reproduces the network ODE to within 0.4% (exp 25), but is kept as a reference, not for calibration.

**Carrying over a homogeneous calibration instead is approximate.** With q = 1 − e^(−b·D):

> b = −ln(1 − `beta`·D/K)/D, and each rr′ = −ln(1 − rr·q)/(b·D), with κ′ the same function of κ

(`tbode4.beta_edge_from_ode`, and `tbode4.rr_edge_from_ode` for each `rr_reinfection_*` and for `trans_asymp`). On exp 20's posterior at monthly steps this gives 0.978 ± 0.006 of the homogeneous ODE's prevalence (0.94–1.02 by set); weekly, 0.996 is predicted (exp 25). No choice of b, rr′ and κ′ is exact, because every factor multiplies inside the exponent. Ratios per prevalent case (the share asymptomatic, incidence) stay within 0.2%, so the error is in the level of the epidemic, not its natural history. Leaving κ unconverted gives about 2% too much prevalence instead.

- **Valid only for networks with fresh contacts every step** (`dur=0`). Networks whose edges persist, such as tb_LAI_TPT's `HouseholdNetFixed`, are untested.
- **Pitfalls:**
  - **v4's own mapping, b = `beta`/K, is low:** about b·D/2 fewer infections, 4% at monthly steps, and 22–26% less prevalence near threshold (exp 12).
  - **The step matters.** A per-contact b means something different at another D; recalibrate, or carry over through the homogeneous `beta`.
  - **The network ODE's b isn't the matched conversion:** in set B at monthly steps it is 0.94× the converted b, because reinfection at a reduced risk and asymptomatic sources saturate less than first infection from symptomatic sources (exp 25).
  - **Don't reuse tb_LAI_TPT's old `beta_ode_from_tb`.** It took −ln(1 − p) per contact, which is K·b, the low mapping above.
  - **At the v4 track's earlier pins (starsim `80c81f53`, exps 12–25), transmission differed:** only rr·b saturated, κ and the other factors scaled the chance linearly, and the network ODE applied the susceptibles' saturation to every transition, so it under-transmitted reinfection (set B ran at 1.20× the ODE's prevalence; exp 18). Converting b and each rr, with κ unchanged, was then exact (1.000 ± 0.006; exp 22). That's why `beta_edge_from_ode` and `rr_edge_from_ode` exist; they are kept so those experiments can be rerun.
- **Evidence:** [exp 25](experiments/25_starsim_network_fix/SUMMARY.md); before Starsim's fix, [exp 12](experiments/12_v4_setup/SUMMARY.md), [exp 18](experiments/18_network_ode/SUMMARY.md) and [exp 22](experiments/22_network_conversion/SUMMARY.md); [proposal 01](proposals/01_network_beta/proposal.md) and [proposal 03](proposals/03_network_reinfection/proposal.md), which Starsim's fix resolves.

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
- **After bumping the pins, also check how results are indexed.** Starsim `d83f4f87` moved a step's flows from the index of the step's end to the index of its start, while states stayed at their times; code that dates flows (e.g. yearly means) then shifts them by a step. Check that a flow at index i equals the change in its source's state from i to i + 1 (exp 26).

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
  - **Every feature added to the model** (child transmission next) goes into `fastode_age` too (and `fastode`, while it's used), and is tested against v4's ODE before it's used for calibration.
  - **After bumping the v4 pins,** rerun `v4/tests/test_fastode.py`.
- **The equilibrium shortcut needs constant rates.** With time-varying inputs (treatment scale-up, demography, trends over time), integrate trajectories with `fastode.run`, which is slower than root-finding but still far faster than Starsim. Or start from an equilibrium at an early year and integrate from there.
- **Check convergence on parameters, not only predictions:** run IMIS again with another seed and compare the marginals. In exp 24, two seeds agreed on every prediction, but not on parameters along a ridge (`k_latent` with `inf_asy`).
- **Check the weight on badly fitting points:** IMIS can give an initial prior draw that fits badly a large weight where the mixture's density is low. In exp 32, about 1% of each seed's weight sat on draws with log-likelihood more than 20 below the best, and one of them, outside 3 target intervals, was among the heaviest points (exp 34 drew it). Before using a posterior, report this weight and the log-likelihoods of resampled sets.
- **With time-varying inputs** (treatment scale-up), start each draw at its untreated equilibrium before the first change, and integrate from there (exp 24: about 90 ms per draw).

**With ages and demography (convention 10): `tbode4.fastode_age`** ([exp 31](experiments/31_fastode_ages/SUMMARY.md)), about 1.1 s per draw against about 15 s through Starsim's quarter-year ODE.

- **The model:** `run_cohort()`, exact aging by cohort, at monthly steps with the force of infection at mid-step (`foi='mid'`). Groups are one step wide and move up one group each step. Within a step, the TB transitions are one matrix exponential for every age.
  - **Accuracy:** within 0.1% of exact aging (exp 30's extrapolated reference), and within 0.06% of weekly steps on every target at 20 prior draws. It converges at second order.
  - **Not the start-of-step force of infection:** that converges at first order and is 2–4% off at monthly steps.
- **Each draw's start:** `settle(years=100)` (1950's demography frozen, no treatment, from the untreated equilibrium without ages), then `start_table()`; then 1950–2041 from that table. 100 years gives the same targets as 200 to 0.002%; 50 doesn't settle every draw.
- **Its copy of v4's ODE,** `run_rate()` (quarter-year groups aged at a rate), matches Starsim's within 0.15% in every year. `v4/tests/test_fastode_age.py` checks it against v4's ODE run live; rerun it after bumping the v4 pins.
- **Pitfalls in exact aging by cohort** (each an error proportional to the step):
  - **Rates by age over a step:** a group ages a whole width during the step, so it takes the mean of its own rate and the next group's. The lower edge's rate alone put population 1.2% low by 2020 at monthly steps.
  - **Newborns:** they need half a step of TB transitions.
  - **The step's births** belong in the predicted end-of-step population, the force of infection's denominator.
  - **Steps that put WPP's age edges on the grid:** quarter-year, monthly, 1/52 year.
- **Comparing runs at different steps:** use stocks on 1 January and flows over calendar years (`calendar()`), not yearly means of end-of-step values. Those shift by half a step, up to 0.4% when prevalence falls fast.
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
- **Dates of results:** a flow result at index i holds the transitions of the step from t_i to t_(i+1), and a state result at index i is the state at t_i (index 0 is the starting state); see convention 7. An `ss.timeseries` change applies from the first step starting on or after its date, and monthly steps (30.44 days) don't line up with calendar months.

**Starting runs:** before treatment starts, use the untreated equilibrium with nobody on treatment. A start at a treated equilibrium needs people in the treatment and `relapsing` states and, with a fixed course, a time already on treatment (uniform over the course). Otherwise start with nobody on treatment and burn in for a few years.

**To add later (user, 2026-11-09):**

- **Loss to follow-up during treatment:** leaving treatment early, back to the starting state.
- **Lower coverage for children,** with age structure.
- **Not planned:** partial benefit from incomplete treatment, and transmission on treatment.

## 10. Demography: births, ages and the starting population

**What:** with age- and sex-specific fertility and mortality (e.g. WPP), agents and v4's ODE grow alike with Starsim's own defaults, at starsim `0089ac0f` (4.0.0rc5) and later. For India, agents are within 0.1% of the ODE's population through 2050 (extrapolated to exact aging), and within about 1% of its share in each broad age group (0–14, 15–29, 30–49, 50–64, 65+). With TB and treatment on top, agents match the exactly aged ODE: population within 0.4%, ratios per prevalent case within 0.6%, and prevalence within 1% overall and by age group at weekly steps (exp 33). Across 12 sets from the calibrated posterior (exp 32), agents ÷ the exactly aged ODE prevalence is 0.996–1.004 by window at weekly steps, ratios per prevalent case within 0.25%, and the calibration targets within about 0.5% on average (exp 34).

**The setup:**

- **Both solvers:** deaths as `ss.death(rate=ss.by(age, sex))`; the starting population from `ss.population(country, year)`, which gives the fraction female at each age.
- **The ODE's births:** `ss.birth(rate=ss.by(age, sex))`, at the rate for each woman's current age group, which is how WPP's rates are defined (by age at delivery).
- **Agents' births:** `ss.Pregnancy` with WPP's table as given. `ss.birth` by age is compartmental only. Its defaults read the rates at the age at delivery (`fertility_age='delivery'`: a woman conceiving at age x takes the rate for x plus the mean gestation, 0.769 years), and start the share of women pregnant that a steady state would have.
- **Births' sex:** both solvers give births the starting population's fraction female at ages 0–4 (49.06% for India 1950), unless `p_female` is given. That stands in for WPP's sex ratio at birth, which isn't used yet. The choice matters: births at 49.06% against 48.63% girls make 1.0% more population by 2020.
- **Code:** `tbode4.demog`: `ode_demographics()`, `death_rate()`, `age_edges()` and `WPP`; the starting states by age below are in `tbode4.start` (tested in `v4/tests/test_start.py`). Its `agent_demographics()`, `shifted_fertility()`, `PregnancyStart` and `joint_sex()` were workarounds for starsim `d83f4f87` and earlier, kept to rerun exps 29–32 at their pins.

**The ODE:**

- **For calibration, `tbode4.fastode_age`'s exact aging by cohort instead** (convention 8, exp 31): no spread of ages, within 0.1% of exact aging at monthly steps. Its births' female share follows Starsim's default (`Demography.birth_female_default = 'young'`). v4's own ODE, with the groups below, stays the reference for comparisons with agents.
- **Age groups of a quarter year** (`age_edges(0.25)`). The ODE ages people by constant-rate exits, so a cohort's ages spread with variance (age × group width). With TB and treatment (exp 30), quarter-year groups are within 0.5% of exact aging in population and prevalence by age, at about 30 s per weekly run (1950–2040, one core). Single years are faster (about 22 s) but overstate population by 1.2–1.9% and put prevalence by age up to 2.7% off; WPP's 5-year groups are worse (up to 7% in population against single years, 14% in the 65+ share, 10% in prevalence by age). Keep the rates on their own groups (map each fine group to its rate group).
- **When a comparison needs better than 0.5%** (e.g. agents against the ODE), also run eighth-year groups (within 0.3%, about 45 s) and extrapolate to exact aging as 2 × eighth − quarter (exps 29–30).
- **Split by age and sex** with `sim.stratify(age=edges, sex=['female', 'male'])`: each age group starts with the population's fraction female at that age.
- **Starsim keeps every step's results by group only up to 100 MB per run;** raise `ss.solvers.Engine.group_bytes` for fine age groups at weekly steps.

**Starting TB states by age and sex** (`init=` with `ss.by(age, sex)` tables; exp 30):

- **One table by single year of age and sex for both solvers.** The ODE takes a table as coarse as its age groups or coarser. Agents need age groups no coarser than the table's bins, or they start with wrong shares when several states are given. So stratify agents by single years of age (and sex); that only affects their results and how they read the table.
- **Give every state but susceptible,** which takes the rest. Agents refuse a share by group for every state ("the initial values sum to 0"), though the ODE accepts it.
- **Draw agents' counts in each group, don't round them.** Starsim gives each state round(share × agents in the group), which loses rare states (e.g. active TB) in small old-age groups: 65+ prevalence started at 0.49–0.73 of the ODE's with 100,000 agents. While agents' states are set, swap `Compartmental.round_preserving` for a multinomial draw seeded by the run's seed: call `tbode4.use_drawn_start()` before initializing agents.
- **For a starting state with ages,** run the ODE with the starting year's demography frozen until the shares by age settle (200 years in exp 30), and average its groups to single years: `tbode4.settled_start(sim)` on the run gives the table and the `init=` tables.

**Pitfalls:**

- **Count agents by age in your own analyzer, leaving out those dying in the step.** Analyzers run before `People.step_die()` removes the dead, so counting `auids` at the end of a step includes that step's deaths (up to 0.3% of the population per month for India 1950; exps 30 and 33 counted them). Leave out agents with `ti_dead` or `ti_removed` at or before the current step (exp 34's `AgeCounts`).
- **Ignore `ss.Pregnancy`'s warning that births are made twice** with `ss.death`: a false positive. Its miscarriage and stillbirth results also miss losses at the mother's death; count them by balance if needed.
- **Results renamed at starsim `0089ac0f`:** the population's deaths are `death` (was `died`), births have a top-level `birth`, and `emigrated` is gone.
- **The ODE smooths sharp features of the age distribution** (e.g. the steps of a 5-year starting population, or the join between those alive at the start and those born after), so its broad age shares cycle by up to about 1% against agents' as these pass group edges.

**Upstream:** [proposal 06](proposals/06_pregnancy/proposal.md)'s changes (the fertility age, the start's share pregnant, the compartmental solver's sex split and births' sex, losses at the mother's death, the warning, 64-bit ages) are in starsim `0089ac0f` (PR #1625). [Proposal 07](proposals/07_starting_states/proposal.md) (draft) asks for drawn starting counts and for shares read per agent to be drawn from all states at once; until then, keep `tbode4.use_drawn_start()` and the single-year table.

**Evidence:** [exp 28](experiments/28_age_demography/SUMMARY.md) (the gap, and single years against 5-year groups); [exp 29](experiments/29_pregnancy_timing/SUMMARY.md) (its parts, and workarounds at earlier pins); [exp 30](experiments/30_tb_corrected_demography/SUMMARY.md) (with TB and treatment, the ODE's widths with TB, and the starting states); [exp 33](experiments/33_starsim_rc5_demography/SUMMARY.md) (Starsim's defaults without workarounds, demography alone and with TB); [exp 34](experiments/34_calibrated_demography_agents/SUMMARY.md) (agents across the calibrated posterior).
