# 17: Summary

**Track:** v4 (tbsim `origin/v4` @ `dcfb522`, starsim `origin/v4` @ `80c81f53`) for question 1 and the two-state cohorts of question 2; v0.10 (tbsim v0.10.3 @ `b92d7e8`, starsim v3.6.1 @ `3521adb4`) for tb_LAI_TPT as run (`versions_v4.json`, `versions_v010.json`).

## Question 1: with latent reinfection on, do v4's agents match v4's ODE?

**Design in brief:**

- **Model:** the two-state latent model: recent infections progress at r₀ and become remote at `k_latent`, and remote infections progress at φ·r₀.
  - **Parameters:** set B (tb_LAI_TPT's best-fit natural history), at φ = 0.003, 0.01 and 0.03. Progression is fitted to tb_LAI_TPT's progression targets.
  - **Reinfection:** every relative risk of reinfection is π = 0.24, including latent reinfection (remote → recent).
- **`beta`:** calibrated at each φ to 500 per 100k active prevalence at the ODE's equilibrium: 22.26, 21.14 and 19.12 per year.
- **Runs:** 50,000 agents with homogeneous mixing, started at the ODE's equilibrium and run 1900–2000. Monthly and weekly `dt`, 40 reps each (seeds 1000–1039). The φ = 0.01 monthly arm was also rerun at 6,250–200,000 agents (seeds from 5000), to test a finite-population explanation for the slightly low prevalence.
- **Measure:** agents ÷ ODE, the mean over 1950–2000.

**Finding: yes. Latent reinfection adds no gap between agents and the ODE, and prevalence and TB mortality match it at every population size.**

| Agents ÷ ODE, 1950–2000 (mean ± SE, 40 reps) | φ = 0.003, monthly | weekly | φ = 0.01, monthly | weekly | φ = 0.03, monthly | weekly |
|---|---|---|---|---|---|---|
| Active prevalence (asymptomatic + symptomatic) | 0.996 ± 0.012 | 0.992 ± 0.010 | 0.990 ± 0.010 | 0.985 ± 0.011 | 0.991 ± 0.007 | 0.991 ± 0.009 |
| Latent, recent + remote | 1.024 ± 0.009 | 0.999 ± 0.007 | 1.019 ± 0.007 | 0.994 ± 0.008 | 1.021 ± 0.005 | 0.998 ± 0.007 |
| Recent latent | 1.329 ± 0.012 | 1.064 ± 0.009 | 1.250 ± 0.010 | 1.044 ± 0.009 | 1.179 ± 0.007 | 1.033 ± 0.008 |
| Remote latent | 0.994 ± 0.009 | 0.993 ± 0.007 | 0.990 ± 0.007 | 0.988 ± 0.008 | 0.992 ± 0.005 | 0.992 ± 0.007 |
| Cleared | 0.997 ± 0.001 | 0.999 ± 0.001 | 0.997 ± 0.001 | 0.998 ± 0.001 | 0.997 ± 0.001 | 0.998 ± 0.001 |
| TB deaths | 0.989 ± 0.012 | 0.994 ± 0.011 | 0.988 ± 0.010 | 0.987 ± 0.010 | 0.988 ± 0.007 | 0.990 ± 0.010 |

All outputs and the 10–90% ranges are in `q1_summary.csv`; `q1_agents_vs_ode.png` shows each rep.

- **Prevalence is consistent with the ODE.** At 50,000 agents, pooled over the six arms, it is 0.991 ± 0.006, averaging each seed across arms because runs with the same seed are weakly correlated (r ≈ 0.17). Every arm is slightly below 1, but at 1–2 SE each. The population-size runs show no shortfall (below).
- **Not a time-step effect:** monthly and weekly arms show the same ratio, and a step effect would be about 4× larger at monthly `dt`.
- **Exp 13's set B arms** (no latent reinfection; 0.989 ± 0.006 and 0.993 ± 0.004, 200 reps each) look similar. But they were extended from 40 to 200 reps because their first 40 came out low (0.972 and 0.979), so pooling them overstates the evidence. Their other 160 reps give about 0.993 and 0.997.

**By population size** (φ = 0.01, monthly `dt`; agents ÷ ODE over 1950–2000, mean ± SE; `q1_n_summary.csv`, `q1_n_population_size.png`):

| Agents | Reps | Active prevalence | TB deaths | Susceptible | Cleared |
|---|---|---|---|---|---|
| 6,250 | 160 | 0.988 ± 0.014 | 0.985 ± 0.013 | 1.038 ± 0.012 | 0.9936 ± 0.0011 |
| 12,500 | 160 | 0.989 ± 0.010 | 0.994 ± 0.011 | 1.020 ± 0.009 | 0.9957 ± 0.0007 |
| 25,000 | 80 | 1.011 ± 0.009 | 1.007 ± 0.010 | 0.996 ± 0.008 | 0.9983 ± 0.0006 |
| 50,000 (main arm) | 40 | 0.990 ± 0.010 | 0.988 ± 0.010 | 1.015 ± 0.008 | 0.9967 ± 0.0006 |
| 200,000 | 40 | 1.008 ± 0.005 | 1.007 ± 0.005 | 0.992 ± 0.004 | 0.9989 ± 0.0003 |

- **Prevalence and TB deaths don't depend on population size.** The four new sizes pooled (inverse variance) give prevalence 1.004 ± 0.004. A finite-population bias would shrink roughly as 1/N; one large enough to give 1% at 50,000 agents would give about 8% at 6,250, against 1.2 ± 1.4% observed.
- **The cleared share is slightly low and the shortfall shrinks with N,** from 0.6% at 6,250 agents to 0.1% at 200,000, with the susceptible share moving the other way. This is most likely a finite-population effect on the long-lived states. The same small shortfall appears in exp 13 (0.997–0.999) and in v0.10's exact ABM (exp 08). Since it doesn't affect prevalence or mortality, it isn't investigated further.
- **The recent state reads high in end-of-step snapshots, as expected.** v4 dates each new infection at the end of its step, so the infections of the last step are all counted as recent, with none having had a chance to leave. The predicted excess is x/(1 − e^(−x)), with x the rate of leaving the recent state × `dt`. That gives 1.33 / 1.26 / 1.19 at monthly `dt` and 1.07 / 1.06 / 1.04 at weekly, against 1.33 / 1.25 / 1.18 and 1.06 / 1.04 / 1.03 observed.
  - **The rate of leaving** is `k_latent` + r₀ + `inf_cle`, 4.2–7.3 per year, so recent infections last 1.6–2.8 months on average.
  - **It's large** because the recent state is short-lived.
  - **The mean over a step is right.** This matters only if a calibration target is a snapshot of recent infection.

## Question 2: how does the two-state rule change what latent reinfection does, compared with tb_LAI_TPT as run?

**Design in brief:**

- **Cohort:** 250,000 agents, all latent at the start, followed for 30 years with no transmission or deaths.
- **Infection:** a constant force of infection λ (0, 0.10 and 0.60 per year), times the reinfection risk ratio π. It applies to people who have cleared or recovered, and, in the **on** arm, to latent people; in the **off** arm latent people aren't reinfected. **On − off** is what latent reinfection adds.
- **Models:**
  - **tb_LAI_TPT's model:** the progression hazard is r₀·e^(−kτ), where τ is time since infection, r₀ the initial progression rate and k the rate at which it falls. Reinfecting a latent person resets the clock (τ = 0). It appears twice, solved two ways; the rule and parameters are the same, and only the time step differs:
    - **tb_LAI_TPT as run (v0.10, monthly `dt`):** v0.10 `tbsim.TB` agents, as tb_LAI_TPT runs them. Each agent's hazard is evaluated at the start of the month and held for the whole month. Just after an infection or reinfection the true hazard falls quickly (to about 0.61 of r₀ within the first month at k = 5.85), so this overstates progression (exp 08).
    - **tb_LAI_TPT's model, solved exactly (continuous time):** the same model as an exact chain, one person's continuous-time Markov chain solved without sampling or time steps. It's what tb_LAI_TPT's agents would give with a very small `dt`.
  - **The two-state model** as v4 agents at φ = 0.003, 0.01 and 0.03, with exp 16's fits to the progression targets. Reinfection moves remote → recent; a recent infection is unchanged. Its own exact chain checks the agents.
- **Parameter sets:** tb_LAI_TPT's 12 face-validity sets. Values are means weighted by `sample_weight`, 4 seeds per set.
- **Counting:** everything is counted from each agent's transitions, so a step in which an agent moves more than once counts in full.

**Check: the v4 agents reproduce the two-state chain.** Agents ÷ chain is 0.996–1.003 for every cumulative curve and latent share, in every arm (`compare_q2.log`). So the per-person results below describe the two-state model.

**At λ = 0 the two-state cohorts are exactly exp 16's arm 1** (the same fits in `fits.csv`, each set's own cascade, `sym_dead` = 0, no reinfection). Exp 16 computed them only as exact chains, to 50 years. Here the chains agree with exp 16's to 2e-14 in curves A, B, C and the latent share, for every set and φ. The v4 agents give 0.992–1.007 of exp 16's curve B by set, at 1, 10 and 30 years (mean 0.999). What differs is tb_LAI_TPT's reference: exp 16 read it from end-of-step snapshots (100k agents, 4 seeds); here it is counted from transitions (250k agents), which is higher (finding 2).

**Finding 1: the two-state rule adds 83–93% as much progression through latent reinfection as tb_LAI_TPT's model solved exactly, and 75–84% as much as tb_LAI_TPT as run. In absolute terms the difference is small.** The first figure compares the reinfection rules themselves; the second adds v0.10's monthly-step overstatement of progression.

| Ever asymptomatic TB by 30 years, added by latent reinfection (on − off, share of the cohort) | λ = 0.10 | λ = 0.60 |
|---|---|---|
| tb_LAI_TPT as run (v0.10, monthly `dt`) | 0.00324 | 0.0366 |
| tb_LAI_TPT's model, solved exactly (continuous time) | 0.00289 | 0.0336 |
| Two-state, φ = 0.003 (agents) | 0.00267 | 0.0303 |
| Two-state, φ = 0.01 | 0.00259 | 0.0295 |
| Two-state, φ = 0.03 | 0.00248 | 0.0281 |

Weighted means over the 12 sets. SEs from seeds are 0.00001–0.00003 at λ = 0.10 and under 0.0001 at λ = 0.60. Across sets, the 25–75% range is about 0.001–0.005 at λ = 0.10. By 10 years, and for any progression, the ratios are the same within ±0.03 (`q2_summary.csv`; `q2_reinfection_adds.png`).

- **tb_LAI_TPT as run adds 6–12% more than its model solved exactly.** That is v0.10's monthly-step bias in progression (exp 08), which applies to reinfections as to first infections.
- **The two-state rule adds 7–17% less than tb_LAI_TPT's model solved exactly,** and less as φ rises. Remote infections already progress at φ·r₀ without being reinfected, so moving them back to recent gains less. Even at φ = 0.003, though, the gap is 7%. The rest comes from the two-state model's less front-loaded progression after each (re)infection (finding 3).
- **This is small against the progression itself.** At λ = 0.10, the difference from tb_LAI_TPT as run is 0.0006–0.0008 of the cohort by 30 years, about 0.5% of cumulative progression (0.13–0.14). At λ = 0.60 it's about 2% (0.006–0.009 of 0.34–0.37). That's well within the difference between the models without reinfection (finding 2).
- **Across the 12 sets, the direction holds; the relative size varies.**
  - **Direction:** the two-state rule adds less than tb_LAI_TPT as run in all 12 sets, at every φ, λ and output. It adds less than tb_LAI_TPT's model solved exactly in 10 of 12. The exceptions are small: set 0 (set B) at φ = 0.003 is 1.00–1.03, and set 10 at φ ≤ 0.01 is up to 1.08.
  - **Relative size:** two-state ÷ as run is 0.58–0.99 across sets, and two-state ÷ solved exactly 0.61–1.08; most sets are 0.65–0.95.
  - **The outlier is set 1, at about 0.6 on both ratios.** It carries the second-largest weight (0.21), and has the weakest front-loading (k = 1.0 per year, against 5.8 for set B) and low progression rates.
  - **Absolute size:** in every set the difference from tb_LAI_TPT as run is small. It's 0.1–1.0% of that set's cumulative progression at λ = 0.10, and 0.1–4.9% at λ = 0.60. In every set it's smaller than the set's own difference without reinfection (finding 2), usually by an order of magnitude.

**Finding 2: with reinfection on, the two-state model's cumulative progression is 4–15% below tb_LAI_TPT as run, as without it. Most of that is v0.10's step bias, not structure.**

| Ever asymptomatic TB, reinfection on (λ = 0.10) | 1 year | 10 years | 30 years |
|---|---|---|---|
| tb_LAI_TPT as run (v0.10, monthly `dt`), counted from transitions | 0.0611 | 0.1027 | 0.1428 |
| tb_LAI_TPT as run (v0.10, monthly `dt`), end-of-step snapshots (its `FirstPassage`) | 0.0542 | 0.0961 | 0.1344 |
| tb_LAI_TPT's model, solved exactly (continuous time) | 0.0526 | 0.0892 | 0.1247 |
| Two-state, φ = 0.003 / 0.01 / 0.03 (agents) | 0.0527 / 0.0526 / 0.0528 | 0.0945 / 0.0953 / 0.0977 | 0.1317 / 0.1329 / 0.1361 |

- **Two-state ÷ tb_LAI_TPT as run** (transitions): 0.85–0.87 at 1 year and 0.91–0.96 at 10 and 30 years, at every λ. Symptomatic TB and any progression follow the same pattern. The latent share is 2–6% higher in the two-state model (`compare_q2.log`; `q2_curves.png`).
- **The weighted gap hides a wide spread across sets, in both directions.** Without reinfection, two-state ÷ tb_LAI_TPT as run is 0.75–1.03 at 1 year, and 0.76 (set 9) to 1.12 (set 5) at 30 years, where five sets are above 1. This spread comes from exp 16's per-set fits to the progression targets, not from the reinfection rule, and reinfection doesn't change it much.
- **Against tb_LAI_TPT's model solved exactly,** the two-state model matches at 1 year and is 6–11% higher at 10–30 years: its latent hazard floors at φ·r₀ instead of declining to zero (exp 16).
- **tb_LAI_TPT's own curve B undercounts its progression.** `FirstPassage` (end-of-step snapshots, which tb_LAI_TPT's progression targets use) misses agents who reach asymptomatic TB and leave it within the same month. v0.10's step ordering lets an agent go latent → asymptomatic → non-infectious in one step. At λ = 0, ProgressionY1 is 0.0604 counted from transitions, against 0.0536 from snapshots and 0.0520 for tb_LAI_TPT's model solved exactly.
  - **Two biases nearly cancel at 1 year:** v0.10's monthly step overstates progression by about 16% (exp 08), and the snapshots lose about 11%.
  - **They don't cancel at 10 years:** 0.0800 (snapshots) against 0.0738 (solved exactly).
  - **Exp 16's comparisons** with tb_LAI_TPT as run used the snapshot curve, which is what its targets were calibrated on. By transitions, which is what drives its dynamics, the two-state model is 4–15% below tb_LAI_TPT, not 2–6%.

**Finding 3: per person, the two-state model behaves like tb_LAI_TPT in kind, but its progression comes later after each (re)infection.**

| Per person, reinfection on (weighted mean; SE ≤ 0.001); tb_LAI_TPT = as run (v0.10, monthly `dt`) | tb_LAI_TPT, λ = 0.10 | Two-state, λ = 0.10 (φ = 0.003–0.03) | tb_LAI_TPT, λ = 0.60 | Two-state, λ = 0.60 |
|---|---|---|---|---|
| Clock resets before the first progression, per person | 0.053 | 0.036–0.037 | 0.670 | 0.449–0.469 |
| Share ever progressing: with a reset ÷ without | 0.90 | 0.90 | 0.68 | 0.67 |
| Mean time from last infection or reset to first progression (years) | 0.31 | 0.42–0.51 | 0.30 | 0.40–0.45 |
| ... share within 1 year | 0.96 | 0.88–0.91 | 0.97 | 0.90–0.92 |
| Mean latent spell ending in clearance (years) | 1.37 | 1.34–1.35 | 1.32 | 1.28–1.29 |
| Share with 2 or more episodes of TB | 0.0085 | 0.0071–0.0075 | 0.082 | 0.068–0.071 |

Distributions: `q2_per_person.png`.

- **Clock resets:** the two-state model has a third fewer, by definition. Reinfecting a recent infection resets nothing, and a person must first become remote. Every latent reinfection in tb_LAI_TPT is a reset.
- **Selection is the same in both models.** People whose clock was reset while latent are less likely to progress than people whose clock never was: 0.90× at λ = 0.10 and 0.67–0.68× at 0.60. Set by set the ratio is 0.86–0.94 (λ = 0.10) and 0.63–0.72 (λ = 0.60) in both models.
  - **This isn't the effect of the reset,** which raises progression (finding 1, on − off). It compares two groups that differ in who gets into them.
  - **To be reset, a person has to stay latent** until a reinfection reaches them, at π·λ, about 0.024 per year at λ = 0.10. Most progression comes within months of infection, so most people who progress do so before any reset, and they fall in the "no reset" group. The "reset" group is people who were still latent, mostly past their high-risk period, when reinfected. In the two-state model they must also have become remote first.
  - **The reset restores some risk to that group, but not as much as the early progressors' risk.**
  - **The ratio is lower at λ = 0.60,** probably (not tested separately) because the "no reset" group also gains progression from new infections after clearance or recovery, which are more frequent at high λ but aren't resets.
  - **That the ratio matches between models** means the two-state rule selects the same kind of people for reinfection as tb_LAI_TPT's clock reset.
- **Timing is the main structural difference.** In the two-state model, progression after a (re)infection is less concentrated in the first year (88–91% against 96–97%), as at λ = 0: 0.44–0.54 years on average, against 0.31 for tb_LAI_TPT as run and 0.27 for its model solved exactly. This is the latent hazard floor exp 16 found, not the reinfection rule.
- **Repeat episodes** are 12–17% less common, in line with the lower overall progression.
- **How long latent spells last is the same in both models** (`q2_latent_spells.png`; Kaplan–Meier, with spells still latent at 30 years censored). A spell runs from a new infection to progression or clearance; a reset doesn't end it.
  - **Duration:** the mean is 1.29–1.30 years at λ = 0, 1.27–1.28 at 0.10 and 1.24–1.26 at 0.60, in both models. The share still latent at 2, 5 and 10 years is 0.19–0.20, 0.026–0.030 and 0.0013–0.0016 in all of them. At λ = 0 the agents' curves lie on their exact chains.
  - **Spells shorten slightly as λ rises** because reinfection raises the risk of progression for the rest of a spell, so some long spells end in progression instead of later clearance.
  - **Why the models are so close:** 87–90% of spells end in clearance, at the set's `inf_cle`, which both models share. Progression ends only 10–13%, and its timing differs, so the overall duration barely moves.
  - **Per person, total time latent over 30 years is distributed almost identically too** (`q2_time_latent.png`: density plots of each person's total time latent, top, and of individual spell durations, bottom). This sums all of a person's spells; a new infection after clearance or recovery starts a new one. It comes from seed 0 of each set, 250,000 agents, rerun with a per-person total added. The reruns reproduce the saved runs exactly.

    | Total time latent per person, over 30 years | λ = 0 | λ = 0.10 | λ = 0.60 |
    |---|---|---|---|
    | Median (years): tb_LAI_TPT as run / two-state | 0.75 / 0.83 | 1.50 / 1.50 | 4.67 / 4.67 |
    | 90th percentile | 3.1 / 3.1 | 4.8 / 4.8 | 9.8 / 9.7 |
    | 99th percentile | 6.8 / 6.8 | 9.6 / 9.3–9.5 | 15.0 / 14.9–15.0 |
    | Share latent for more than 5 years in total | 0.030 / 0.028–0.030 | 0.092 / 0.088–0.090 | 0.46 / 0.45 |
    | Mean (years) | 1.29 / 1.29–1.30 | 2.13 / 2.11–2.12 | 5.26 / 5.20–5.23 |

    Two-state values span φ = 0.003–0.03 where they differ. Medians and percentiles are in monthly steps.
    - **Dating:** tb_LAI_TPT as run infects at the start of a month, before that month's transitions, and its spells are dated from then. v4 dates a reinfection at the end of its month. So spells that start with a reinfection, rather than at time 0, are recorded about half a month longer in tb_LAI_TPT than in the two-state model. That's at most about 0.03 years in the means; the first spell is dated the same way in both.
    - **The shape is the same in both models:** roughly exponential at λ ≤ 0.10, with a tail of 1 in 100,000 people at about 20 years. At λ = 0.60, repeated infection keeps people latent, and the distribution peaks at 3–4 years.
    - **The visible difference is in the first months:** tb_LAI_TPT as run has slightly more very short latent spells. Its monthly step overstates progression right after infection, so more people leave latent infection for TB within a month.
    - **The number of latent spells per person is identical** (at λ = 0.10, 52% have two or more).
    - **Means agree within 1–2% at every λ,** within the half-month dating difference above.
  - **Where the models differ is how spells end:** tb_LAI_TPT's spells that end in progression do so earlier, and slightly more of them end that way. The share ending in progression is 0.109 / 0.112 / 0.125 in tb_LAI_TPT as run at λ = 0 / 0.10 / 0.60, against 0.101–0.105 / 0.103–0.107 / 0.114–0.117 in the two-state model. That's v0.10's monthly step overstating progression (finding 2).

**What this means for option B.** The two-state rule is not a material departure from tb_LAI_TPT's clock reset in a cohort. At set B's force of infection, it changes cumulative progression by about 0.5%, against a 4–15% difference without reinfection that is mostly v0.10's monthly-step bias. Latent reinfection gives no reason to reopen option B. Still open: reinfection at tb_LAI_TPT's equilibrium, with its network, demography and treatment, which a cohort can't show (README, "Not covered").

**For interpreting the φ = 0.03 arm:** that fit for set B misses ProgressionRatio by 1.65σ. With the cascade fixed, the model can't meet both targets at that φ for this set.

## Files

- `run_v4.py q1` (question 1) → `q1.obj` (calibrations, ODE runs, yearly means per rep), `q1_summary.csv`, `versions_v4.json`, `run_q1.log`
- `run_v4.py q1n` (population sizes) → `q1_n.obj` (yearly means per rep), `q1_n_summary.csv`
- `compare.py q1` → `q1_agents_vs_ode.png`; `compare.py q1n` → `q1_n_population_size.png`
- `cohort.py`: question 2's per-person outputs, used by both environments
- `run_lai_tpt.py` (question 2, v0.10) → `lai_tpt.obj` (per-run curves, distributions and counts), `versions_v010.json`, `run_lai_tpt.log`
- `run_v4.py q2` (question 2, v4) → `q2_v4.obj` (per-run outputs and the chains), `run_q2_v4.log`
- `compare.py q2` → `q2_summary.csv` (agents), `q2_chains.csv`, `q2_summary.obj` (with the mean curves and distributions behind the figures), `compare_q2.log`, `q2_reinfection_adds.png`, `q2_curves.png`, `q2_per_person.png`
- `compare.py latent` → `q2_latent_spells.obj` (the weighted Kaplan–Meier curves and exact chains behind the figure), `compare_latent.log`, `q2_latent_spells.png`; `compare.py latent plot` redraws the figure only
- `run_lai_tpt.py perperson` and `run_v4.py q2p` (seed 0 of each set's reinfection-on runs again, with each person's total time latent and number of spells) → `lai_tpt_perperson.obj`, `q2_v4_perperson.obj`; `compare.py perperson` → `q2_time_latent.obj`, `compare_time_latent.log`, `q2_time_latent.png` (add `plot` to redraw only)
