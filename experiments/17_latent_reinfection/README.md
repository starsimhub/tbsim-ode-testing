# 17: Latent reinfection with v4's two-state latent model

**Track:** both environments; results pass between them as saved files.

- **v4** (`v4/.venv`; tbsim `origin/v4` @ `dcfb522`, starsim `origin/v4` @ `80c81f53`): the two-state model, as agents, as an ODE and as a deterministic chain (`tbode4.latent`).
- **v0.10** (repo root `.venv`; tbsim v0.10.3 @ `b92d7e8`, starsim v3.6.1 @ `3521adb4`): tb_LAI_TPT as run.

## Background

- **tb_LAI_TPT** runs tbsim v0.10's agent-based `tbsim.TB` at monthly `dt`.
  - **Progression:** latent infections progress at a rate r₀·e^(−kτ) that falls with time since infection τ.
  - **Reinfection:** reinfecting a latent person resets τ to 0.
- **v4 has no time since infection.** The latent structure chosen for it, tentatively, is two states, **recent** and **remote**:
  - recent infections progress at r₀ and become remote at `k_latent`;
  - remote ones progress at φ·r₀;
  - both clear at `inf_cle`.
- **v4's reinfection rule:** a reinfected remote infection becomes recent again. A recent one stays recent; the state has no memory, so that is a reset, but there is no declining hazard within it to restart.

Latent reinfection hasn't been tested in v4.

## Questions

With latent reinfection on:

1. **Do v4's agents match v4's ODE?** They should, since every rate is constant.
2. **How does the two-state model change what reinfection does, compared with tb_LAI_TPT as run?** This compares the two agent-based models.
   - **Why agents:** they show what happens to each person (how often they're reinfected, and what follows), which an ODE can't.
   - **Why it's valid:** if question 1 holds, v4's agents run the same model as its ODE.

## Design

### Common to both questions

- **φ = 0.003, 0.01 and 0.03,** a range in which the two-state model meets tb_LAI_TPT's progression targets with tb_LAI_TPT's other rates unchanged (exp 16). φ is swept, not calibrated.
- **Progression targets:** the share of a latent cohort that has reached asymptomatic TB by 1 year is 0.0523, and by 1 year ÷ by 10 years is 0.663.
  - **The fits:** at each φ, the recent rates and `k_latent` are fitted to these, keeping the parameter set's other rates.
  - **Source:** they're exp 16's fits for each set (`../16_two_state_vs_lai_tpt/fits.csv`).
- **Reinfection:** every relative risk of reinfection is π (0.24), as in tb_LAI_TPT. That covers people who have cleared latent infection, recovered from non-infectious TB or been treated, and latent people (`rr_reinfection_inf`, remote → recent). Non-infectious TB isn't reinfected: v4 doesn't, and in v0.10 it changes nothing.

### Question 1: v4 agents vs v4 ODE

- **Parameters:** set B (tb_LAI_TPT's best-fit natural history), with `beta` calibrated at each φ to 500 per 100k active prevalence (asymptomatic + symptomatic) at the ODE's equilibrium.
- **Population:** births replace every death (`mu` = 1/70); homogeneous mixing; no treatment.
- **Runs:**
  - 50,000 agents, started at the ODE's equilibrium, 1900–2000;
  - monthly and weekly `dt`, 40 reps each: 240 runs, under an hour on 4 cores.
- **Measure:** agents ÷ ODE, the mean over 1950–2000, with mean ± SE and the 10–90% range across reps. For active prevalence, latent infection (recent and remote), the cleared and recovered shares, and TB mortality.
- **If there's a gap,** rerun that arm with latent reinfection off, to separate reinfection from remote progression.
- **Population size,** added after the first results: prevalence came out about 1% below the ODE in every arm, at 1–2 SE each. To test for a finite-population effect, which would shrink roughly as 1/N, the φ = 0.01 monthly arm was rerun at 6,250, 12,500, 25,000 and 200,000 agents (160, 160, 80 and 40 reps; new seeds), against the same ODE run.

### Question 2: the two agent-based models, in a cohort

**Cohort:** everyone latent at the start, followed for 30 years, with no transmission and no deaths. People are infected at a constant force of infection λ, times their relative risk:

- **cleared or recovered people:** at π·λ, in every arm;
- **latent people:** at π·λ in the **on** arm, by each model's own rule, and not at all in the **off** arm.

The difference, on − off, is what latent reinfection adds: the extra progression caused by reinfecting latent people. The two models' progression without reinfection already differs slightly (exp 16), so comparing this difference between them isolates the reinfection rule. The two arms share random seeds, so most run-to-run noise cancels in the difference.

λ takes three values:

- **0:** no reinfection at all, so the on and off arms are the same and only one is run. This is the baseline for the per-person outputs, and the setting of exp 16, whose results the v4 agents should reproduce;
- **0.10 per year:** about set B's force of infection at 500 per 100k with the two-state model (0.092–0.108 across φ, question 1);
- **0.60:** a stress test.

**Parameter sets:** tb_LAI_TPT's 12 face-validity sets, the 12 most likely samples of its posterior (set B is the first). Results are weighted by `sample_weight`, renormalised over the 12.

**Models:**

1. **tb_LAI_TPT as run:** v0.10 `tbsim.TB`, monthly `dt`, built as tb_LAI_TPT's `cohort_progression.py` builds its cohort, with each set's calibrated parameters and π.
   - **Infection:** a subclass of `tbsim.TB` (`TBConstantFOI`) infects each eligible agent each step with probability 1 − e^(−π·λ·`dt`), at the point in the step where transmission happens, through `tbsim.TB`'s own `set_prognoses`. So a reinfected latent agent's clock resets, as in transmission.
   - **Check:** at λ = 0 the cohort is identical, draw for draw, to tb_LAI_TPT's own (`cohort_progression._run_once`).
2. **The two-state model as v4 agents,** at each φ, monthly `dt`, with each set's exp 16 fit (recent rates and `k_latent` fitted to the progression targets, the set's own other rates). Infection is declared as constant-rate transitions at π·λ (cleared → recent, recovered → recent, and remote → recent in the on arm), in a subclass of `tbsim.TB` (`TBCohort`).
3. **Exact chains** (`tbode4.latent`): one person's continuous-time Markov chain, solved exactly, for the two-state model at each φ (an exact check of the v4 agents), and for **tb_LAI_TPT's model, solved exactly (continuous time)**: r₀·e^(−kτ) with a clock reset, the same rule and parameters as tb_LAI_TPT as run. The two differ only in the time step; tb_LAI_TPT's agents approximate it at monthly `dt`. So it separates v0.10's step effects from the difference in structure.

**Counting:** every output is counted from the transitions each agent makes, so a step in which an agent moves more than once is counted in full. v0.10's `tbsim.TB` can move an agent through several states in one step (e.g. latent → asymptomatic → non-infectious); `TBConstantFOI` records each move. For v4, each transition's `.uids` gives the agents who took it. tb_LAI_TPT's own curve B (`FirstPassage`, from end-of-step snapshots) is recorded alongside. **Timing:** each event is dated at the end of its step, except infections in v0.10, which `TBConstantFOI` applies at the start of the step, before the step's transitions, and which are dated there. v4 dates infections at the end of their step, so spells that start with a reinfection are recorded about half a month shorter in v4 than in v0.10.

**Size:** 250k agents × 4 seeds per set, λ, arm and (for v4) φ; the on and off arms share seeds (in v4, same-seed on and off runs agree on whether 93% of agents ever progress, against 53% between seeds). That's 240 v0.10 runs (5 minutes on 4 cores) and 720 v4 runs (40 minutes).

**Outputs,** for each set, model, φ and λ:

- **What latent reinfection adds** (λ > 0): the share that has reached asymptomatic TB by 10 and 30 years, on − off. Two-state ÷ tb_LAI_TPT.
- **Progression with reinfection on:** cumulative curves to 30 years, two-state ÷ tb_LAI_TPT, for any progression (non-infectious TB or worse), asymptomatic or worse, and symptomatic.
- **Latent infection:** the share still latent over time (recent and remote for v4).
- **Per person,** at every λ, so λ = 0 is the baseline. Both models use the same definitions, written once (`cohort.py`) and used in both environments. "Latent" is recent + remote in v4. A **clock reset** is a latent reinfection that restarts the latent course: in tb_LAI_TPT every latent reinfection; in the two-state model remote → recent (reinfecting a recent infection changes nothing).
  - clock resets before the first progression (within 30 years);
  - progression among people with a reset, against those without;
  - the time from the last infection or reset to the first progression (distribution);
  - the duration of latent spells, by how they end (progression or clearance);
  - episodes per person: progressing again after recovery and a later reinfection.
- **Chain check:** v4 agents ÷ chain, for the cumulative curves and the latent shares, at every λ. At λ = 0 this also checks v4's agents against exp 16's results.
- **Summaries:** weighted mean and 25–75% range across the 12 sets; SEs across seeds.

### Figures

- **Question 1:** agents ÷ ODE by φ and `dt`.
- **What latent reinfection adds,** against φ, with tb_LAI_TPT as a reference line, at each λ.
- **Cumulative curves with reinfection on,** one panel per φ, against tb_LAI_TPT.
- **The per-person distributions,** two-state against tb_LAI_TPT.

The raw outputs behind each figure are saved with `sc.save`.

### What would answer the questions

- **Question 1:** agents ÷ ODE within about 2 SE of 1 at every φ and `dt` (SE about 1%). A larger gap points to a bug or a solver problem.
- **Chain check:** v4's agents match the chain within their SE. Otherwise the per-person results don't describe the two-state model.
- **Question 2:** how far the two-state model is from tb_LAI_TPT, in what latent reinfection adds and in the per-person outputs. The differences are judged against the effect of reinfection itself, and against the 2–6% difference in cumulative progression without reinfection (exp 16).
  - **If they're material,** for example in the direction or size of the effect at plausible φ, that counts against the two-state model on consistency with tb_LAI_TPT ([options note](../../proposals/02_latent_progression/options.md)).

## Files

- **`cohort.py`** (both environments): the per-person outputs, from each step's events.
- **`run_lai_tpt.py`** (v0.10 environment): `TBConstantFOI`, the analyzer and the cohorts (`check` for the λ = 0 check). Writes `lai_tpt.obj`, `versions_v010.json`, `run_lai_tpt.log`.
- **`run_v4.py`** (v4 environment): question 1 (`q1` → `q1.obj`, `q1_summary.csv`; `q1n`, the population sizes → `q1_n.obj`, `q1_n_summary.csv`), and question 2 (`q2`: `TBCohort`, the agent cohorts and the chains → `q2_v4.obj`, `run_q2_v4.log`). Uses `sc.parallelize` with `ncpus=tbode4.ncpus`. Writes `versions_v4.json`.
- **`compare.py`** (v4 environment): summaries and figures from the saved files (`q1` → `q1_agents_vs_ode.png`; `q1n` → `q1_n_population_size.png`; `q2` → `q2_summary.csv` (agents), `q2_chains.csv`, `q2_summary.obj` (with the mean curves and distributions behind the figures), `compare_q2.log`, `q2_reinfection_adds.png`, `q2_curves.png`, `q2_per_person.png`; `latent` → `q2_latent_spells.obj`, `compare_latent.log`, `q2_latent_spells.png`, the duration of latent spells and how they end; `perperson` → `q2_time_latent.obj`, `compare_time_latent.log`, `q2_time_latent.png`, each person's total time latent).
- **Per-person reruns,** added after the first analysis: seed 0 of each set's reinfection-on runs again, with each person's total time latent and number of latent spells added to `cohort.py` (`run_lai_tpt.py perperson` → `lai_tpt_perperson.obj`; `run_v4.py q2p` → `q2_v4_perperson.obj`). Each reproduces its saved run exactly (curve B checked).
- **`SUMMARY.md`:** findings, written after running.

## Not covered

- **Reinfection at tb_LAI_TPT's equilibrium:** that would need its full model, with its network, demography and treatment.
- **Multiple strains:** later.
- **The household network,** which concentrates reinfection within households: later.
