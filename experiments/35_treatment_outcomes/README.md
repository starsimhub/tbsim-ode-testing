# 35: Death and loss to follow-up on treatment

**Status:** designed, not yet run.
**Track:** v4 (`v4/.venv`; tbsim `origin/v4` @ `a2ae583`, starsim `origin/v4` @ `0089ac0f`).
**Run on:** zebra for the agents (about 1 hour on 100 workers, as exp 34); the rest on the laptop.

## Why this experiment

**What treatment is now.** It follows [CONVENTIONS.md 9](../../CONVENTIONS.md#9-treatment); [exp 23](../23_treatment/SUMMARY.md) built it.

- **Two outcomes:** each course is drawn at its start as one that will cure (probability s, "efficacy") or one that will fail.
- **Failing courses** face TB death at `sym_dead` for the whole 6 months, then return to symptomatic TB.
- **Missing:** nobody who would be cured dies on treatment, and nobody leaves treatment early.
- **What efficacy absorbs:** calibrated efficacy is low (0.4–0.65 in 1995 across exp 32's posterior, against reported treatment success near 85%). So the "fail" outcome stands in for deaths, losses to follow-up (LTFU) and true failures together.

**What it should become.** Death and LTFU were set aside as "to add later" (user, 2026-11-09). This experiment adds them so that:

- **They match reported data:** deaths and losses on treatment match WHO's reported treatment outcomes for India.
- **People lost to follow-up can benefit partly** from the part of the regimen they took.
- **Regimens can be compared later:** a more effective or shorter regimen should change cures, losses and deaths through the model's structure, not through outcome shares fixed by hand.

## The design (user, 2026-11-22 and 2026-11-23)

### A course: two tracks, in phases

**The track is drawn when a course starts.**

- **Cure track:** a completed regimen would cure, with probability c, the regimen's efficacy among people who complete it.
- **Fail track:** a completed regimen would fail, with probability 1 − c.

**The course is split into phases.** The fixed course of length D (6 months) becomes M consecutive states of fixed length w (the phase width), with M = D/w. Here w = 2 months, so M = 3. Phases let the hazards below, and the benefit of a partly taken regimen, depend on how far into the course someone is. Each is exact in both solvers.

- **Set by width, not by count** (user, 2026-11-23). Regimens of different lengths are then resolved alike. With 2-month phases, a 6-month regimen has 3 and a 4-month one 2. With a fixed count, their phases would differ in width, and that difference would leak into a comparison of partial benefit or of front-loaded deaths. Two months also matches the intensive phase of current regimens.
- **Limits:** a phase must be at least one step long (v4's ODE runs fixed durations in whole steps), so 1-month phases are the finest possible at monthly steps.

**Three things can happen in phase k,** as competing exits:

| Exit | Rate | Then |
|---|---|---|
| Death on treatment from TB | μ_c,k on the cure track, μ_f,k on the fail track | Death, counted in TB mortality |
| Death from other causes | WPP's background rate for the person's age and sex, as in every state | Death, counted as a background death |
| Loss to follow-up | h on both tracks | Cure track: cured with probability q_k, otherwise back to the starting state. Fail track: back to the starting state |
| End of the phase | after D/M | The next phase. At the end of the last phase, the cure track is cured (to `treated`, or `relapsing` with share `p_relapse`); the fail track returns to the starting state |

### Death: a hazard by track and phase

Each track has its own death hazard in each phase. This replaces TB death at `sym_dead` on failing courses; `sym_dead` still applies to symptomatic TB off treatment. Deaths from other causes on treatment stay at the background rate (user, 2026-11-23): every excess death on treatment is in the TB hazards, and so in TB mortality.

**Why a hazard and not a share drawn at the start.** The model should let a better regimen reduce deaths through its structure. Two simpler forms don't:

- **A died share drawn at the start, independent of the track:** this matches WHO's outcome categories directly. But nothing about the regimen changes it: neither higher efficacy nor a shorter course would reduce deaths.
- **Deaths as a share of failures** (died = δ(1 − c)): a more effective regimen does reduce deaths here. But every death falls among people whose regimen would fail, and that stops fitting when efficacy is high.
  - **Example:** with c ≈ 0.95, failures are 5% of courses, so WHO's ~5% died would need nearly every failure to die.
  - **In reality,** many deaths on treatment are early deaths among people with advanced disease whom the drugs would otherwise have cured.
  - **Course length** has no effect either.

**Hazards by track and phase include both as special cases:**

- **μ_c = 0** gives deaths only among people whose regimen would fail.
- **A cure-track hazard in the first phase only** gives early deaths among people the drugs would cure, before they act.
- **The data and calibration can choose between them.**

**How deaths then respond to a regimen:**

- **Higher efficacy:** fewer courses on the fail track, so fewer deaths from μ_f.
- **A shorter course:**
  - **Deaths on treatment fall,** through less time exposed.
  - **For someone whose course will fail,** a shorter course mainly moves their death risk to after the course, as symptomatic TB at `sym_dead`. A shorter regimen then saves lives mainly through completion: fewer people lost, so more cures, and earlier retreatment after failure.
  - **Beyond that,** unless its hazards are assumed to differ. We think this is the right behaviour for the model to show, not a flaw.
- **Front-loaded cure-track deaths** (phase 1 only) are the same for any regimen with the same first phase. Neither efficacy nor duration reduces them, which is plausible for deaths in people who are already very ill when treatment starts.

### Loss to follow-up: a hazard, with partial benefit by phase

**LTFU is a hazard h, the same on both tracks,** not a share drawn at the start.

- **The share lost depends on the course length,** about 1 − e^(−hD), so a shorter regimen loses fewer people at the same hazard, and those lost have completed more of it.
- **"Shorter regimens have less LTFU"** can then be a lower h for that regimen, as an input.
- **A share drawn at the start** would ignore D.

**Partial benefit:** someone on the cure track lost in phase k is cured with probability q_k = q(f_k), where f_k is the fraction of the course completed by the middle of phase k.

- **This experiment's rule:** q(f) = f, a placeholder, with q = 0 (no partial benefit) as the other bound.
- **On the fail track,** loss to follow-up never cures.

### Back to the state treatment started from

As CONVENTIONS.md 9 sets out, each state treatment can start from has its own set of treatment states, so nobody needs a per-person memory of where they came from.

- **Each set:** M cure-track phases and M fail-track phases (2M states).
- **Every exit that isn't a cure or a death** returns to that starting state: failing at the end of the course, and loss to follow-up without cure.
- **Relapse after cure** goes to symptomatic TB, as now.
- **The sets are generated from a list of starting states** (`tx_from`). Today the list is symptomatic TB only. Arm B tests two more, so the return paths are checked before anything relies on them.

**A course's outcomes depend on whether the starting state has TB disease** (user, 2026-11-23):

- **From a state with TB disease** (symptomatic, asymptomatic or non-infectious TB): as above, two tracks, death and loss to follow-up hazards, and cure to `treated` or `relapsing`.
- **From a state without TB disease** (susceptible, latent infection, cleared, recovered, treated), e.g. a false-positive diagnosis: one track, with no cure.
  - **During the course:** no TB death hazard, only background deaths, and loss to follow-up at h.
  - **Afterwards:** everyone returns to their starting state, with no relapse and no change in reinfection risk. Sending them to `treated` would give them protection against reinfection (π) they haven't earned, and `relapsing` a relapse risk they can't have.
  - **Counting:** they still count as notifications and as on treatment, as false positives would.
- **Treating latent infection on purpose** (preventive therapy) is a separate product with its own effect, not covered by this rule.

### Unchanged

- **Who starts treatment:** symptomatic TB at θ(t) = −ln(1 − cov(t)) per year, where cov is the annual coverage.
- **On treatment:** nobody transmits or is reinfected.
- **After cure:** `treated`, or `relapsing` with share `p_relapse`, then symptomatic TB at 1/1.5 per year.
- **Notifications:** every start.
- **Incidence:** relapses count, and first entries are counted as CONVENTIONS.md 5.

## Parameters

| Symbol (code) | Definition | Source in this experiment | Source later (after reparameterising) |
|---|---|---|---|
| cov(t) (`tx_rate` = −ln(1 − cov)) | Annual treatment coverage of symptomatic TB | Calibrated (exp 32's posterior) | Calibrated |
| D (`tx_dur`) | Course length | Regimen definition: 0.5 years | Regimen definition, per regimen |
| w (`tx_phase_dur`) | Phase width; M = D/w phases per course | Modelling choice: 2 months (M = 3); checked against finer phases | Modelling choice, the same for every regimen |
| c(t) | Cure share among people who complete the regimen (the cure track) | Derived: set so that overall cure among people starting treatment equals the calibrated efficacy s(t) | **Calibrated**, prior from trials (about 0.9–0.97 for drug-susceptible TB) |
| s(t) (`efficacy_1995`, `efficacy_2024`) | Overall cure among people starting treatment | Calibrated (exp 32); held fixed by the carve-out | Output, checked against WHO's reported treatment success |
| d(t) | Share of starts who die on treatment | WHO outcome data (India, by cohort year; target of the solve below) | WHO outcome data, as a target or to fix the hazards |
| ℓ(t) | Share of starts lost to follow-up | WHO outcome data (as d) | WHO outcome data, as d |
| μ_f,k(t) | Death hazard on the fail track, phase k | Solved so the died share equals d(t); the same in every phase | From d(t) given c, or calibrated |
| μ_c,k(t) | Death hazard on the cure track, phase k | Arm assumption: 0 (baseline), or phase 1 only, carrying a share ω of deaths on treatment | Literature on the timing and causes of deaths on treatment; or calibrated |
| ω | Share of deaths on treatment on the cure track | Arm assumption: 0, or 0.5 | Follows from μ_c |
| h(t) (`tx_ltfu`) | Loss-to-follow-up hazard, both tracks, every phase | Solved so the lost share equals ℓ(t) | From ℓ(t), or calibrated; per regimen |
| q_k (`tx_ltfu_cure`) | Chance that someone on the cure track lost in phase k is cured | Placeholder q(f) = f, and q = 0 | Literature (outcomes after incomplete treatment, e.g. relapse by doses taken) |
| `p_relapse` | Share of cures that relapse | Calibrated (exp 32) | Calibrated |
| `dur_relapse_tx` | Mean time from cure to relapse | Fixed, 1.5 years (tb_LAI_TPT) | Fixed, or literature |
| `sym_dead`(t) | TB death rate in symptomatic TB | Calibrated (exp 32); now off treatment only | Calibrated |
| `tx_from` | States treatment can start from; a state without TB disease gets one track with no cure | Symptomatic TB (arm B: also asymptomatic TB and susceptible, a test) | Design; screening or false positives would add states |

**WHO's outcome data.** The source is the treatment outcomes CSV, a public download from WHO's TB data portal. It gives India's cohort outcomes from the mid-1990s: success, failure, died, lost to follow-up, and not evaluated.

- **Which cohorts:** new smear-positive cohorts before 2012, and new and relapse cohorts after. The definitions change over time, which is reported.
- **Public sector only, in effect:** these are notified, mostly public-sector cohorts. Private-sector outcomes, often worse, were largely not reported, especially in earlier years. They are used as published here, labelled as public-sector outcomes. That's enough for this experiment, because the carve-out holds overall cure fixed, so they only decide how failures split into deaths, losses and true failures. Adjusting them for the private sector belongs with the reparameterisation (see "Not covered").
- **"Not evaluated"** is left out of the denominator (user, 2026-11-23); its share is reported.
- **Background deaths:** WHO's "died" counts deaths from any cause on treatment. So the model's deaths on treatment from TB and from background together are matched to d(t), not TB deaths alone. Background deaths are about 0.5–1% of courses at typical treated ages, against a died share of about 4–5%.
- **Smoothing:** d(t) and ℓ(t) are smoothed onto tb_LAI_TPT's knots (1995, 2010, 2024; linear between, constant after), as coverage and efficacy are.

**The carve-out (this experiment).** Exp 32's calibrated posterior is kept, and the new outcomes are carved out of its failure share, so overall cure among people starting treatment stays at each set's calibrated s(t). At each knot, three conditions fix three quantities:

- **The fixed outcome shares:** overall cure = s(t), the died share = d(t) (deaths from TB and background together), and the lost share = ℓ(t). These are solved numerically for c, μ_f and h, given q and ω, with the phases' competing exits computed exactly.
- **Background deaths are a competing exit in the solve,** not subtracted afterwards: someone who dies of other causes can't then die of TB, be lost or be cured.
- **Their rate is the age-weighted background death rate of the people starting treatment** in that year, taken from the model's run (people with symptomatic TB are mostly adults, and India's mortality at those ages falls over time). It is not the population's average.
- **c > 1** can happen where s is high and d and ℓ are not small. It is capped at 1 for the runs, and the draws where it happens are counted and reported.

**Later, the reparameterisation** (agreed, 2026-11-22; not in this experiment):

- **Calibrated:** c, with a prior from trials. It stays at or below 1 by construction.
- **From data or literature:** the hazards and q.
- **Overall success becomes an output.**
- **Likely to expose a tension:** the calibration currently uses low efficacy as a lever, and overall success would then be about 0.85.

## Question

1. **Does the structure run identically in both solvers?** v4's ODE (`TBTreat`) against v4's agents, and `fastode_age` against v4's ODE. The structure here is the two tracks, phases, death and LTFU hazards, partial benefit, and a second starting state.
2. **How much does it change the calibrated model?** At exp 32's posterior, with overall cure held at calibrated efficacy:
   - How do the eight run-based targets and the log-likelihood move?
   - How often would c exceed 1?
   - Where do the deaths fall: on which track, and how the fail-track hazard compares with `sym_dead`?

## Design

**Stage 0, data:** `data.py` downloads WHO's outcome CSV into `data/` and builds d(t) and ℓ(t) for India. Outputs: `outcomes.csv`, and `outcomes.png` with each year's raw shares and the smoothed curves.

**Stage 1, code:**

- **`tbode4.TBTreat`:** new options.
  - **Structure:** `tx_from`, `tx_phase_dur`, and pathways: a start of treatment can branch into several pathways, each with its own set of treatment states, outcomes and hazards (one pathway here). Stratifying by sector later is then a setting, not a rewrite.
  - **Hazards and benefit:** `tx_death_fail` and `tx_death_cure` (by phase; numbers or `ss.timeseries`), `tx_ltfu`, `tx_ltfu_cure` (by phase).
  - **With the defaults** (`tx_from=['symptomatic']`, one pathway, one phase of the whole course, no LTFU, `tx_death_fail` = `sym_dead`, `tx_death_cure` = 0), it is the current model. A regression test in `v4/tests/test_treat.py` checks it reproduces exp 34's ODE run step for step.
  - **The carve-out solve:** in `tbode4` (a function of s, d, ℓ, q, ω and D), with a test that a cohort run reproduces the three shares.
- **`tbode4.fastode_age`:** the same structure, tested against v4's ODE in `v4/tests/test_fastode_age.py` before any use for calibration (CONVENTIONS.md 8). `tbode4.fastode` (no ages) is no longer used for calibration and is left as it is.
- **Phases at weekly steps:** v4's ODE runs a fixed duration as whole steps (nearest). So three 2-month phases at weekly steps are 9 steps each, 27 weeks (0.52 years) in total, against 26 weeks for one 6-month phase. Agents run the exact duration.
  - **Measured in stage 2.**
  - **If it matters,** the phases are given lengths in whole steps that add up to the course.
- **The phase width:** with an exponential loss to follow-up, the expected partial benefit for a smooth q(f) has a closed form. It is compared at w = 6, 3, 2 and 1 months (M = 1, 2, 3, 6) against that continuous limit, at the solved h (`phases.log`); no simulation is needed. Stage 3 also compares the targets at w = 2 and 1 month.

**Stage 2, agents against the ODE:**

- **Setup:** exp 34's, with India's demography from 1950, each set's settled start by single year of age and sex, and 100,000 agents at the start.
- **Sets:** 3 of exp 34's sets, spread along the ridge (lowest, middle and highest `k_latent`), with their carve-out schedules (c capped at 1 if needed).
- **Arm A:** the design with every exit in use, to exercise the code: ω = 0.5 (cure-track deaths in phase 1), q(f) = f, treatment from symptomatic TB.
- **Arm B (the return paths):** as A, but treatment also starts from two more states at placeholder rates, each with its own set of treatment states. Not a model of screening or of false positives; it checks the return paths.
  - **Asymptomatic TB,** at a tenth of symptomatic TB's rate: a state with TB disease.
  - **Susceptible,** at a small rate (e.g. 0.001 per year): a state without TB disease, checking that nobody is cured, dies of TB on treatment or relapses, and that everyone returns to susceptible.
- **Runs:** v4's ODE at quarter- and eighth-year age groups (extrapolated to exact aging), and agents, 20 reps each, at monthly and weekly steps.
- **Measured:**
  - **As in exp 34:** agents ÷ ODE by window for prevalence (all ages and 15+), for incidence, notifications, the number on treatment and TB mortality per prevalent case, and for prevalence by age group.
  - **New:** deaths on treatment by track, losses to follow-up and partial cures, per prevalent case.
- **Extended pitfall checks** (`tbode4.TxChecks`):
  - everyone leaving treatment without cure or death returns to the state they started from;
  - nobody treated from a state without TB disease is cured, dies of TB on treatment, or enters `treated` or `relapsing`;
  - each phase lasts w (agents);
  - the shares who die and who are lost, against d(t) and ℓ(t), and the deaths' split by track against ω;
  - partial cures only from the cure track, at q_k by phase;
  - no TB death at `sym_dead` on treatment.
- **What would answer question 1,** as exps 27 and 34:
  - at weekly steps, each set's prevalence within about 2 SE of the ODE in each window;
  - ratios per prevalent case within about 0.5% on average;
  - every pitfall check passing;
  - `fastode_age` within 0.2% of v4's ODE.

**Stage 3, the calibrated model** (`fastode_age` only, monthly steps, as exp 32):

- **Draws:** exp 32's 200 posterior draws from each seed (`real_draws.csv`, `real_seed2_draws.csv`).
- **Each draw is run five ways:**
  - the current model;
  - the carve-out baseline: ω = 0 and q(f) = f;
  - as the baseline, with q = 0;
  - as the baseline, with front-loaded cure-track deaths, ω = 0.5;
  - as the baseline, with 1-month phases (the phase width check).
- **Measured:**
  - each target, new ÷ current (median and range over draws), and the change in log-likelihood;
  - how many draws have c > 1, and by how much;
  - the solved hazards: μ_f against `sym_dead`, h, and μ_c,1;
  - the shares of starts cured, dying (by track) and lost, against WHO's.
- **What would answer question 2:**
  - **Targets move by less than their uncertainty and c ≤ 1 almost everywhere:** the posterior stands with the new structure until the reparameterisation.
  - **Otherwise:** a recalibration is needed with the new structure, best done after the reparameterisation.

## Not covered

- **The reparameterisation** to clinical efficacy, and the recalibration (next, after this).
- **The private sector** (next, with the reparameterisation, when the outcome data start to matter):
  - **First, adjustment factors:** overall death and loss shares as a mix of public outcomes (WHO's data) and worse private outcomes (literature), weighted by the private sector's share of treatment over time, with priors.
  - **Stratification by sector** (the pathways above) only if sectors need to differ in more than outcomes. Notification is the likely case: much private-sector treatment wasn't notified, the model ties notifications to every start through `over_detection`, and notifications are among the targets the calibration only just meets.
- **Other regimens** (shorter courses, different hazards): the structure allows them, but nothing is compared here.
- **Lower coverage for children,** and death or LTFU varying by age.
- **Treatment from states other than symptomatic TB** as a model; arm B is a test only.
- **Transmission on treatment.**

## Files (planned)

- **`data.py`:** WHO's outcomes → `data/`, `outcomes.csv`, `outcomes.png`.
- **`run.py`:**
  - `check`: stage 1's regression and copy checks → `check.log`;
  - `phases`: the phase-width check → `phases.log`;
  - `ode`, `agents`, `summary`: stage 2 → `ode.obj`, `runs.obj`, `by_set.csv`, `pooled.csv`, `pitfall_checks.csv`, `summary.log`, figures;
  - `targets`: stage 3 → `targets.csv`, `targets.log`, figures.
- **Changed in `tbode4`:** `treat.py` (`TBTreat`, `TxChecks`, the carve-out solve), `fastode_age.py`; tests in `v4/tests/`.
- **Reuses:** exp 34's run setup (sets, starting tables, `AgeCounts`, `ode_incidence`), exp 32's `calib.py` (parameter mapping, targets).
- **`SUMMARY.md`:** findings, written after running.
