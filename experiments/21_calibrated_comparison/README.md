# 21: Do v4's agents match its ODE across the calibrated posterior?

**Status:** design confirmed by the user (2026-11-07); run (see SUMMARY.md).
**Track:** v4 (`v4/.venv`; tbsim `origin/v4` @ `dcfb522`, starsim `origin/v4` @ `80c81f53`).

## Background

- **What earlier experiments found:** v4's agents and ODE agree within about 1%, with homogeneous mixing. Exps 13 and 17 showed this for prevalence; exp 19 showed it for incidence counted as the first entry to asymptomatic TB per episode. But those experiments used two hand-picked parameter sets: set A (tbsim's defaults) and set B (tb_LAI_TPT's best fit, with the two-state latent model at φ = 0.01).
- **What exp 20 added:** a calibrated posterior for the ODE without treatment. It fits five of tb_LAI_TPT's targets: prevalence, the share asymptomatic, prevalence ÷ incidence, and the two progression targets. Its parameters range widely: `k_latent` 2.3–11.6, φ 0.008–0.095, π 0.21–0.66, `beta` 10.7–34.1. It also pushes the disease cascade towards keeping people asymptomatic.
- **What exp 20 doesn't show:** whether agents still match the ODE across that range, including the posterior's corners (fast or slow front-loading, large φ, high π).

## Question

At parameter sets drawn from exp 20's posterior, do v4's agents reproduce v4's ODE? In particular:

1. **The targets:** prevalence, the share asymptomatic, and prevalence ÷ incidence. If agents matched the ODE, they would fit the targets as well as the ODE does.
2. **The other outputs:** latent infection, TB mortality, and mortality ÷ incidence. Mortality ÷ incidence isn't fitted.

## Design

**Parameter sets:** 12 draws from exp 20's posterior (`real.obj`), by weighted resampling without replacement, so each is a distinct point and together they represent the posterior.

- **Each set counts equally,** as resampled draws.
- **The seed is fixed** and the draws are saved, so the sets can be reused, e.g. for the treatment work.
- **Why 12:** it matches tb_LAI_TPT's face-validity sets, and is enough to cover the posterior's range. The table of drawn sets is checked against the posterior's 95% intervals before running.

**Model:** v4's `tbsim.TB` with the two-state latent model, homogeneous mixing, no treatment, births replacing deaths (`mu` = 1/70), every reinfection risk at π. These are exp 20's settings, mapped to tbsim's parameters with `calib.to_pars`.

**Runs for each set:**

- **The ODE:** v4's own ODE (the reference), started at the standalone ODE's equilibrium (`tbode4.fastode`), 1900–2000. It should stay at the equilibrium, which checks the mapping. First-entry incidence comes from the standalone ODE, because v4's ODE doesn't count it (exp 19).
- **Agents:** 50,000, started at the same equilibrium, 1900–2000, at monthly and weekly `dt`, 20 reps each (new seeds). Incidence is counted with exp 19's `Incidence` analyzer.
- **Size:** 12 sets × 2 steps × 20 reps = 480 agent runs, about 20 minutes on 4 cores.

**Measured:** agents ÷ ODE, the mean over 1950–2000, for each set and step, with mean ± SE across reps. Then averaged over the 12 sets, with SE.

- **Outputs:** active prevalence; the share asymptomatic; incidence (first entry) per prevalent case and prevalence ÷ incidence; TB mortality and mortality ÷ incidence; latent infection, recent and remote; the cleared share.
- **Expected reading high:** recent latent infection in end-of-step snapshots at monthly `dt`, which reads ×1.2–1.3 when recent infections are short-lived (exp 17). It's reported, but not counted as a mismatch.
- **Against the targets:** the agents' prevalence, share asymptomatic and prevalence ÷ incidence, against each target's 95% interval (exp 20's README).

**Figures:**

- agents ÷ ODE by set, for each output, at both steps, with each set's posterior parameters shown alongside;
- the agents' target outputs against the target intervals, next to the ODE's.

The raw outputs behind each figure (per-rep yearly means, the ODE runs, the drawn sets) are saved with `sc.save`.

## What would answer the question

- **Agents match the ODE across the posterior** if every set's agents ÷ ODE is within about 2 SE of 1, and the average over sets is within about 1%. 20 reps give an SE of about 1–2% of prevalence for one set, and the average over 12 sets about 0.5%.
- **If some sets don't match,** check what sets them apart (their position in the posterior, e.g. a high `k_latent`, high φ or low `beta`) before calling it a mismatch. Rerun those sets with more reps, and at a daily step, to separate step effects from noise.

## Not covered

- Treatment, demography, age structure and networks, which aren't in exp 20's model.
- The posterior's tails beyond what 12 draws reach.

## Files

- `run.py`: `sets` draws the 12 sets (seed 21) → `sets.obj` (each set's parameters, posterior index and weight, standalone-ODE equilibrium and outputs), `sets.log` (the drawn sets against the posterior's 95% intervals). With no argument it runs everything → `runs.obj` (v4 ODE runs and agents' per-rep yearly means), `by_set.csv`, `pooled.csv`, `agents_vs_ode.png`, `targets.png`, `run.log`, `versions.json`. `summary` redraws the summaries and figures.
- Uses exp 20's `calib.py` (parameter mapping, targets) and exp 19's `Incidence` analyzer.
- `SUMMARY.md`: findings, written after running.
