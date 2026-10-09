# 03 Baseline calibrated to a prevalence target, started at equilibrium

## Question

With `beta` tuned so the ODE equilibrium hits the burden target, do the ABM and ODE agree at steady state, at daily, weekly and monthly ABM `dt`? Are the implied incidence, mortality, % asymptomatic and progression reasonable?

## Design

- **Target:** active (ASY+SYM) prevalence of 500/100k, all ages, with an acceptable range of 200–800. This is the only calibration target.
- **Natural history:** tbsim defaults (`TBODE` defaults), with `k_asy = k_non = 0`.
  - tb_LAI_TPT's calibrated rates were not used. Without time-since-infection progression they put 42% of a latent cohort into active disease within 10 years. They are deferred to the time-since-infection experiment.
  - Background mortality μ = 1/70 per year. Births = μN. No treatment or other interventions.
- **ODE:** `beta` found by root-finding on the equilibrium prevalence (`tbode.calibrate_beta`). The `beta` values for 200 and 800/100k are also reported.
- **ABM:** `tbode.make_abm` at the calibrated `beta`, with states assigned from the ODE equilibrium (`tbode.init_from_ode`), so no burn-in is needed. Runs cover 100 years with 50,000 agents and 10 seeds each, at `dt` = 1 day, 1 week and 1 month. The daily arm was added in a rerun (2026-10-09) to check that the drift disappears as `dt` shrinks. The weekly and monthly arms were rerun with the same seeds. Any drift away from the equilibrium measures the inconsistency.
- **Outputs:**
  - prevalence;
  - incidence: first ASY entry per episode, as in tb_LAI_TPT's `IncidenceTracker` with `sym_only=False`;
  - TB mortality;
  - % asymptomatic (ASY/(ASY+SYM));
  - state proportions;
  - progression: ProgressionY1 and ProgressionRatio in a no-reinfection latent cohort, from the ODE (exact) and from the ABM at each `dt`. Targets are 0.0523 and 0.663.
- **Event rates** (incidence, mortality) are compared as calendar-year averages.
- **Extra daily reps (`run_more_reps.py`, 2026-10-09):** with 10 reps, the daily mean (469 ± 21/100k SE) could not be distinguished from the ODE's 500. 30 more daily reps (seeds 10–39) bring the SE down to about 2%. They are compared using the mean and SE across reps, as well as the median and the 10–90% range.
- **Saved data:** `outputs.obj` holds the per-seed ABM time series, the ODE run and the cohort results behind the figure.

## Versions

tbsim main @ `b92d7e8` (v0.10.3), starsim @ `3521adb4` (v3.6.1).
