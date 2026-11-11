# Starter: running TBsim as agents

`run_tbsim.py` runs TBsim v4 as an agent-based model and plots active TB prevalence, incidence and TB mortality over time:

- **Agents:** 50,000, in 4 runs with different random seeds.
- **Period:** 1950–2000, at monthly steps.

It builds in the choices worked out in this repo (listed below), so you can start from a setup that is consistent with TBsim's ODE version and with the calibration.

## Running it

From the repo root, set up the pinned environment once, then run the script:

```sh
cd v4 && uv sync && cd ..
v4/.venv/bin/python starter/run_tbsim.py
```

If another Python environment is active, run `unset VIRTUAL_ENV` first. A run takes a few seconds per seed once Starsim has compiled its kernels (`v4/.venv/bin/starsim compile`, about 2 minutes, once).

**Outputs**, saved next to the script:

- `tbsim_burn_in.csv`: yearly averages for each run, per 100,000 people;
- `tbsim_burn_in.obj`: the same, with the parameters and starting state (load with `sciris.load`);
- `tbsim_burn_in.png`: the figure. Each run is shown, with their mean, and the dashed line is the ODE's equilibrium.

## What's built in

- **TBsim v4 on Starsim v4,** pinned to exact commits (`v4/pyproject.toml`), so results are reproducible. Later versions may behave differently.
- **Parameters:** one draw from a posterior calibrated to India's prevalence survey (prevalence, share asymptomatic, prevalence ÷ incidence) and to data on progression after infection (`experiments/20_ode_calibration`). It has prevalence 312 per 100,000, against a target of 310. Other draws are in `experiments/21_calibrated_comparison/sets.obj`.
- **Latent infection** has two stages, recent and remote. Most progression to TB happens soon after infection, and remote infection progresses at a small fraction (`phi`) of the recent rate.
- **Reinfection:** anyone infected before can be reinfected, at a relative risk `pi` compared with someone never infected. Reinfection of remote latent infection makes it recent again.
- **Homogeneous mixing:** everyone mixes with everyone, with no contact network. This matches TBsim's ODE directly and runs about twice as fast as a random network.
  - **If you need a random contact network,** calibrate the per-contact transmission parameters on Starsim's ODE of the network, which is the same model as agents on it (`CONVENTIONS.md`, section 1). Dividing `beta` by the number of contacts is not enough.
- **Starting at equilibrium:** agents start in each state in the proportions of the ODE's endemic equilibrium, computed from the parameters. If you change the parameters, the script recomputes the starting state.
  - **The average needs no burn-in.** Averaged over runs, prevalence stays within about 2% of the equilibrium from the first year (240 runs, `experiments/21_calibrated_comparison`).
  - **The spread between runs takes time to build up.** Every run starts from exactly the same proportions, so the runs are nearly identical at first. They take about 20–25 years to reach their full variability (about ±14% in a single year's prevalence). Starting in 1950 leaves enough time before 2000. If you want uncertainty bands for an earlier period, start earlier.
  - **This works only while the model has no memory** beyond each person's current state. Age structure, fixed treatment or relapse delays, or parameters that change over time would need a realistic starting population, or a run through history, instead.
- **Incidence** is counted as the first entry to asymptomatic TB in each episode, from the transitions in each step (`CONVENTIONS.md`, section 5). It's not counted from the number of people with TB at the end of each step, which misses short episodes. People who move between asymptomatic and non-infectious TB within one episode are counted once.
- **TB mortality** is TB deaths per 100,000 per year. Births replace every death, so the population stays at 50,000.

## Limitations

- **No demography:** the population is constant, with no age structure, no age-specific death rates and no population growth. Everyone dies at the same background rate, and births replace every death.
- **No treatment:** with no diagnosis or treatment, burden stays at its endemic level. The model can't show the declines over time that treatment scale-up would bring.
  - **TB mortality is too high relative to incidence** (about 0.25 deaths per new case, against about 0.11 in the data), because nobody is cured. The calibration left this target out until treatment is added.
- **No household network,** and none of the other structure in tb_LAI_TPT's model: no lower transmission from children, no time-varying parameters, no drug resistance.
- **One parameter draw:** the results don't reflect uncertainty in the calibration. Rerun with other draws to see it.
- **Noise:** at 50,000 agents, a run's prevalence averaged over 1950–1999 varies by about 9% from run to run (about 14% in a single year). The mean of 4 runs is uncertain by about ±5% (one standard error). Use more runs, or more agents, for precise comparisons.
