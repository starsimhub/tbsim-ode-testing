# Starter: running TBsim as agents

`run_tbsim.py` runs TBsim v4 as an agent-based model and plots active TB prevalence, incidence and TB mortality over time:

- **Agents:** 50,000, in 4 runs with different random seeds.
- **Period:** 1900–2000, at monthly steps.

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
  - **If you need a random contact network,** the transmission parameters must be converted (`CONVENTIONS.md`, section 1), not just divided by the number of contacts.
- **Starting at equilibrium:** agents start in each state in the proportions of the ODE's endemic equilibrium, computed from the parameters. So no long burn-in is needed, and the lines stay flat around the dashed equilibrium. If you change the parameters, the script recomputes the starting state.
- **Incidence** is counted as the first entry to asymptomatic TB in each episode, from the transitions in each step (`CONVENTIONS.md`, section 5). It's not counted from the number of people with TB at the end of each step, which misses short episodes. People who move between asymptomatic and non-infectious TB within one episode are counted once.
- **TB mortality** is TB deaths per 100,000 per year. Births replace every death, so the population stays at 50,000.

## Limitations

- **No demography:** the population is constant, with no age structure, no age-specific death rates and no population growth. Everyone dies at the same background rate, and births replace every death.
- **No treatment:** with no diagnosis or treatment, burden stays at its endemic level. The model can't show the declines over time that treatment scale-up would bring.
  - **TB mortality is too high relative to incidence** (about 0.25 deaths per new case, against about 0.11 in the data), because nobody is cured. The calibration left this target out until treatment is added.
- **No household network,** and none of the other structure in tb_LAI_TPT's model: no lower transmission from children, no time-varying parameters, no drug resistance.
- **One parameter draw:** the results don't reflect uncertainty in the calibration. Rerun with other draws to see it.
- **Noise:** at 50,000 agents, one run's prevalence varies by about 9% from run to run. The mean of 4 runs is uncertain by about ±5% (one standard error). Use more runs, or more agents, for precise comparisons.
