# 06 Profiling: how much slower is `TBExact` than `tbsim.TB`?

## Question

`tbode.TBExact` is a version of `tbsim.TB` that removes tbsim's step-ordering bias. It draws each agent's end-of-step state from the exact one-step probabilities of the continuous-time model, `expm(Q·dt)`, instead of applying transitions one state at a time. In exp 05 it matched the ODE at daily, weekly and monthly `dt`. How much longer does it take to run than `tbsim.TB`, at each `dt`? And which parts of the step account for the time?

Run times from exps 03 and 05 suggest a modest difference: roughly 440 s vs 510 s per 100-year daily run. But those runs were 4 at a time in parallel and mixed `dt`s, so they don't give a clean comparison.

## Design

- **Setup:** same as exps 03 and 05. tbsim default natural history, ODE `beta` = 9.917 (500/100k), 50,000 agents, the ABM started at the ODE equilibrium (`init_from_ode`), with `ss.RandomNet(n_contacts=10, dur=0)`, constant births and deaths, and the `StateCounts` analyzer.
- **Arms:** `tbsim.TB` and `TBExact`, at `dt` = 1 day, 1 week and 1 month.
- **Length:** 10 simulated years per run. Time per step doesn't depend on run length, so results scale to 100-year runs.
- **Reps:** 3 per arm, each with a different seed, plus one warm-up run that is discarded.
- **Timing conditions:** all runs in a single process, one after another (no parallel runs, to avoid competing for cores), with the arms interleaved so slow drift in machine load affects them equally.
- **Measured:**
  - wall time for `sim.run()`, and separately for `sim.init()` plus `init_from_ode`;
  - time per step;
  - the per-module breakdown from `sim.run(profile=True)` (`sim.loop.cpu_df`), e.g. the TB step vs the network vs the analyzer.
- **Outputs:**
  - `timings.csv`: one row per run;
  - `cpu_breakdown.csv`: time per loop entry, per run;
  - `summary.csv`: per arm, the median time per step and per 100-year run, and the `TBExact` ÷ `tbsim.TB` ratio;
  - `timings.png`, drawn from `timings.csv`.
- **Caveats:** timings are specific to this laptop and to 50k agents. Other load on the machine adds noise, which the 3 reps help gauge.

## Versions

tbsim main @ `b92d7e8` (v0.10.3), starsim @ `3521adb4` (v3.6.1).
