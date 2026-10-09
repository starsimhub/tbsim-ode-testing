# 06 Profiling: summary

## Background

`tbode.TBExact` is a version of `tbsim.TB` that removes tbsim's step-ordering bias. It draws each agent's end-of-step state from the exact one-step probabilities of the continuous-time model, `expm(Q·dt)`, instead of applying transitions one state at a time. In exp 05 it matched the ODE at daily, weekly and monthly `dt`, where `tbsim.TB` fell 5–33% short. This experiment measures what it costs in run time.

**Setup:** as in exps 03 and 05 (50k agents, started at the ODE equilibrium, `ss.RandomNet`, constant births and deaths, the `StateCounts` analyzer). Each run covers 10 simulated years, 3 reps per arm, all run one after another in a single process on this laptop.

## Result: `TBExact` costs nothing extra

| `dt` | `tbsim.TB` ms per step (range over reps) | `TBExact` ms per step (range over reps) | `TBExact` ÷ `tbsim.TB` | Implied 100-year run, either model |
|---|---|---|---|---|
| 1 day | 10.3 (10.0–10.6) | 10.4 (10.2–10.8) | 1.01 | ~6.3 min |
| 1 week | 11.2 (10.9–12.3) | 10.9 (10.8–11.3) | 0.98 | ~1 min |
| 1 month | 10.9 (10.5–11.2) | 10.8 (10.5–10.9) | 0.99 | ~13 s |

- **The ratios are 0.98–1.01, within the rep-to-rep spread,** so there's no measurable difference.
  - `TBExact` computes one 7×7 matrix exponential per step and draws one random number per infected agent.
  - `tbsim.TB` applies a sequence of state-by-state draws instead.
  - The two cost about the same.
- **Time per step is about 10.5 ms at every `dt`,** so run time is simply proportional to the number of steps. Daily is about 30× the cost of monthly, and weekly about 4.3×.
- **Initialisation is negligible:** 0.03–0.17 s, including `init_from_ode`.

**Where the time goes** (ms per step, median over reps; nearly identical for both models and all `dt`s):

| Loop entry | ms per step | Share |
|---|---|---|
| TB module step (`tb.step` / `tbexact.step`: transmission and transitions) | 4.5–4.6 | ~43% |
| People updates (`people.finish_step`) | 1.9–2.5 | ~20% |
| Random network (`randomnet.step`) | 2.0 | ~19% |
| TB results (`update_results`) | 0.65 | ~6% |
| `StateCounts` analyzer | 0.46 | ~4% |
| Births | 0.23–0.26 | ~2% |

The TB module step takes the same time in both models, so the transitions are a small part of it.

**Earlier run times weren't a clean comparison.** Run times from exps 03 and 05 (about 440 s vs 510 s per daily 100-year run) came from 4 runs at a time in parallel, with other `dt`s mixed in. Measured cleanly, a 100-year daily run takes about 6.3 minutes for either model.

## Implication

Removing the step-ordering bias with `TBExact` doesn't cost any run time. Run time is driven by `dt`: daily `dt` with `tbsim.TB` is 30× more expensive than monthly `TBExact`, and still 5% biased. This applies to the constant-rate natural history only. With time-since-infection progression, exact sampling needs each agent's hazard integrated over the step, which may cost more and should be profiled again then.

Outputs: `timings.csv` (one row per run), `cpu_breakdown.csv` (per loop entry, per run), `summary.csv`, and `timings.png` (drawn from `timings.csv`).
