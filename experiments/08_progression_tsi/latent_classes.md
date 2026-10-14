# Reference: Poisson-mixture latent classes

This note explains how the ODE (`tbode.TBODE`) and the corrected ABM (`tbode.TBExact`) represent tbsim's progression by time since infection, as used in this experiment. The method is exact, not an approximation: it gives the same distribution of progression times as tbsim's declining hazard.

## What tbsim does

In `tbsim.TB.progression_rates`, a latent agent (state INFECTION) at time τ since infection leaves INFECTION by progression at the rates

- `inf_non·e^(−k_non·τ)` to NON_INFECTIOUS
- `inf_asy·e^(−k_asy·τ)` to ASYMPTOMATIC

Here τ is in years, counted from `ti_infected`, so reinfection restarts the clock. `k_asy` and `k_non` are decay rates per year. This experiment uses `k_asy = k_non` = 5.846, tb_LAI_TPT's best fit, so the hazard falls by a factor of e every 1/`k` ≈ 2 months. Below, `k` stands for this common value.

Clearance (`inf_cle`) and background death (μ) don't depend on τ.

With `k_asy = k_non = k`, the total progression hazard is r₀·e^(−kτ), where r₀ = `inf_non` + `inf_asy`. A fraction f = `inf_non` / r₀ of progressions go to NON_INFECTIOUS, whatever τ is.

## Notation: Poisson(m; λ)

A Poisson distribution has a single parameter, λ, which is its mean. **Poisson(m; λ)** means the probability that a Poisson(λ) random variable equals m:

Poisson(m; λ) = e^(−λ) · λᵐ / m!,  for m = 0, 1, 2, …

m isn't a second parameter; it's the value where the probability is evaluated. In code this is `scipy.stats.poisson.pmf(m, λ)`.

Here λ = r₀ / `k`, and m is the index of a latent class.

## Why progression time is a mixture of exponentials

Ignoring clearance and death for now, the probability of not yet having progressed by time τ is

S(τ) = exp(−∫₀^τ r₀·e^(−ks) ds) = exp(−λ·(1 − e^(−kτ))),  with λ = r₀ / `k`.

Write this as e^(−λ) · exp(λ·e^(−kτ)) and expand the second factor as a power series:

S(τ) = Σₘ [e^(−λ) · λᵐ / m!] · e^(−m·k·τ) = Σₘ Poisson(m; λ) · e^(−m·k·τ).

Each term e^(−m·k·τ) is the survival function of an exponential distribution with rate m·`k`. So the progression time is exactly a mixture of exponentials, with Poisson weights. This is an identity: the only approximation in practice is truncating the sum (see below).

## The classes

INFECTION is split into latent classes m = 0, 1, 2, …:

- **Assignment.** At infection, an agent joins class m with probability Poisson(m; λ). Equivalently, draw a Poisson(λ) random number; that number is the agent's class. The class is fixed for as long as the agent stays in INFECTION.
- **Progression.** Class m progresses at the constant rate m·`k`: m·`k`·f to NON_INFECTIOUS and m·`k`·(1 − f) to ASYMPTOMATIC. Class 0 never progresses.
- **Clearance and death.** Every class clears at `inf_cle` and dies at μ, as in tbsim.

### Why the population hazard still declines

Each class has a constant rate. But among agents still latent at τ, the fast classes have left first, so the survivors become concentrated in low classes. The average progression hazard over survivors is

Σₘ wₘ · m·`k` · e^(−m·k·τ) / Σₘ wₘ · e^(−m·k·τ),  with wₘ = Poisson(m; λ).

With Poisson weights, the numerator equals `k`·λ·e^(−kτ) times the denominator. So the ratio is exactly `k`·λ·e^(−kτ) = r₀·e^(−kτ), which is tbsim's hazard.

### Why clearance and death don't break it

Clearance and death act at the same rate in every class. So they remove survivors without changing the mix of classes among those still latent, and the result above holds with these competing risks included.

## Which models use the classes

- **ODE (`TBODE`):** yes. INFECTION is split into compartments INFECTION_0, INFECTION_1, …, and each new infection is spread across them by the Poisson weights.
- **`TBExact`:** yes. Each agent gets a class, drawn from the same weights the first time it is stepped after infection. Its end-of-step state then comes from that class's constant-rate `expm(Q·dt)`. `TBExact` needs the classes because exact sampling with `expm(Q·dt)` requires rates that are constant over the step: within a class they are, but tbsim's declining hazard isn't.
- **Standard `tbsim.TB`:** no. It applies the hazard r₀·e^(−kτ) directly, using each agent's actual time since infection.

The ODE and `TBExact` share the same class construction, so `TBExact` matching the ODE doesn't test the mixture itself. Checks that don't use classes do: the comparison with the directly integrated hazard, and daily `tbsim.TB` cohorts matching the exact progression curve (see "Where it is in the code" below).

## Why `TBExact` needs the classes

**What `expm(Q·dt)` does.** Q is the generator matrix: entry (i, j) is the rate from state i to state j. If the rates stay fixed for the whole step, `expm(Q·dt)` gives the exact probability of ending the step in each state, from each starting state, over every path within the step (including several moves). That depends on Q being the same throughout the step.

**Why tbsim's hazard breaks this.** A latent agent's progression rate, r₀·e^(−kτ), falls during the step: at monthly `dt`, by a factor of e^(−5.846/12) = 0.61. So there's no single Q for the step. Two further problems:

- **Q differs between agents.** Each agent has its own τ, so each would need its own calculation every step, rather than one matrix shared by everyone.
- **Averaging the rate over the step isn't exact.** `expm` of the step-averaged rate is still an approximation, because when an agent progresses within the step affects where it can get to next. An early progressor has time to move on, e.g. from ASYMPTOMATIC to SYMPTOMATIC.

Exact per-agent probabilities would need a numerical integral over the progression time for each latent agent in each step. tbsim instead holds the start-of-step rate for the whole step. Because the rate is actually falling, this overstates progression: at monthly `dt`, the expected number of progressions per agent is about 1.26× the exact value. That is the `tbsim.TB` bias this experiment found.

**How the classes solve it.** Within class m, the progression rate m·`k` is constant. So each class has a fixed Q, and `expm(Q·dt)` is exact for it. There are only 7–8 classes, so `TBExact` computes 7–8 matrix exponentials per step and shares them across all agents, whatever their τ.

**The class takes the place of τ.** A new infection's class is Poisson(λ). Because high classes progress first, an agent still latent at τ has class distributed as Poisson(λ·e^(−kτ)). This is tbsim's falling hazard, expressed as a shift in which classes remain. Clearance and death don't change this, because they act the same on every class. Within a class the process is memoryless: what happens next depends on the class and the current state, not on how long the agent has been latent. That is what `expm` needs. It is also why `init_from_ode` gives each latent agent a class when starting at equilibrium.

## Probability of progression

Within class m, rates are constant, so progression is an ordinary competing-risks calculation:

P(progress | m) = m·`k` / (m·`k` + `inf_cle` + μ).

Overall, P(progress) = Σₘ Poisson(m; λ) · P(progress | m). Here "progress" means leaving INFECTION for NON_INFECTIOUS or ASYMPTOMATIC rather than clearing or dying. It isn't the same as ever reaching ASYMPTOMATIC, which ProgressionY1 and ProgressionRatio measure.

Without clearance and death, P(progress) would be 1 − e^(−λ): all the weight outside class 0. So Poisson(0; λ) = e^(−λ) is the probability that the hazard dies away before progression ever happens.

## Values in this experiment

Both sets use `k` = 5.846/yr and μ = 1/70/yr. Set A uses tbsim's default natural history; set B uses tb_LAI_TPT's best fit (see the [README](README.md)).

| | Set A | Set B |
|---|---|---|
| `inf_non`, `inf_asy` (/yr) | 0.16, 0.06 | 0.410, 0.274 |
| r₀ (/yr) | 0.22 | 0.684 |
| λ = r₀ / `k` | 0.0376 | 0.117 |
| f (share of progressions to NON_INFECTIOUS) | 0.73 | 0.60 |
| `inf_cle` (/yr) | 1.9 | 0.537 |
| P(progress), mixture | 0.0279 | 0.1012 |
| P(progress), by integrating tbsim's hazard directly | 0.0279 | 0.1012 |
| P(progress) without clearance or death, 1 − e^(−λ) | 0.0369 | 0.110 |
| Classes used (cutoff 1e-10) | 7 (m = 0–6) | 8 (m = 0–7) |
| Classes with weight above 1e-6 | 4 (m = 0–3) | 5 (m = 0–4) |

The mixture and the direct integral agree to 12 significant figures in both sets.

Per class, set B:

| Class m | Assignment probability, Poisson(m; 0.117) | Progression rate m·`k` (/yr) | → NON_INFECTIOUS (/yr) | → ASYMPTOMATIC (/yr) | P(progress \| m) | Mean time in INFECTION |
|---|---|---|---|---|---|---|
| 0 | 0.890 | 0 | 0 | 0 | 0 | 22 months |
| 1 | 0.104 | 5.85 | 3.50 | 2.34 | 0.914 | 1.9 months |
| 2 | 0.0061 | 11.7 | 7.01 | 4.68 | 0.955 | 1.0 month |
| 3 | 0.00024 | 17.5 | 10.5 | 7.02 | 0.970 | 0.7 months |

Per class, set A:

| Class m | Assignment probability, Poisson(m; 0.0376) | Progression rate m·`k` (/yr) | → NON_INFECTIOUS (/yr) | → ASYMPTOMATIC (/yr) | P(progress \| m) | Mean time in INFECTION |
|---|---|---|---|---|---|---|
| 0 | 0.963 | 0 | 0 | 0 | 0 | 6.3 months |
| 1 | 0.036 | 5.85 | 4.25 | 1.59 | 0.753 | 1.5 months |
| 2 | 0.00068 | 11.7 | 8.50 | 3.19 | 0.859 | 0.9 months |
| 3 | 8.6e-6 | 17.5 | 12.8 | 4.78 | 0.902 | 0.6 months |

## Truncation: how many classes

The sum over m is infinite, so it is cut off at the first class beyond which the remaining Poisson weight is below `latent_tol` (default 1e-10). That leftover weight is added to the last class. This gives 7 classes for set A and 8 for set B; both counts include class 0. Set B needs one more class because its λ is larger.

The count depends only on λ = r₀ / `k`, which has no time units. So `TBExact`, which computes the classes from per-step rates, gets the same classes as the ODE at every `dt`.

A looser cutoff of 1e-6 would give 4 and 5 classes, which also matches tbsim's cohort progression to about 1e-6.

## Where it is in the code

- `tbode.latent_classes` ([ode.py](../../tbode/ode.py)) returns the class weights and each class's rates to NON_INFECTIOUS and ASYMPTOMATIC.
- **ODE.** `TBODE` splits INFECTION into compartments INFECTION_0, INFECTION_1, … when `k_asy = k_non` > 0. New infections are split across them by the weights. `outputs` sums them back into INFECTION.
- **ABM.** `TBExact` ([exact.py](../../tbode/exact.py)) draws each agent's class when it is first stepped after infection, stores it in `latent_class`, and samples the end-of-step state from that class's `expm(Q·dt)`. `ti_infected` is still set, but isn't used for progression.
- **Equilibrium start.** At the ODE equilibrium, each class has constant inflow and a constant per-capita exit rate, so time since infection within class m is exponential with rate inflow / stock (`TBODE.latent_age_rates`). `init_from_ode` gives each latent agent a class from the ODE's class proportions and a τ from its class's distribution. `tbsim.TB` reads τ; `TBExact` reads the class.
- **Checks.**
  - `tests/test_progression.py`: `TBODE`'s cohort progression matches `ode_cohort_progression`, which integrates the time-varying hazard directly, without classes.
  - `tests/test_exact.py`: `TBExact`'s one-step end states match the class-weighted `expm(Q·dt)`, and its cohort progression matches within sampling error.
  - Daily `tbsim.TB`, which uses the actual τ, matches the exact cohort curve (`run_cohort.py`).

## Reinfection of latent agents

In tbsim, reinfecting a latent agent resets τ to 0. In the mixture, the matching step is to redraw the agent's class with the same Poisson weights. That is also exact, because a fresh infection's class is independent of everything before it. `TBExact.set_prognoses` already does this by resetting `latent_class` to −1 on reinfection. It only takes effect once `rr_reinfection_inf` > 0, which is off in this experiment.

## If `k_asy` ≠ `k_non`

The single mixture needs `k_asy = k_non`, so that the two hazards combine into one; `TBODE` and `TBExact` raise an error otherwise. tb_LAI_TPT sets them equal.

If they differed, the same expansion would still give an exact model with two class indices, (m₁, m₂):

- assignment probability Poisson(m₁; `inf_non`/`k_non`) × Poisson(m₂; `inf_asy`/`k_asy`);
- progression rates m₁·`k_non` to NON_INFECTIOUS and m₂·`k_asy` to ASYMPTOMATIC.

That needs about as many compartments as the product of the two class counts. It isn't implemented.

## Why not a chain of time-since-infection stages

A chain of stages (0–1 month since infection, 1–2 months, …) was considered and rejected. Agents leave each stage after an exponentially distributed time, which blurs the timing of progression. With 1-month stages, cohort progression came out 4–7% too low; with 0.25-month stages, 1–2% too low. Near threshold, prevalence amplifies errors in progression, so that's too much. The mixture is exact with 7–8 compartments.
