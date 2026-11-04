# Latent progression and reinfection in v4: options

*A comparison to support a decision, not a recommendation. Numbers are for set B (tb_LAI_TPT's best-fit natural history, `k` = 5.846), from [`compare_structures.py`](compare_structures.py) (results in `compare_structures.json`) and Cliff's v4 redo ([`v4/v4_findings.qmd`](../../v4/v4_findings.qmd)), confirmed at pinned versions in exp 13.*

## What has to be decided

tb_LAI_TPT (on tbsim v0.10) makes progression from latent infection front-loaded. The rates to non-infectious and asymptomatic TB decline as r·e^(−kτ) with time since infection τ, and reinfection resets τ to 0.

tbsim v4 replaced this with `k_latent`: two latent states, recent and remote, with constant rates. Reinfection moves a remote infection back to recent, and does nothing to a recent one. We want to propose a structure for v4 that the ODE and agents both run, and that is consistent with the earlier work. This note compares the candidates against four criteria.

## Criteria

1. **Shape of progression.** Must it reproduce r·e^(−kτ) exactly, or only the data it was fitted to? The data are the progression targets (ProgressionY1 = 0.0523, ProgressionRatio = 0.663) and the household-contact studies behind them.
2. **Reinfection.** Should a reinfection reset the clock (as in v0.10), or something else? How much does the choice change outcomes?
3. **Exactness under both solvers.** Can the ODE and agents run the same model with no approximation?
4. **Size and identifiability.** How many compartments it adds (multiplied by any stratification, e.g. age), and how many parameters a calibration has to identify.

## Candidates

### A. Poisson-mixture latent classes (exact r·e^(−kτ))

A new infection joins class m with probability Poisson(m; r₀/k), r₀ = `inf_non` + `inf_asy`, and progresses at the constant rate m·k. Class 0 never progresses.

- **Shape:** exactly the r·e^(−kτ) model, for any `k`.
- **Reinfection:** a reinfected latent person is reassigned across the classes with the same weights. That is exactly a clock reset.
- **Solvers:** exact under both. Cliff declared it in v4 as an ordinary model, `TBClasses` (about 60 lines). Agents match the ODE: set B 1.003 ± 0.008 / 0.993 ± 0.004 at monthly / weekly `dt` (exp 13, 40 / 200 reps).
- **Size:** 8 latent states for set B at a truncation of 1e-10 (4–5 are enough to 1e-6), against 1–2 for the others. One parameter, `k`, as in v0.10. So tb_LAI_TPT's `k` carries over unchanged.
- **Interpretation.** It is a mixture of people who will and won't progress, rather than one person whose rate declines. The two are indistinguishable in every distribution of progression times, with or without reinfection.

### B. v4's two states (`k_latent`), as they are

Recent infections progress at `inf_non` / `inf_asy` and become remote at rate `k_latent`. Remote infections progress at `inf_non_remote` / `inf_asy_remote`; set to 0 here, to mimic r·e^(−kτ).

- **Shape:** close but not equal. A recent infection's progression density is r₀·e^(−(r₀ + `k_latent`)·τ), against about r₀·e^(−kτ) for the exact model. With `k_latent` = k, about 4–5% fewer people progress at every time (table below).
- **Reinfection:** remote → recent. Reinfecting a recent infection does nothing, where the exact model resets the clock. Recent infections are a small share of latent ones, since they last about 2 months, so this matters little (table below).
- **Solvers:** exact under both (set B 1.01 ± 0.01 / 1.01 ± 0.01 at monthly / weekly `dt`; exp 13, 40 reps).
- **Size:** 2 latent states. One parameter (`k_latent`), or three if progression from the remote state is used.
- **Pitfall:** `inf_non_remote=None` means "the same as `inf_non`", so by default only progression to asymptomatic TB is front-loaded.
- **An option the exact model doesn't have.** Progression from the remote state (`inf_*_remote` > 0) gives lifelong reactivation at a low rate. In r·e^(−kτ), the rate goes to 0, so after a few years nobody progresses.

### C. v4's two states with `k_latent` refitted

The same as B, with `k_latent` chosen so that ProgressionY1 matches the exact model: **`k_latent` = 5.533** rather than 5.846.

- **Shape:** within about 1% of the exact model at every time from 3 months to 50 years (table). ProgressionRatio is 0.675 vs 0.674.
- **Reinfection:** as B. The progression it adds matches the exact model's clock reset to 0.0001 at set B's force of infection, and to 0.0001 even at 5× (table).
- **Solvers, size:** as B. `k_latent` is a refit of `k`, so tb_LAI_TPT's `k` doesn't carry over as is.

### D. A chain of recent stages (Erlang or Coxian)

Several recent stages with declining progression rates, before remote. This generalises B and C: with more stages it can approach any shape, including r·e^(−kτ). It is exact under both solvers, since every rate is constant. It adds stages and parameters: either fitted rates, or rates tied to `k` by a formula. Not computed here; given how close C already is, more stages may not be needed for set B.

### E. A non-exponential dwell time for progression

v4 lets a transition take a dwell-time distribution (`ss.Transition(..., dur=...)`). Agents sample it exactly, but the ODE approximates it with Erlang stages (or `solver='ide'`), so it fails criterion 3. It also fits poorly:
- **Defective.** Progression is defective: most people never progress, and clearance competes with it.
- **Close to exponential.** The time to progression among those who progress is close to exponential (coefficient of variation ≈ 1), which Erlang stages represent badly.

## Numbers

**Cohort: the share that has ever reached asymptomatic TB.** Set B, everyone latent at t = 0, no reinfection or deaths, as tb_LAI_TPT's targets are defined:

| Structure | 3 mo | 6 mo | 1 y | 2 y | 5 y | 10 y | 50 y | ProgressionY1 | ProgressionRatio |
|---|---|---|---|---|---|---|---|---|---|
| A. Exact (classes) | 0.0346 | 0.0439 | 0.0511 | 0.0597 | 0.0715 | 0.0757 | 0.0764 | 0.0511 | 0.674 |
| B. Two states, `k_latent` = 5.846 | 0.0340 | 0.0424 | 0.0489 | 0.0570 | 0.0683 | 0.0723 | 0.0729 | 0.0489 | 0.676 |
| C. Two states, `k_latent` = 5.533 | 0.0349 | 0.0441 | 0.0511 | 0.0596 | 0.0714 | 0.0756 | 0.0763 | 0.0511 | 0.675 |

Targets: 0.0523 / 0.663. None of the three hits them exactly; set B's `k` was fitted with v0.10's biased monthly `tbsim.TB` (exp 08).

**Reinfection: what each structure's reinfection rule adds.** The same cohort, with background deaths (1/70 per year) and a constant force of infection λ. People no longer infected are reinfected at π·λ (π = 0.24), as new infections. The table shows the extra share ever reaching asymptomatic TB when latent people are also reinfected at π·λ by the structure's own rule, compared with not at all:

| Structure | λ = 0.12/y (set B at 500/100k), by 10 y | by 30 y | λ = 0.60/y (5×), by 10 y | by 30 y |
|---|---|---|---|---|
| A. Exact: clock reset | +0.0035 | +0.0046 | +0.0222 | +0.0398 |
| B. Two states: remote → recent | +0.0033 | +0.0045 | +0.0216 | +0.0390 |
| C. Two states, refitted | +0.0035 | +0.0046 | +0.0223 | +0.0399 |

So the two reinfection rules give the same answer, to within the difference in progression itself. That isn't surprising. Under the exact model, reinfecting someone infected 2 months ago changes little, because their rate has barely fallen. Under the two-state model, those are the recent infections whose reinfection does nothing. The bulk of latent reinfection acts on long-standing infections, and both rules restart their progression risk in the same way.

## Summary against the criteria

| | A. Classes | B. Two states, as is | C. Two states, refitted | D. Stage chain | E. Dwell time |
|---|---|---|---|---|---|
| 1. Shape | Exact | 4–5% low | Within ~1% (set B) | Any, with enough stages | Poor fit |
| 2. Reinfection | Exact clock reset | Remote → recent; matches clock reset to the progression difference | Matches clock reset (to 0.0001) | Choice of rule | — |
| 3. Both solvers exact | Yes | Yes | Yes | Yes | No (ODE approximates) |
| 4. Size, parameters | 8 states; `k` | 2 states; `k_latent` | 2 states; `k_latent` (refit of `k`) | 3+ states; more | 1 transition |
| Late reactivation possible | No (rate → 0) | Yes, with `inf_*_remote` | Yes, with `inf_*_remote` | Yes | — |
| Already in tbsim v4 | No (Cliff's `TBClasses`) | Yes | Yes | No | Partly |

## Questions this leaves for the decision

- **Is "within ~1% of r·e^(−kτ)" consistent enough with the earlier work?** That is option C. Or does the proposal need the exact form (A), e.g. because tb_LAI_TPT's `k` should carry over as is, or for other parameter sets?
- **Is C's closeness general?** It depends on r₀/k being small (0.12 for set B, 0.04 for set A). A check over tb_LAI_TPT's posterior samples, not just the best fit, would show whether it holds where calibrations put the parameters.
- **Should the model allow late reactivation?** Only the remote-state structures can.
- **Reinfection of recent infections.** Is v4's "does nothing" acceptable as a modelling statement, given that it makes no measurable difference here? Or should the proposal ask for a clock reset in principle?
- **What does the size cost once the model is stratified,** e.g. by age with 16 groups: 8 × 16 latent compartments for A, against 2 × 16 for B or C?
