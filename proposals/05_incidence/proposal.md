# Proposal 05: an incidence result in tbsim v4

**For:** tbsim v4 (`tbsim.TB`).
**Status:** draft, not sent.
**Track:** tbsim `origin/v4` @ `dcfb522`, starsim `origin/v4` @ `80c81f53`.
**Evidence:** [exp 19](../../experiments/19_v4_incidence/SUMMARY.md) (v4, both solvers); [exp 11](../../experiments/11_incidence_counting/SUMMARY.md) (v0.10 agents, snapshot against transition counts).

## Summary

This moves the incidence count we already use (CONVENTIONS.md 5; exps 11, 19–21) into tbsim. Today it lives in our own code: the ODE uses the same split of non-infectious TB in `tbode4.fastode`, and agents use a per-person flag in exp 19's `Incidence` analyzer.

- **The gap:** `tbsim.TB` has no incidence result. Users add up transitions themselves, and the obvious sum gives about twice the incidence that calibration targets mean, because it counts every entry to asymptomatic TB.
- **The proposal:** add a declared `incidence` result, counting new episodes of infectious TB. Split `non_infectious` into two states with identical rates, by whether the person has been infectious yet in the current episode.
- **Why this way:** it's exact in both solvers, counted from transitions, and needs no per-person flags. The split doesn't change the model's dynamics.

## The problem

People in tbsim's natural history move back and forth between asymptomatic and non-infectious TB within one episode of disease: `asy_non` is about 1.7 per year in tbsim's defaults and in calibrated sets. So an episode often enters asymptomatic TB more than once.

- **Every entry isn't incidence.** The sum `latent_asy + remote_asy + worsening` counts every entry to asymptomatic TB. About half of those entries are re-entries within the same episode: 46% with tbsim's defaults, 50% with a calibrated set, at 500 per 100,000 prevalence (exp 19).
  - Counting them doubles incidence, and halves prevalence ÷ incidence.
  - With tb_LAI_TPT's target of 1.57 years for prevalence ÷ incidence (1.20–2.05), the first-entry count gives 1.82, and the every-entry count about 0.9, far outside the target's range.
- **It's easy to get wrong in other ways too.** The output function in our rerun of the earlier comparisons on v4 (`outputs()` in `v4/common.py`) counts every entry, and with front-loaded latent progression (`k_latent` > 0) it also misses remote latent → asymptomatic TB (2% of entries in a calibrated set).
- **Agent-only workarounds don't carry over to the ODE.** A per-person "has been infectious this episode" flag works in agents (exp 19 uses one), but has no ODE form. Counting from end-of-step snapshots also misses brief episodes: 3.3% at monthly steps in v0.10 (exp 11).

## The proposal

**Definition:** an episode of TB disease starts at the first entry to asymptomatic TB, from latent infection or from non-infectious TB in someone not yet infectious in this episode. It ends at recovery, cure or death. Moving between asymptomatic and non-infectious TB within an episode doesn't start a new one. Incidence is new episodes per person-year.

**States:** split `non_infectious` into two states with identical rates:

| State | Who | Entered from |
|---|---|---|
| `non_infectious` | Not yet infectious in this episode | Latent infection (`latent_non`, `remote_non`) |
| `non_infectious_post` | Infectious earlier in this episode | Asymptomatic TB (`improvement`) |

- Both have the same exits at the same rates: `recovery` to `recovered`, `worsening` to `asymptomatic`, and background death.
- Only `worsening` from `non_infectious` starts an episode.
- A `non_infectious_total` result (or keep the name `non_infectious` for the result and rename the states) keeps existing outputs unchanged.

**The result:**

```python
self.define_results(incidence = self.latent_asy + self.remote_asy + self.worsening + self.relapse)  # relapse: with p_relapse > 0
```

A result that sums transitions runs in both solvers and is counted from transitions (checked: `latent_asy + worsening` declared this way equals the sum of the two transitions' own results, as agents and as an ODE). `worsening` here is only the transition out of `non_infectious`; the one out of `non_infectious_post` has its own name, e.g. `worsening_post`.

**Relapse** (with `p_relapse` > 0): `relapsing` → `symptomatic` counts as a new episode, so `self.relapse` is added to the sum. tb_LAI_TPT counts relapses in its incidence (`IncidenceTracker`), and WHO's incidence estimates include relapses. A failed course of treatment that returns the person to the state they started from continues the same episode, and isn't counted.

**Why it's exact:** the two non-infectious states have identical rates, so their total behaves exactly as the single state does now. Only the bookkeeping changes. Our standalone copy of v4's ODE with this split reproduces v4's ODE to 1e-14 at equilibrium. As agents, a flag-based count of the same definition matches the split ODE to within 0.4% per prevalent case, at monthly and weekly steps (exp 19).

## Consequences, and choices for the maintainers

- **Starting states:** `init=dict(non_infectious=...)` would need a rule. Options:
  - accept both keys, with `non_infectious` meaning "not yet infectious";
  - or take a total and split it at the equilibrium ratio when starting at equilibrium.

  Starting everyone in `non_infectious` slightly overcounts incidence early in a run, until those who were already infectious recover or progress.
- **Treatment from non-infectious TB:** if treatment can start from both non-infectious states (e.g. active case finding), each needs its own treatment state, so failure returns people to the right one. This matches [proposal 04](../04_treatment/proposal.md)'s "failure back to the starting state".
- **Drug resistance:** the split carries over to strains like any other state. Incidence by strain then comes from the same transitions.
- **An alternative implementation:** a two-level dimension ("infectious this episode: no / yes"), changed by `move=` on `improvement` and reset on recovery. It would avoid naming a second state, but we haven't tested it.
- **A symptomatic-based definition** (first entry to symptomatic TB per episode, tb_LAI_TPT's `sym_only` option) would need asymptomatic and non-infectious TB each split by "symptomatic yet in this episode". We don't propose it as the default, but the same approach extends to it.

## Not proposed

- **Leaving out brief subclinical episodes** from incidence. If wanted, define the rule by duration in both solvers, not by step boundaries.
- **Notifications:** with treatment declared as transitions, they are the transitions into treatment, and need no new states.
