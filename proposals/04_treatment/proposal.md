# Proposal 04: treatment outcomes and relapse in tbsim v4

**For:** tbsim v4 (`tbsim.TB`, `tbsim.regimen`).
**Status:** draft, not sent. We are implementing the same design in our own code first ([CONVENTIONS.md, section 9](../../CONVENTIONS.md#9-treatment)), and will add evidence from [exp 23](../../experiments/23_treatment/README.md), which tests it.
**Track:** tbsim `origin/v4` @ `dcfb522`, starsim `origin/v4` @ `80c81f53`.

## Summary

Two changes to how `tbsim.TB` ends a course of treatment, and a documentation fix:

1. **A failed course returns the person to the state they started treatment from,** not always to symptomatic TB. Recommended.
2. **An option to draw the outcome when treatment starts,** rather than from competing completion and failure rates. For discussion: it matters only in some uses.
3. **Documentation:** `dur_relapse` is described as "the time from the end of treatment to relapse", but as a plain duration Starsim takes it as the mean of an exponential, so relapse is at a rate (mean 1.5 years) in both solvers. That's the behaviour we want; the docstring should say so, and that `ss.constant()` gives a fixed time.

## How treatment works now

- **Start:** `diagnosis` moves non-infectious, asymptomatic or symptomatic TB into `treatment`, at rate `theta` (0 by default) or as `tbsim.care()` delivers it (symptomatic TB by default, via `where=`).
- **End:** `completion` (to `treated`, or a share `p_relapse` to `relapsing`) and `failure` (to `symptomatic`) compete, at rates `delta` = s/D and `phi` = (1 − s)/D, with success probability s = cure × adherence and mean course length D (`tbsim.regimen`).
- **Relapse:** `relapsing` → `symptomatic` at rate 1/`dur_relapse` (mean 1.5 years).

## 1. Failure back to the starting state

**The problem:** whenever treatment starts from non-infectious or asymptomatic TB (with `theta` > 0, `care(where=...)` wider than symptomatic TB, or active case finding), a failure moves the person to symptomatic TB. Failing treatment then makes people worse than if they hadn't been treated: infectious at the full rate, and at risk of TB death. In the ODE this is part of the model's structure, not a detail of delivery. tbsim v0.10's `Tx` returned people to their prior state, and so does v4's own reference resistance model, `LSHTMResistant` ("a failure returns people to where they were").

**The proposal:** one treatment state for each state treatment can start from, as `LSHTMResistant` does with `treat_asy` and `treat_sym`:

- `diagnosis` becomes one transition per starting state (e.g. `diagnosis_non`, `diagnosis_asy`, `diagnosis_sym`), and `failure` from each treatment state goes back to its starting state.
- A `treatment` result can still report the total.
- This is exact in both solvers and needs no per-person memory, so agents and the ODE stay the same model.
- Only states that treatment can start from need a treatment state. A model that treats symptomatic TB only has one, as now.

## 2. The outcome drawn at the start

**What it would be:** a person starting treatment enters `treatment_success` with probability s, or `treatment_fail` otherwise. Each ends when the course ends: success to `treated` (or `relapsing`), failure back to the starting state.

**When it matters:**

- **Time-varying cure or adherence:** drawn at the start, the outcome uses the efficacy when treatment started (as v0.10's `Tx` does). With competing rates, it uses the efficacy at each moment of the course.
- **Different course lengths by outcome,** e.g. people lost to follow-up leaving early. Competing rates force both outcomes to share one course length.
**When it doesn't:** with constant parameters, it gives the same outcomes and times as competing rates (for a course ending at a rate) or as a branch at completion, one transition to `treated` or the starting state with probability s (for a fixed course).

**The cost:** two treatment states per starting state instead of one.

## Not proposed yet

- **Death and loss to follow-up during treatment.** We plan to add both in our own code, and may propose them later.
- **Whether the course should be fixed or a rate.** Our experiment compares a fixed 6 months (`ss.constant`, which v4's ODE runs as a conveyor of whole steps) with a rate. At equilibrium, with no death on treatment, we expect a difference of order 10⁻⁵ of prevalence. During scale-up, we expect a difference only in the timing of returns from treatment.
