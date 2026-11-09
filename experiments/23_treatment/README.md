# 23: Treatment in the v4 ODE and agents: the design, and a fixed course against a rate

**Status:** done; see [SUMMARY.md](SUMMARY.md). Preliminary checks of how v4 and tb_LAI_TPT behave are in [checks/](checks/README.md). The agents' runs were restarted once, with Starsim's kernels limited to one thread per process so the run stays within 4 cores.
**Track:** v4 (`v4/.venv`; tbsim `origin/v4` @ `dcfb522`, starsim `origin/v4` @ `80c81f53`). One check also runs v0.10 `Tx` as tb_LAI_TPT ran it.

## Background

**The treatment design** (agreed with the user, 2026-11-09; the rules are in [CONVENTIONS.md, section 9](../../CONVENTIONS.md#9-treatment)):

- **Starting treatment:** a declared transition from symptomatic TB at θ(t) = −ln(1 − c(t)) per year, with c the annual coverage.
- **The outcome:** drawn at the start, by entering `treatment_success` with probability s = adherence × efficacy, or `treatment_fail` otherwise.
- **Failure:** back to the state treatment started from (symptomatic TB).
- **The course:** either a fixed 6 months or a rate with mean 6 months. This experiment compares the two.
- **Relapse:** a share `p_relapse` of successes goes to `relapsing`, then to symptomatic TB at rate 1/1.5 per year.
- **On treatment:** nobody transmits or is reinfected. People whose treatment will fail keep their TB death rate (`sym_dead`, since they started from symptomatic TB) for the course (user, 2026-11-09). People who will be cured don't die of TB.

**How it compares with tb_LAI_TPT and with tbsim v4 as it stands.** tb_LAI_TPT is the agent-based model this repo aims to approximate. It runs tbsim v0.10's `Tx` through its own `SymptomaticTxDelivery` and `TxVaryEfficacy`, and its calibration `hm_real_20260728` ran tbsim `85a977a` with starsim `8e7868a`. Rows marked † are checked by running the model ([checks/](checks/README.md)); the rest are from reading the code.

| | tb_LAI_TPT (v0.10) | tbsim v4 as it stands | This design |
|---|---|---|---|
| Who starts | Symptomatic, not on treatment; chance per step 1 − (1 − c)^dt † | Whoever `tbsim.care()` reaches (symptomatic by default); the bare `theta` rate applies to non-infectious, asymptomatic and symptomatic TB | Symptomatic only, at θ = −ln(1 − c), which reproduces tb_LAI_TPT's per-step chance at any `dt` |
| Outcome | Drawn at the start, with efficacy read at the start † | Completion and failure compete, at rates s/D and (1 − s)/D (`tbsim.regimen`) | Drawn at the start † |
| Course | Fixed 180 days, completing at the 6th monthly step † | A rate, mean D | Fixed 6 months † or a rate with mean 6 months |
| After success | `CLEARED`, reinfection risk `rr_reinfection_treat` | `treated`, `rr_reinfection_treat` | `treated` (risk π) |
| After failure | The state treatment started from | Symptomatic, whatever the starting state | The state treatment started from |
| Relapse | A share `p_relapse` of successes. Meant to be 1.5 years after treatment ends, but it fires at the next step, about 1 month after cure † | A share `p_relapse` → `relapsing` → symptomatic at rate 1/`dur_relapse`, mean 1.5 years † | As v4 |
| Restart in the step of failure | Not possible: eligibility is read before failures return | Possible, at the continuous-time rate † | As v4 |
| On treatment | Not infectious, no TB death, not reinfected | The same | Not infectious, not reinfected; TB death at `sym_dead` for those whose treatment will fail † |

**tb_LAI_TPT's relapse delay** comes from a unit error in v0.10's `Tx`. It divides the delay by the days per step as if it were in days, but Starsim has already converted the default `ss.constant(v=ss.years(1.5))` to the sim's unit (18 months). So its calibrated `p_relapse` (0.034) means 3.4% of cures return to symptomatic TB a month later, which acts much like a slightly lower efficacy. v4 doesn't have the error, and it isn't pursued further here.

**tb_LAI_TPT's best fit** (row 0), used as the treatment values here:

- **Coverage:** 0 to 1994, 0.548 in 2010 and 0.677 in 2024, linear between.
- **Efficacy:** 0.45 in 1995 and 0.53 in 2024.
- **Adherence:** 1.
- **`p_relapse`:** 0.034.

**Fixed course or rate: what to expect.** Nobody on treatment is infectious. For people who will be cured, only background death (μ = 1/70 per year) competes with the end of the course. For people whose treatment will fail, TB death (d = `sym_dead`) competes too.

- **Cures, almost no difference:** the share finishing is e^(−μD) = 0.99288 for a fixed course and 1/(1 + μD) = 0.99291 for a rate.
- **Failures, a real difference:** with d = 0.396 per year (exp 21's set 1), the share dying of TB on treatment is 1 − e^(−dD) = 18.0% for a fixed course and d/(d + 1/D) = 16.5% for a rate, ignoring background death. With efficacy about 0.5, about half of all courses fail, so about 9.0% against 8.3% of all courses end in TB death. So TB mortality, and mortality ÷ incidence, should differ at equilibrium by a little under 1% of courses. This experiment measures the effect on mortality and prevalence.
- **During scale-up, the timing of returns differs:** with a rate, 28% of failures return within 2 months; with a fixed course, all return at 6 months. Coverage rises over 15–30 years, so the lag should be a small share of the change.
- **Loss to follow-up on treatment** would make the course length matter more, because it competes with the course too. It's to be added later, not here.

## Questions

1. **Does the design behave as intended in agents?** Check every pitfall listed below.
2. **Do v4's agents match its ODE with treatment,** with the course either fixed or a rate, at equilibrium and through tb_LAI_TPT's coverage scale-up?
3. **How much does a fixed 6-month course differ from a rate with the same mean?** Compare at equilibrium and during scale-up, in prevalence, incidence, notifications and mortality. A small difference would let either serve, and the choice could follow convenience.

## Implementation (before running)

- **Agents and v4's ODE:** a `tbsim.TB` subclass in `tbode4`, `TBTreat`. It replaces `diagnosis`, `failure` and `completion` with the design's transitions, with options `course='fixed'|'rate'`. Coverage and efficacy vary over time as `ss.timeseries` with `ss.modify`.
- **TB death on treatment:** `ss.death(treatment_fail, rate=p.sym_dead)`, a transition of its own (`tb_death_tx`), so a time-varying `sym_dead` applies to it too. Two things go with it:
  - **Replacement births:** `tbsim.TB` replaces TB deaths with births at `sym_dead × symptomatic`. That must also include people on failing treatment, or the population shrinks.
  - **Mortality:** TB mortality is `tb_death + tb_death_tx`.
- **`fastode`:** add the treatment states and flows, for calibration later. Test it against v4's ODE with treatment, at both course options, before using it (CONVENTIONS.md 8).
- **Incidence:** the first entry to asymptomatic TB per episode (CONVENTIONS.md 5), plus relapses (user, 2026-11-09).
  - **Why relapses need adding:** a relapse goes straight from `relapsing` to symptomatic TB, so the first-entry rule alone would never count it. tb_LAI_TPT counts relapses as new episodes, and WHO's incidence estimates include them.
  - **Agents:** extend exp 19's `Incidence` analyzer to count `relapse` transitions and set the person's "has been asymptomatic this episode" flag.
  - **ODE:** add the relapse flow to `fastode`'s incidence. People relapsing go to symptomatic TB, so on remission and improvement they join NON1, as now.
  - **Not new episodes:** a failed course returns the person to symptomatic TB within the same episode. Cure ends the episode, so a reinfection that progresses later is counted from latent infection, as now.
  - **Check:** agents ÷ ODE incidence per prevalent case. Report the share of incidence from relapse.
- **Notifications:** counted from transitions into treatment.

## Design

**Parameters:** exp 21's set 1 natural history (the starter script's set), with tb_LAI_TPT's treatment values (above).

**Scenarios:**

- **Constant coverage, at equilibrium (`eq`):** tb_LAI_TPT's 2024 values (c = 0.677, s = 0.53, `p_relapse` = 0.034), 2000–2049, compared over 2005–2049.
  - **The starting state:** `fastode`'s treated equilibrium, with the course as a rate for the `rate` arm, and as 50 Erlang stages, nearly fixed, for the `fixed` arm.
  - **The transient:** v4's `init=` can't set a time already on treatment, so everyone on a fixed course starts at its beginning. The ODE and agents start from the same state, so they share that transient, and the first 5 years are left out.
- **Scale-up (`ramp`):** start at the untreated equilibrium in 1980, then follow tb_LAI_TPT's coverage and efficacy trajectories from 1995 to 2024, and hold them to 2035. Compared over 1995–2035.
  - **Coverage:** linear from 0 at the end of 1994 to 0.548 in 2010 and 0.677 in 2024, given to the model as θ = −ln(1 − c) at monthly knots.
  - **Efficacy:** 0.45 from 1995, linear to 0.53 in 2024.

**Arms**, each in both scenarios, at monthly and weekly `dt`:

| Arm | Course | Solver |
|---|---|---|
| `rate` | Rate, mean 6 months | ODE, and agents |
| `fixed` | Fixed 6 months (`ss.constant`; the ODE runs it as whole-step slots, exact to the step) | ODE, and agents |
| `stages_k` (ODE only) | Erlang, k = 1, 2, 4 and 12 stages, mean 6 months | ODE |

The `stages_k` arms show how the difference between a fixed course and a rate grows with the spread of course lengths. k = 1 is the rate.

**ODE equilibria:** every ODE arm is also run for 400 years at the constant 2024 values, and the last year is taken as its equilibrium. That gives fixed ÷ rate at equilibrium without the starting transient.

**Agents:** 200,000, 40 reps per course, scenario and `dt`, with new seeds (23000 + 1000·scenario + 100·arm + rep). The treated equilibrium's prevalence is about 80 per 100,000, so 200,000 agents hold about 160 people with active TB. At 50,000, a 2-rep trial run varied 39–118 per 100,000 between reps.

**Outputs**, per 100,000 per year:

- prevalence (asymptomatic + symptomatic TB);
- incidence (first entry to asymptomatic TB per episode);
- notifications, and notifications ÷ incidence;
- TB mortality, and mortality ÷ incidence;
- the number on treatment, and the time on treatment.

Report agents ÷ ODE (mean ± SE, 10–90% interval), and `fixed` ÷ `rate` in the ODE, over time for the scale-up.

**Figures** (yearly values per 100,000; the data behind each saved with `sc.save`):

- **Over time, one figure per scenario and `dt`:** five panels:
  - prevalence;
  - incidence, with the part from relapse shaded;
  - notifications: courses started per year, including restarts after failure or relapse;
  - the number on treatment: a stock, at any time about notifications × 0.5 years at equilibrium, and lagging notifications during scale-up;
  - TB mortality, with the part on treatment shaded.

  Each panel shows the ODE for `rate` and `fixed` as lines, and agents as the mean of the reps with a 10–90% band. The scale-up figure marks 1995 and 2024.
- **Fixed ÷ rate, ODE, over time:** each of the five outputs, including the `stages_k` arms, to show where and when a fixed course differs from a rate.
- **Agents ÷ ODE:** each output by arm, mean ± SE, at monthly and weekly `dt`.

**Pitfall checks** (agents, `rate` and `fixed`, from per-step transition records). v0.10's delivery machinery holds per-person flags such as `diagnosed`, `tested`, `pending_success`, `pending_relapse` and `prior_state`, and most of its pitfalls came from them going stale. This design has no such flags; each check confirms that.

| Pitfall | Check |
|---|---|
| Starting treatment again while on it, or every step | No transitions into treatment from a treatment state. Starts ÷ symptomatic person-time = θ |
| Restarting straight after failing | The share of failures that restart in the same step is 1 − e^(−θ·dt/2), as in continuous time (about 5% at θ = 1.1 per year, monthly), and no more |
| Outcomes or relapse carried into a later course | `relapsing` is entered only from `treatment_success`. The success share of each course doesn't depend on earlier courses |
| Time-varying efficacy applied at the wrong time | The success share by year of starting matches s(t) at the start |
| Treating people without symptomatic TB | No starts from other states |
| TB death or transmission on treatment | TB deaths only from `treatment_fail`, at `sym_dead`, matching the ODE's share per course. None from `treatment_success`. The force of infection excludes both. The population stays constant |
| Annual coverage used as a rate | Starts per symptomatic person-year match θ = −ln(1 − c), not c |
| Course length | The time on treatment is exactly 6 months (`fixed`), or exponential with mean 6 months (`rate`) |
| Results read a step out | Results at index i are dated at the end of step i. Coverage changes apply from the first step starting on or after their date |
| A restart counted as a new episode | Notifications ÷ starts = 1. Incidence is the same with and without restarts |

**Then:** recalibrate with treatment, bringing back mortality ÷ incidence (a later experiment).
