# 33: Summary

**Track:** v4, trial environment in this folder: starsim `origin/v4` @ `0089ac0f` (4.0.0rc5, the merge of Starsim PR #1625) and tbsim `origin/v4` @ `a2ae583`, against `v4/`'s pins (starsim `d83f4f87`, tbsim `619c296`). Run on the laptop (4 cores), 2026-11-20 to 2026-11-22. `versions.json` holds the trial environment's versions; `pins_ode_*.obj` hold each pin's.

**Question:** does tbsim `a2ae583` work with starsim `0089ac0f`, and with Starsim's fixes to its demography (proposal 06), do agents match the ODE without our workarounds?

**Result:** yes on every count.

- **The bump:** the pins work together. v4's ODE gives the same numbers as at the current pins to 8e-16, with flows dated as before. One `tbode4` line needed changing, for a stricter check in `ss.Pregnancy`.
- **Agents without workarounds:** with Starsim's defaults, agents' demography matches the ODE's as closely as exp 29's corrected setup did, and TB with demography matches as closely as in exp 30. Every proposal 06 item is fixed.
- **What changes:** births are now female at the 1950 fraction at ages 0–4 (49.06%), not the overall 48.63%. That raises the population by 1.0% by 2020, but moves exp 32's calibration targets by at most 0.4%.

## Stage 1: the bump

- **Building:** tbsim `a2ae583` with starsim `0089ac0f` builds and runs.
  - **One change to `tbode4`:** `ss.Pregnancy` now refuses a fertility rate above 1 per person per year without `rate_per`, so `tbode4.gestation()` gives `rate_per=1000`. This works at both pins.
  - **Renamed results:** `died` is now `death`, there is a top-level `birth`, and `emigrated` is gone. Exp 28's demography outputs read `died`, so this experiment reads either name.
- **`v4/tests/`:** 20 of 21 pass (`tests_starsim_4.0.0rc5.log`).
  - **The failure is the new default for births' sex:** `test_copy_matches_v4_ode` (`fastode_age` against v4's ODE) had population 0.4% apart by 1994. It passes with v4's ODE given `ss.birth(p_female=0.4863)`, the share `fastode_age` used, and with `fastode_age` given Starsim's new default (stage 4). The ODE now honours `ss.birth`'s `p_female`, which it ignored before.
- **The same ODE at both pins:** exp 30's v4 ODE (TB, treatment from 1995, WPP demography, quarter-year groups, monthly steps, 1950–2041) was run in each environment with births' sex fixed to the same share.
  - **Result:** every state and flow at every step agrees to within 8e-16 (`pins_compare.csv`).
  - **The comparison would catch a re-dating of flows:** comparing with the run shifted by one step gives differences of up to 380%. So flows are still indexed by the step's start (CONVENTIONS.md 7).

## Stage 2: demography alone, Starsim's defaults against exp 29's corrected setup

Exp 29's setting: India's WPP 2024 rates, 1950–2050, monthly steps, deaths by age and sex, no TB. Agents are 10 reps of 100,000 at the start (`demog_summary.log`).

| Arm | Setup |
|---|---|
| Q, E | ODE with quarter- and eighth-year groups, split by sex as the 1950 population is at each age (Starsim's default), births at the default share |
| Qj | ODE, quarter-year: exp 29's R4x (`joint_sex()`, births at 48.635%) |
| D | Agents with Starsim's defaults: `ss.population`, `ss.Pregnancy` on WPP's table as given (`fertility_age='delivery'`, default start of pregnancies and births' sex) |
| W | Agents as exp 29's P2b (shifted table, births at 48.635%), read at the age at conception so it isn't shifted twice, with Starsim's new start of pregnancies |

**Population ÷ 1950, ratios between arms** (mean ± SE over 10 reps; "exact" is the ODE extrapolated to exact aging, 2 × E − Q):

| Year | D ÷ exact | W ÷ exp 29's P2b | Qj ÷ exp 29's R4x | Q ÷ Qj (births' sex) |
|---|---|---|---|---|
| 1960 | 0.9994 ± 0.0009 | 1.0008 ± 0.0008 | 1.000000 | 1.000 |
| 1980 | 1.0005 ± 0.0015 | 1.0009 ± 0.0015 | 1.000000 | 1.001 |
| 2000 | 1.0003 ± 0.0025 | 1.0003 ± 0.0017 | 1.000000 | 1.006 |
| 2020 | 1.0003 ± 0.0028 | 0.9998 ± 0.0021 | 1.000000 | 1.010 |
| 2050 | 0.9993 ± 0.0032 | 0.9986 ± 0.0022 | 1.000000 | 1.017 |

- **Starsim's defaults match without workarounds:**
  - **Population:** agents with every default (D) are within 0.1% of the ODE with exact aging through 2050, as exp 29's corrected agents (P2b) were.
  - **Age shares:** within 1% in each broad age group: 0.993–1.012 from 1970 (`demog_shares.csv`).
  - **Fertility:** total fertility is 1.000 of WPP's, realized fertility at ages 15–44 within about 2% of WPP's rates, and the mean age at childbearing within 0.03 years (`demog_fertility.csv`).
- **Nothing else changed:**
  - **The ODE:** Qj equals exp 29's R4x exactly.
  - **Agents:** W equals exp 29's P2b within noise. So the new start of pregnancies (h·d with the unrounded gestation) is exp 29's `PregnancyStart`, and fertility read at the age at conception is the old behaviour.
- **One difference in age shares has a likely explanation:**
  - **The size:** W's share at ages 0–14 is 0.993–0.995 of P2b's from 1970 on, with 15–29 correspondingly higher.
  - **The explanation:** that is about one monthly birth cohort (1/180 of ages 0–14). It fits exp 29's 32-bit ages, which put each cohort that had just reached 15 just below it, so it was counted in 0–14.
  - **Not checked:** D's 0–14 share is similarly 0.993–0.996 of the exact ODE's from 1970, which is within the 1% criterion but not explained.
- **The other fixes all hold:**
  - **Losses:** `ss.Pregnancy`'s loss results equal the losses counted by balance exactly: 15,579 in D and 15,629 in W over all reps.
  - **The warning:** it no longer appears.
  - **Ages:** agents' ages are 64-bit. Of about 500,000 born in each run, none sits just below a whole year, and adding 1e-3 years moves none between 5-year groups.
- **Births' sex:** the ODE's births are 49.06% girls by default (Q), the 1950 fraction at ages 0–4, against 48.64% in Qj. This alone raises the population by 1.0% by 2020 and 1.7% by 2050 (Q ÷ Qj).

## Stage 3: TB with demography, against exp 30

Exp 30's model and starting tables (`TBTreat`, sets 0 and 6, agents' drawn start), with Starsim's default demography in both solvers. 20 agent reps per set and step, 100,000 agents (`tb_summary.log`, `tb_by_window.csv`).

**Agents ÷ the ODE with exact aging** (mean ± SE of window means; weekly steps, on which exp 30's criteria were set):

| Output | Set 0, 1995–2024 | Set 0, 2025–40 | Set 6, 1995–2024 | Set 6, 2025–40 |
|---|---|---|---|---|
| Prevalence | 1.002 ± 0.004 | 0.999 ± 0.008 | 0.999 ± 0.006 | 0.995 ± 0.007 |
| Incidence per case | 0.999 ± 0.002 | 0.999 ± 0.002 | 0.999 ± 0.002 | 1.000 ± 0.002 |
| Notifications per case | 0.999 ± 0.004 | 0.999 ± 0.004 | 0.998 ± 0.003 | 1.002 ± 0.003 |
| Mortality per case | 1.003 ± 0.003 | 1.002 ± 0.005 | 1.004 ± 0.004 | 1.004 ± 0.007 |
| Population | 0.997 ± 0.001 | 0.997 ± 0.001 | 0.997 ± 0.001 | 0.996 ± 0.002 |
| Prevalence, five age groups | 0.999–1.010 | 0.996–1.002 | 0.996–1.004 | 0.990–0.999 |

- **Exp 30's criteria are met:**
  - **Prevalence and ratios per case:** all-ages prevalence and every ratio per prevalent case are within 1%, and within about 2 SE in the later windows at weekly steps.
  - **Population:** within 0.4%, though 2–3 SE low (exp 30: 0.997–0.999).
  - **Prevalence by age group:** within 1% at weekly steps.
- **Monthly steps:**
  - **As in exp 30:** age groups are within 2.3%, recent latent infection reads high in end-of-step snapshots (1.6 for set 0), and agents drift slowly from the ODE's start (exp 27).
  - **Before 1995:** set 0's prevalence is 0.983 ± 0.006, against exp 30's 1.009 ± 0.005, a difference of about 3 SE between two runs with different seeds. Weekly agents show no such gap (1.003 ± 0.005).
- **The pitfall checks are as in exp 30, per rep, with one change in what they see.**
  - **The change:** background deaths of people on a curing course (`bg_death_from_success`) were 0 in exp 30 and are 320–475 here.
  - **Why:** the check reads deaths from the TB module's `deaths` transition. With `mu` = 0 that read nothing at `d83f4f87`; at rc5, deaths from `ss.death` evidently show there. How Starsim routes them wasn't checked.
  - **Not a pitfall:** they come to 0.005 per person-year on a curing course at both steps, a plausible background death rate for that population.
- **Run time:**
  - **The compute:** a monthly agent run takes about 40 s (timed with the laptop kept awake).
  - **The wall clock:** the stage took 41 hours, because the laptop slept repeatedly during the run. The run times in `tb_summary.log` include that sleep.

## Stage 4: births' female share in `fastode_age`, and exp 32

- **`fastode_age`** gains `Demography(birth_female=)`: `'start'` (the 1950 population's overall share, as before), `'young'` (its share at ages 0–4, Starsim's default from rc5), or a number. With `'young'`, `test_copy_matches_v4_ode` passes against v4's ODE with Starsim's default births.
- **Exp 32's targets,** 'young' ÷ 'start', at its best point and 20 posterior draws (`birth_female_targets.csv`):

  | Target | Range |
  |---|---|
  | Prevalence 15+, 2020 | 1.0019–1.0043 |
  | Prevalence 15+ ÷ incidence | 0.9977–0.9985 |
  | Prevalence 2040 ÷ 2020 | 1.0021–1.0040 |
  | Incidence and mortality 2024 ÷ 2015 | 1.0009–1.0018 |
  | Others | within 0.15% |
  | Population, 2020 | 1.0100 |

  The log-likelihood changes by −0.04 to +0.10. Every target moves by under 0.5%, against the 1% threshold and σ of about 10%, so **exp 32's posterior stands**.

## What this means

- **The v4 pins can move** to starsim `0089ac0f` and tbsim `a2ae583`. The only model change is births' female share, which follows Starsim's new default once `fastode_age` uses `'young'`.
- **The workarounds can go** for new work: the shifted fertility table, `PregnancyStart`, `joint_sex()`, births' share set by hand, counting losses by balance, and adding 1e-3 to ages. `tbode4.demog` keeps them for rerunning exps 29–32 at their pins.
- **Still needed:** `tbode4.use_drawn_start()` (proposal 07, unchanged upstream).
- **Proposal 06** is resolved upstream.
- **Still open:** WPP's sex ratio at birth. The ages 0–4 fraction stands in for it.

## Files

- **Environment:** `pyproject.toml`, `uv.lock` (`.venv` not committed).
- **Code:** `run.py`: `tests`, `pins_ode`, `pins_compare`, `demog`, `demog_summary`, `tb`, `tb_summary`, `fastode_age`.
- **Stage 1:** `tests_starsim_4.0.0rc5.log`, `pins_compare.csv`, `pins_compare.log`, `pins_ode_starsim_*.obj`.
- **Stage 2:** `demog.obj`, `demog.out`, `demog_summary.log`, `demog_ratios.csv`, `demog_shares.csv`, `demog_fertility.csv`.
- **Stage 3:** `tb_ode.obj`, `tb_runs.obj`, `tb.out`, `tb_summary.log`, `tb_by_window.csv`, `tb_pitfall_checks.csv`.
- **Stage 4:** `fastode_age.log`, `birth_female_targets.csv`.
- **Raw runs** (`*.obj`) are not committed.
