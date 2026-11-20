# 33: Moving the v4 pins to Starsim 4.0.0rc5: does Starsim's own demography now match without our workarounds?

**Status:** designed, not yet run. Runs after exp 32 finishes, so that exp 32 runs entirely at the current pins.
**Track:** v4, moving from starsim `origin/v4` @ `d83f4f87` (4.0.0rc4) and tbsim `origin/v4` @ `619c296` to starsim `origin/v4` @ `0089ac0f` (4.0.0rc5; the merge of Starsim PR #1625) and tbsim `origin/v4` @ `a2ae583` ("fixes for rc5"), if that pairing builds (stage 1).
**Run on:** the laptop (stages 1, 2 and 4) and zebra (stage 3).

## Why this experiment

Exps 29 and 30 matched v4's agents to its ODE with India's WPP demography, but only with workarounds for behaviours of Starsim's ([proposal 06](../../proposals/06_pregnancy/proposal.md), sent to Starsim on 2026-11-19; CONVENTIONS.md 10):

| Starsim behaviour at `d83f4f87` (proposal 06 change) | Our workaround (`tbode4.demog`) |
|---|---|
| `ss.Pregnancy` applies fertility rates at the age at conception, not at delivery as WPP gives them (1) | A fertility table shifted younger by the mean gestation, 0.769 years (`shifted_table()`) |
| Its start seeds too few pregnancies: 1 − e^(−h·d) of women, not the steady state h·d (2) | `PregnancyStart`, overriding `burnin_plan()` |
| Its miscarriage and stillbirth results miss pregnancies lost when the mother dies (3) | Counted by balance |
| A false warning that births are made twice (4) | Ignored |
| The compartmental solver starts each age group half female (5) | A joint split by age and sex (`joint_sex()`); births then follow the ODE's starting female share (48.63%), so agents' `ss.Pregnancy` is given the same `p_female` |
| Agents' ages are 32-bit floats, putting ±1.7% bands in results by age group (6) | Ages counted in our own analyzer with 1e-3 years added |

Starsim's `origin/v4` has since taken up every one of these, mostly in PR #1625 (merged 2026-11-20, `0089ac0f`). From its code and tests:

1. **Fertility:** `ss.Pregnancy(fertility_age='delivery')`, the default, shifts the rates and `min_age`/`max_age` by the mean gestation.
2. **The start:** `burnin_plan()` uses h·d with the unrounded mean gestation.
3. **Losses at the mother's death:** counted as miscarriages or stillbirths.
4. **The warning:** it now appears only if another module also makes births.
5. **Sex by age:** `ss.population()` gives the fraction female by age. Agents take it at their age, and a compartmental model split by age and sex starts with it in each age group. **Births in both solvers default to the fraction female at ages 0–4** of the named population, unless the birth sets `p_female`.
6. **Ages:** 64-bit, and each step adds 1e-12 years more than the step, so ages reach whole years on time.

Proposal 07 (agents' starting states by age group; drafted, not sent) is unchanged: counts by group are still rounded to whole agents, and shares read per agent are still drawn one state at a time. So `tbode4.use_drawn_start()` stays.

Moving the pins is a change in its own right: there are 82 Starsim commits between the two pins, several of them refactors. As at exp 26's bump, `v4/tests/` must pass and results must be indexed as before (CONVENTIONS.md 7).

## Question

1. Does tbsim `a2ae583` build and pass `v4/tests/` with starsim `0089ac0f`, with flows and states indexed as before?
2. With Starsim's defaults and none of our demography workarounds, do agents match the ODE's demography as closely as exp 29's corrected setup did?
3. With TB and treatment, do agents still match the ODE as in exp 30?
4. Births' female share changes from the ODE's starting share (48.63% in exp 30 and `fastode_age`) to the fraction at ages 0–4. How much does that move the calibration targets, and so exp 32's posterior?

## Design

### Stage 1: the bump

- **A trial environment first:** this folder gets its own `pyproject.toml`, `uv.lock` and `.venv` with the new pins, as exp 18 had. `v4/` moves to the new pins only after stages 1–3 pass. Until then, exps 12–32 rerun unchanged in `v4/.venv`.
- **`v4/tests/`** run in the trial environment (with `PYTHONPATH=v4`): the anchor (v4's ODE reproduces `TB_ODE` and `TBODE`), `fastode`, the latent calculations, `TBTreat` against `fastode`, the starting states, and `fastode_age` against v4's ODE.
  - `tbode4` changes needed only to build (Starsim's API) are made in `tbode4` and listed in the SUMMARY.
  - `fastode_age`'s check against v4's ODE runs with births' female share set explicitly to 48.63%, the current value, so it tests the copy rather than the change of default (stage 4 treats that).
- **Indexing:** flows at index i are the step from t_i (CONVENTIONS.md 7), checked as exp 26 did, on exp 23's ODE run.
- **If the pairing doesn't build:** stop, record what fails, and wait for a matching tbsim commit.

### Stage 2: demography alone, Starsim's defaults against exp 29's corrected setup

Exp 29's setting: India, WPP 2024 rates as bundled with Starsim, 1950–2050, monthly steps, deaths by age and sex (`ss.death`) in both solvers, no TB.

| Arm | Solver | Setup | What it tests |
|---|---|---|---|
| Q, E | ODE, quarter- and eighth-year age groups | `ss.population('India', 1950)`, `stratify(age=, sex=)` with Starsim's default split by sex, `ss.birth` by age and sex with the default `p_female` | The reference at the new pins; E extrapolated with Q to exact aging (2 × E − Q), as in exp 29 |
| Qj | ODE, quarter-year groups | As Q, with `joint_sex()` and `p_female` = 48.63% | Whether the new default split by sex at the start is exp 29's `joint_sex()` (Qj = Q, apart from births' female share), and against exp 29's saved R4x (Starsim's other changes) |
| D | Agents | `ss.population('India', 1950)`; `ss.Pregnancy` with WPP's table as given and every default (`fertility_age='delivery'`, `burnin=True`, `p_female`) | The workarounds are no longer needed |
| W | Agents | Exp 29's P2b (the shifted table, `PregnancyStart`, `p_female` = 48.63%), with `fertility_age='conception'` so the table isn't shifted twice | Against exp 29's saved P2b: Starsim's other changes to agents |

- **Agents:** 10 reps of 100,000 at the start, new seeds.
- **Measured:** as exp 29. Population ÷ 1950 every 10 years; the share in 0–14, 15–29, 30–49, 50–64 and 65+; crude birth and death rates; realized fertility by age at delivery ÷ WPP's and total fertility, by decade; births' female share in each solver; mean ± SE across reps.
- **And the other fixes:**
  - **Losses:** miscarriage + stillbirth results against exp 29's balance count, every year.
  - **The warning:** whether it appears.
  - **Ages:** results by 5-year group from Starsim's own results against our analyzer with 1e-3 years added. The two should now agree.

### Stage 3: TB with demography, a reduced rerun of exp 30

Exp 30's model and starts (`tbode4.TBTreat`; sets 0 and 6; the settled starting table by single year and sex; agents drawn with `use_drawn_start()`; 1950–2041). Demography as stage 2's D (agents) and Q and E (the ODE), with no workarounds. Births' female share follows the new default in both solvers.

- **Agents:** 20 reps per set at weekly steps (exp 30's criteria were set at weekly steps) and 20 at monthly steps; 100,000 agents; new seeds.
- **Measured:** exp 30's agents ÷ ODE (extrapolated to exact aging): all-ages prevalence, prevalence in the five age groups, ratios per prevalent case (incidence, notifications, the number on treatment, TB mortality), population size and age shares, over 1950–94, 1995–2024 and 2025–40; mean ± SE. Side by side with exp 30's saved ratios.
- **Expected effects, reported but not counted as mismatches:** as exp 30 (recent latent infection high in monthly snapshots; treatment starts low in 1995 at monthly steps; monthly agents' slow drift from the ODE's start).

### Stage 4: births' female share in `fastode_age` and exp 32

- **`fastode_age`:** births take the fraction female at ages 0–4 in WPP's 1950 population, as Starsim's ODE now does. `v4/tests/test_fastode_age.py` then runs against v4's ODE with the default.
- **Effect on the calibration:** exp 32's targets at its posterior median and at 20 posterior draws, with each share, as a ratio. If any target moves by more than about 1% (against σ of about 10%), exp 32 is rerun at the new pins. Otherwise its posterior stands, with the difference noted.

## What would answer the question

- **The bump passes** if every `v4/tests/` check passes and indexing is unchanged (stage 1). The ODE's equilibria and exp 23's run should agree with the current pins to about 1e-7: Starsim's ODE solver isn't meant to have changed.
- **The workarounds can go** if D is within about 0.1% of the exact-aging ODE's population through 2050, and within about 1% of its share in each broad age group, as P2b was in exp 29 (stage 2). W should equal exp 29's P2b within noise; if it doesn't, something else changed in agents.
- **TB still matches** if stage 3 meets exp 30's criteria: agents ÷ ODE within about 2 SE and 1% for all-ages prevalence, each ratio per prevalent case and population size in the two later windows at weekly steps, and prevalence in each broad age group within about 2%.
- **Then:**
  - **The repo:** move `v4/`'s pins; drop the workarounds from CONVENTIONS.md 10 (`tbode4.demog` keeps them for rerunning exps 29–32); update CLAUDE.md's demography facts and the README; mark proposal 06 resolved upstream.
  - **Still open:** WPP's sex ratio at birth. The fraction at ages 0–4 stands in for it.

## Not covered

- Proposal 07's starting-state behaviours, unchanged at `0089ac0f`.
- Starsim's other new features since `d83f4f87`, beyond what `v4/tests/` and these runs touch.
- Child transmission (the next model row).

## Files

- `pyproject.toml`, `uv.lock`: the trial environment (`.venv`, not committed).
- `run.py`: `tests` (stage 1), `demog` (stage 2), `tb` (stage 3), `fastode_age` (stage 4), `summary`. Every run's raw outputs saved with `sc.save`, with `versions.json`.
- `SUMMARY.md`: findings, written after running.
