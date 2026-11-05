# 13 Cliff's v4 redo at pinned versions and full reps: do its findings hold?

**Track: v4.** tbsim `origin/v4` @ `dcfb522`, starsim `origin/v4` @ `80c81f53` (the pins in `v4/pyproject.toml`; exact versions in `versions_v4.json`).

## Question

On 2026-10-27, Cliff reran the key findings of exps 01–11 under tbsim v4 ([`v4/v4_findings.qmd`](../../v4/v4_findings.qmd), scripts at the top of `v4/`). He found:

- v4's agents match the v4 ODE without anything like `TBExact`;
- progression by time since infection matches across solvers, both with Poisson-mixture latent classes declared in v4 (`TBClasses`) and with v4's own `k_latent`;
- a few small timing effects remain.

Two things limit how far those numbers can be relied on:

- **Versions.** The runs used tbsim `v4` @ `7b2e644` and starsim `v4-oct27` @ `fe6c6f9`, both with uncommitted changes, so they can't be traced to a commit.
- **Reps.** They were cut for time to 12–40 reps, with no daily step. That leaves SEs of 1–3.5%. Set A at monthly `dt` came out at 1.04 ± 0.02 of the ODE, which can't rule out a bias of ~0.5% in effective `beta` (~4% of prevalence near threshold).

**Do Cliff's findings hold at our pinned versions, with 40 reps per arm and a daily step for set A?**

## Design

Cliff's code is used as it is: models, calibration, outputs and analysis (`v4/common.py`, `calibrate.py`, `run_*.py`, `analyze.py`). Two changes were made to `common.py`, neither of which changes the model or any result:

- `results/` can be redirected with the environment variable `TB_V4_FOLDER` (default: `v4/`, as before), so that this experiment's results go to its own folder. Worker processes inherit the variable.
- `versions()` falls back to `tbode4.versions` when the paths of Cliff's checkouts don't exist.

**Calibration.** Each arm's ODE `beta` is recalibrated to 500/100k active (ASY + SYM) prevalence at the v4 ODE's equilibrium, at the pinned versions (`calibrate.py`, unchanged). Births at `mu` = 1/70 replace every death, including TB deaths. That differs from exps 01–12, which didn't replace TB deaths, so `beta` differs slightly from exp 12's (e.g. A0 9.93 vs 9.917).

**Arms** (from Cliff's `calibrate.py`):

- **A0:** tbsim's default rates, constant progression; near threshold.
- **B0:** tb_LAI_TPT's best-fit natural history, reinfection risk ratio π = 0.24.
- **Ak, Bk:** the same with `k` = 5.846, as exact Poisson-mixture classes (`TBClasses`).
- **Bremote:** B with v4's `k_latent` = 5.846 and no progression from the remote state.

**Equilibrium runs** (`run_equilibrium.py`'s `run_one`, unchanged):

- 50,000 agents, started at the ODE equilibrium (`init=` its shares), 1900–2000;
- homogeneous mixing unless marked;
- each run reduced to yearly means. Agents ÷ ODE is the mean over 1950–2000 of each run, divided by the ODE's.

Tables give the mean ± SE over reps, with the 10–90% range across reps.

| Arm | `dt` | Transmission | Reps here | Cliff's reps |
|---|---|---|---|---|
| A0 | day | homogeneous | 40 | — (new) |
| A0 | week | homogeneous | 40 | 20 |
| A0 | month | homogeneous | 200 | 40 |
| A0 | month | `RandomNet(10)`, `beta/K` per contact | 40 | 20 |
| A0 | month | `RandomNet(10)`, matched `beta_edge` | 40 | 20 |
| B0 | week | homogeneous | 200 | 12 |
| B0 | month | homogeneous | 40 | 20 |
| Ak | month | homogeneous | 40 | 16 |
| Bk | week | homogeneous | 200 | 12 |
| Bk | month | homogeneous | 40 | 20 |
| Bremote | week | homogeneous | 40 | 12 |
| Bremote | month | homogeneous | 40 | 20 |

**About the seeds and the 200-rep arms.** Seeds run from 0, as in Cliff's runs. Since v4's draws are keyed by seed, step and round, a rep with one of Cliff's seeds repeats his run whenever the code paths it uses are unchanged between his versions and ours. So the first 12–40 reps of each arm test the versions, not chance. Three arms came out about 2 SE from the ODE after 40 reps: A0 monthly (+4%, also Cliff's result) and B0 and Bk weekly (−2 to −3%). For these, seeds 40–199 were added, to resolve differences of ~1.5% in prevalence.

**Unchanged from Cliff's design and rerun as is, at the pinned versions:**

- **One step against `expm(Q·dt)`** (`run_onestep.py`): one million agents, 2 reps, sets A and B, monthly and weekly.
- **Growth from a seed** (`run_growth.py`): Bk, 100k agents, 8 reps, monthly and weekly.
- **Run time** (`run_timing.py`): run last and alone, so other runs don't slow it.

**Comparisons.** Each result is set beside Cliff's (`v4/results/summary.md`), and the baseline beside v0.10's (`tbsim.TB` 0.95 / 0.84 / 0.67, `TBExact` 1.00 / 1.00 / 0.97 at daily / weekly / monthly `dt`; exps 03, 05).

## How to run

```bash
v4/.venv/bin/python experiments/13_v4_pinned_rerun/run.py        # from the repo root; about 1–1.5 hours on 4 cores
```

Outputs go to `results/` in this folder: `calibration.json`, `progression.json`, `equilibrium.obj`, `onestep.csv`, `growth.obj`, `timing.csv`. The script also writes `summary.csv` (agents ÷ ODE for every arm and output), `run.log` and `versions_v4.json`.
