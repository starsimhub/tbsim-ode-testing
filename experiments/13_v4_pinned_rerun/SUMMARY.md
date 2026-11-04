# 13 Cliff's v4 redo at pinned versions and full reps: summary

**Track: v4.** tbsim `origin/v4` @ `dcfb522`, starsim `origin/v4` @ `80c81f53` (`versions_v4.json`).

On 2026-10-27, Cliff reran the key findings of exps 01–11 under tbsim v4 ([`v4/v4_findings.qmd`](../../v4/v4_findings.qmd)). His versions had uncommitted changes, and his runs had 12–40 reps and no daily step. This experiment reran his code unchanged at our pinned versions: 40 reps per arm, a daily step for set A, and 200 reps for the three arms that were about 2 SE from the ODE after 40.

**How to read the tables.** Agents ÷ ODE is each rep's mean over 1950–2000, divided by the ODE's; 1.00 = match. Values are mean ± SE over reps. Agents are 50,000, started at the ODE's equilibrium. **Every run uses homogeneous mixing (v4's built-in mass action, no network), except the two network arms in finding 5.** So the ODE's `beta` is used as it is, and the results outside finding 5 don't depend on how `beta` is converted for a network. `beta` is calibrated to 500/100k at each arm's ODE equilibrium, with births replacing every death.

The arms:
- **A0:** tbsim's default rates, near threshold.
- **B0:** tb_LAI_TPT's best-fit natural history.
- **Ak, Bk:** the same with `k` = 5.846 as exact Poisson-mixture classes (`TBClasses`).
- **Bremote:** set B with v4's `k_latent` = 5.846 and no progression from the remote state.

## Findings

**1. Cliff's results reproduce exactly at the pinned versions.** Every rep with one of his seeds gives his result. In all 10 arms he ran, the first 12–40 reps match his means to three decimals (e.g. A0 monthly 1.043, Bk weekly 0.979), and the one-step and growth results are identical. v4's draws are keyed by seed, so the code paths these runs use didn't change between his versions and ours. The calibrated `beta`s also match his (A0 9.926, B0 1.952, Ak 36.898, Bk 24.799, Bremote 25.115).

**2. The baseline matches the ODE at every step, including near threshold.** Prevalence, agents ÷ ODE:

| Arm | Daily | Weekly | Monthly |
|---|---|---|---|
| A0, v4 agents | 0.99 ± 0.02 (40) | 1.01 ± 0.02 (40) | **1.000 ± 0.009** (200) |
| B0, v4 agents | — | 0.989 ± 0.006 (200) | 0.99 ± 0.02 (40) |
| A0, v0.10 `tbsim.TB` (exp 03) | 0.95 | 0.84 | 0.67 |
| A0, v0.10 `TBExact` (exps 05, 10) | 1.00 | 1.00 | 0.97–0.98 |

Reps in brackets.
- **The monthly +4% was chance.** Cliff's A0 monthly result (1.043 ± 0.019) came from the first 40 seeds; seeds 40–199 give 0.990 ± 0.010, and all 200 give 1.000 ± 0.009.
- **Weekly B0** is within about 1% (0.989 ± 0.006, z = −1.8). Its first 40 reps were 0.972 ± 0.013.

**3. Progression by time since infection matches the ODE, with either structure.** Prevalence, agents ÷ ODE:

| Arm | Weekly | Monthly |
|---|---|---|
| Ak (set A, classes) | — | 1.00 ± 0.02 (40) |
| Bk (set B, classes) | 0.993 ± 0.004 (200) | 1.003 ± 0.008 (40) |
| Bremote (set B, `k_latent`) | 1.01 ± 0.01 (40) | 1.01 ± 0.01 (40) |
| Bk, v0.10 `TBExact` (exp 08) | 1.00 | 1.01 ± 0.01 |
| Bk, v0.10 `tbsim.TB` (exp 08) | 1.18 | 1.83 |

Weekly Bk's first 40 reps were 0.979 ± 0.010; with 200 it is within 1% (z = −1.7). Progression targets from the ODE are unchanged from Cliff's: Bk 0.0511 / 0.674 and Bremote 0.0489 / 0.676, against 0.0523 / 0.663.

**4. Other outputs.** Agents ÷ ODE, same reps as above:

| Arm | Latent infection | No longer infected | Incidence per prevalent case | TB deaths |
|---|---|---|---|---|
| A0 daily | 0.99 ± 0.02 | 0.996 ± 0.004 | 1.000 ± 0.002 | 0.99 ± 0.02 |
| A0 weekly | 1.03 ± 0.02 | 1.001 ± 0.003 | 0.997 ± 0.002 | 1.01 ± 0.02 |
| A0 monthly | **1.089 ± 0.009** | 0.995 ± 0.002 | 0.999 ± 0.001 | 0.999 ± 0.009 |
| B0 weekly | 1.002 ± 0.006 | 0.996 ± 0.003 | 1.000 ± 0.001 | 0.988 ± 0.006 |
| B0 monthly | 1.04 ± 0.02 | 1.000 ± 0.005 | 0.999 ± 0.002 | 1.00 ± 0.02 |
| Bk weekly | 1.001 ± 0.003 | 0.999 ± 0.000 | 1.001 ± 0.001 | 0.995 ± 0.004 |
| Bk monthly | 1.027 ± 0.006 | 0.998 ± 0.000 | 1.001 ± 0.002 | 1.004 ± 0.008 |
| Bremote monthly | 1.033 ± 0.007 | 0.997 ± 0.001 | 0.998 ± 0.002 | 1.004 ± 0.009 |

- **Latent infection reads high in end-of-step snapshots at monthly `dt`:** +9% in A0 (Cliff, 40 reps: +13%) and +3–4% in set B. The excess is roughly halved at weekly `dt` and gone at daily. As Cliff noted, infections are dated at the end of their step, so a snapshot includes that step's new infections in full. It matters only if a target is a snapshot of latent infection.
- **Incidence counted from transitions** is within 0.3% of the ODE per prevalent case in every arm, so exp 11's snapshot problem doesn't arise in v4.

**5. On a network, `beta/K` per contact gives 0.74 of the ODE's prevalence; the matched `beta_edge` gives 1.03.** These are the experiment's only network runs: A0 monthly on `ss.RandomNet(n_contacts=10, dur=0)`, 40 reps each. The two arms differ only in the per-contact rate:
- **`beta/K`:** starsim's suggested conversion. It gives 0.737 ± 0.017 (Cliff, 20 reps: 0.780).
- **Matched:** −ln(1 − `beta`·`dt`/K)/`dt`, the same formula as `beta_edge_from_ode`. Its per-step probability equals the ODE's `beta`·`dt`/K. It gives 1.03 ± 0.02. This is the per-contact step effect of exp 12 (4% of effective `beta`), magnified near threshold.

**6. One step matches `expm(Q·dt)`, and growth lags by 0.1–0.3 years (both identical to Cliff's).**
- **One step:** from each state, one million agents × 2 reps. The largest difference from `expm(Q·dt)` is 0.0005; the ODE's is 6e-6.
- **Paths agents don't take:** the only clear misses (|z| = 2.9–5.9) are five paths needing three moves in a month, with exact probabilities of 4e-6 to 2e-5, which agents never take.
- **Growth:** from a 0.1% seed (Bk, 8 reps), agents reach 250/100k and 400/100k 0.09–0.27 years after the ODE (± 0.09–0.15). This is about the size of the generation-interval lag of end-of-step dating, but 8 reps can't separate monthly from weekly.

**7. Run time: 1.3–2.4 ms per step at 50,000 agents** (A0 and Bk, daily to monthly; best of 3, one process). Cliff measured 3.8–6.1 ms on his laptop. The ODE takes 0.01–0.2 s for 100 years.

## What this means

- **v4 agents are consistent with the v4 ODE** for the baseline and for progression by time since infection, at daily, weekly and monthly `dt`: within 1%, at 200 reps where the question was close. v0.10's step-ordering and start-of-step-hazard gaps are gone, without `TBExact`.
- **What remains is small and understood:**
  - end-of-step latent snapshots at monthly `dt` (+3 to +9%);
  - a growth lag of a few months;
  - the per-contact `beta` on networks, which has a fix (the matched `beta_edge`) and a proposal ([`proposals/01_network_beta/`](../../proposals/01_network_beta/proposal.md)).
- **Both latent structures run exactly under both solvers.** So choosing between them is a modelling question; see [`proposals/02_latent_progression/options.md`](../../proposals/02_latent_progression/options.md).
- **Seeds:** with v4, a rep with the same seed and unchanged code gives the same result. New reps need new seeds to add information.

## Files

- `run.py`: the design. Rerunning resumes, and skips arms already saved.
- `run_first40.log`: the first pass with 40 reps everywhere. `run_more.log`: seeds 40–199 for A0 monthly, B0 weekly and Bk weekly.
- `summary.csv`: agents ÷ ODE for every arm and output: mean, SE, and the 10th and 90th percentiles over reps.
- `results/`: Cliff's scripts' outputs (`calibration.json`, `progression.json`, `equilibrium.obj` with per-rep yearly series and the ODE references, `onestep.csv`, `growth.obj`, `timing.csv`).
- `versions_v4.json`.
