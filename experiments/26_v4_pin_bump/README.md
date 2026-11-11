# 26: Does treatment still match between agents and the ODE after moving the v4 pins to Starsim's network fix?

**Status:** done; see [SUMMARY.md](SUMMARY.md). Added after the first run: `fail_death.py` and `timeseries_check.py` (TB deaths on treatment read about 1% low) and `extra_reps.py` (one arm read low). The first run also showed that Starsim now indexes a step's flows by its start, so `tbode4.ode_incidence` and exp 23's `yearly()` were fixed and every run repeated.
**Track:** v4, at the v4 track's new pins (`v4/pyproject.toml`): starsim `origin/v4` @ `d83f4f87` (4.0.0rc4, which includes `6dcad9fd`, Cliff's network fix) and tbsim `origin/v4` @ `619c296`. These are [exp 25](../25_starsim_network_fix/SUMMARY.md)'s pins. The previous pins (starsim `80c81f53`, tbsim `dcfb522`) are those of exps 12–25.

## Background

[Exp 25](../25_starsim_network_fix/SUMMARY.md) found that at these pins, agents on `ss.RandomNet(dur=0)` match Starsim's ODE of the network, reinfection included, and that homogeneous agents are unchanged. It didn't check treatment, because `tbode4.TBTreat` didn't build at these pins.

`TBTreat` is [exp 23](../23_treatment/SUMMARY.md)'s treatment design (CONVENTIONS.md 9): a `tbsim.TB` subclass that adds treatment as declared states and transitions, so the same model runs as v4's ODE or as agents. To build at the new pins it needed these changes, all to how the model is declared, not to what it does:

- `ss.State` takes the label as its first argument, so the state's name is passed as `name=`;
- a transition with several destinations takes them as a dict of each destination and its probability, not a list with `p=`;
- births that replace TB deaths come from `ss.unborn`, not `ss.null`.

With these, and with `ss.timeseries` taking `interpolation=` rather than `mode=`, all of `v4/tests/` passes at the new pins: the anchor (v4's ODE reproduces tbsim's `TB_ODE` and our `TBODE`), `fastode`, the latent-structure calculations, and `TBTreat`'s ODE against `fastode`.

## Question

Do agents still match v4's ODE with treatment, and does treatment still have none of exp 23's pitfalls, at the new pins?

## Design

Exp 23's `run.py` and `diagnose.py`, with `TB_EXP23_FOLDER` sending their outputs here, and flows dated by Starsim's new indexing (`flows_at_start`; SUMMARY finding 1):

- **Parameters:** exp 21's set 1 natural history, with tb_LAI_TPT's best-fit treatment values. Homogeneous mixing, no network.
- **Scenarios:** 2024's treatment held constant (`eq`, compared over 2005–2049), and tb_LAI_TPT's scale-up from 1995 (`ramp`, 1980–2035, compared over 1995–2035).
- **Arms:** a fixed 6-month course and a rate, at monthly and weekly steps; 200,000 agents and 40 reps each, with exp 23's seeds; the ODE alongside.
- **Outputs:** agents ÷ ODE for prevalence, and per prevalent case for incidence (first entry to asymptomatic TB per episode, plus relapses), notifications, the number on treatment and TB mortality; exp 23's pitfall checks; the ODE's fixed course against a rate.

**Pass:** agents ÷ ODE within the same bounds as exp 23 (per prevalent case within 1%; prevalence within about 2 SE per arm), every pitfall check passing, and the ODE's equilibria close to exp 23's (any change there comes from Starsim, not from agents).

Run from the repo root:

```bash
TB_EXP23_FOLDER=experiments/26_v4_pin_bump PYTHONPATH=v4 v4/.venv/bin/python experiments/23_treatment/run.py
TB_EXP23_FOLDER=experiments/26_v4_pin_bump PYTHONPATH=v4 v4/.venv/bin/python experiments/23_treatment/diagnose.py
PYTHONPATH=v4 v4/.venv/bin/python experiments/26_v4_pin_bump/extra_reps.py   # 60 more reps (new seeds) of the equilibrium, fixed-course, weekly arm
PYTHONPATH=v4 v4/.venv/bin/python experiments/26_v4_pin_bump/fail_death.py   # the share of failing courses that end in TB death, 8 × 400,000 agents
PYTHONPATH=v4 v4/.venv/bin/python experiments/26_v4_pin_bump/timeseries_check.py   # TB deaths on treatment, time-varying against constant treatment
```
