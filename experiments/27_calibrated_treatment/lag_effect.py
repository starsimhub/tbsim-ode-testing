"""
Exp 27: how much of agents' monthly drift in prevalence could the half-step timing of time-varying inputs explain?

Agents apply a time-varying rate as if half a step earlier than v4's ODE (ramp_timing.py). This runs v4's ODE at each set with every
time-varying input (coverage, efficacy, sym_dead) delayed by half a month, ÷ the same ODE without the delay, monthly steps: prevalence,
and notifications per prevalent case, by year. Run from the repo root:

    PYTHONPATH=v4 v4/.venv/bin/python experiments/27_calibrated_treatment/lag_effect.py
"""
import os
os.environ.setdefault('NUMBA_NUM_THREADS', '1')
import sys
import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss

outdir = sc.thispath(__file__)
sys.path.insert(0, str(outdir))
import run as r27 # noqa: E402
tbode4 = r27.tbode4


def ode_run(s, lag):
    """ v4's ODE for one set, monthly, with the time-varying inputs delayed by `lag` years """
    init, non0, _ = r27.starting_state(s)
    t = s.trend
    lai = tbode4.lai_tpt_treatment(row=dict(tx_covg_2010=t.tx_covg_2010, tx_covg_2024=t.tx_covg_2024, efficacy_1995=t.efficacy_1995,
                                            efficacy_2024=t.efficacy_2024, p_relapse=t.p_relapse), end=r27.stop)
    grid = np.concatenate([[r27.cb.t0], np.round(np.arange(1995, r27.stop + 1/24, 1/12), 6)])
    tx_rate = ss.timeseries({float(g + lag): ss.peryear(float(tbode4.theta_from_coverage(lai.coverage(g)))) for g in grid}, interpolation='linear')
    tx_success = ss.timeseries({float(g + lag): float(lai.efficacy(g)) for g in grid}, interpolation='linear')
    sym_dead = ss.peryear(ss.timeseries({1995.0 + lag: t.sym_dead_1995, 2024.0 + lag: t.sym_dead_2024}, interpolation='linear'))
    nat = {k: v for k, v in s.nat.items() if k != 'sym_dead'}
    tb = tbode4.TBTreat(init=init, course='fixed', tx_rate=tx_rate, tx_success=tx_success, tx_p_relapse=t.p_relapse, dur_relapse_tx=1.5,
                        sym_dead=sym_dead, **nat)
    sim = ss.Sim(tb, solver='ode', start=r27.start, stop=r27.stop, dt=ss.months(1), total_pop=1, verbose=0).run()
    out = r27.per_step(sim, tbode4.ode_incidence(sim, non0=non0).first)
    return r27.yearly(out, sim.t.dt_year)


if __name__ == '__main__':
    sets = sc.load(r27.fn.sets)
    rows = []
    for i, s in enumerate(sets):
        base, lagged = ode_run(s, 0.0), ode_run(s, 1/24)
        d = pd.DataFrame(dict(set=i, year=base.index, prevalence=lagged.prevalence/base.prevalence,
                              notifications_per_case=(lagged.notifications/lagged.prevalence)/(base.notifications/base.prevalence)))
        rows.append(d)
    df = pd.concat(rows, ignore_index=True)
    df.to_csv(outdir / 'lag_effect.csv', index=False)
    years = [1995, 1996, 1997, 2000, 2005, 2010, 2015, 2020, 2024, 2030, 2040]
    with pd.option_context('display.width', 200, 'display.precision', 4):
        print('ODE with inputs delayed by half a month ÷ ODE, monthly; mean over the 12 sets (min–max):')
        g = df[df.year.isin(years)].groupby('year')
        print(pd.DataFrame({'prevalence': g.prevalence.mean(), 'prev min': g.prevalence.min(), 'prev max': g.prevalence.max(),
                            'notif per case': g.notifications_per_case.mean()}).to_string())
