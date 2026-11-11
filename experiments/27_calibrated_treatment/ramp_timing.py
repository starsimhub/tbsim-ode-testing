"""
Exp 27: when a rate varies in time (ss.timeseries), do agents and v4's ODE apply it at the same point in each step?

People leave Y at a rate rising linearly from 0 at the end of 1994 (as treatment coverage does), and the flow per person in Y each
year is compared, agents ÷ ODE, at monthly and weekly steps; also a falling rate (as sym_dead does). If agents read the rate at the
start of each step and the ODE over the step, agents fall short by about (step/2) × (rate's slope ÷ rate). Run from the repo root:

    PYTHONPATH=v4 v4/.venv/bin/python experiments/27_calibrated_treatment/ramp_timing.py
"""
import os
os.environ.setdefault('NUMBA_NUM_THREADS', '1')
import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss
import tbode4

outdir = sc.thispath(__file__)
t0 = 1994 + 364/365
n_agents, seeds = 1_000_000, range(4)
shapes = dict(rising={t0: 0.0, 2010.0: 0.5}, falling={1995.0: 0.5, 2024.0: 0.1})


def run(shape, step, solver, seed=0):
    dt = dict(month=ss.months(1), week=ss.weeks(1))[step]
    Y, Z = ss.states(2)
    rate = ss.peryear(ss.timeseries(shapes[shape], interpolation='linear'))
    dis = ss.Disease(name='d', Y=Y, Z=Z, init=dict(Y=1.0), go=ss.Transition(Y, Z, rate=rate), back=ss.Transition(Z, Y, rate=ss.peryear(5.0)))
    kw = dict(n_agents=n_agents, rand_seed=seed) if solver == 'abm' else dict(total_pop=1)
    sim = ss.Sim(dis, solver=solver, start=1990, stop=2026, dt=dt, verbose=0, **kw).run()
    r = sim.results.d
    year = np.asarray(sim.t.yearvec, dtype=float)
    y, f = np.asarray(r['Y'].values, dtype=float), np.asarray(r['go'].values, dtype=float)
    per = f/y # Flow per person in Y, each step (flows indexed by the step's start)
    df = pd.DataFrame(dict(year=np.floor(year + 1e-9).astype(int), per=per)).groupby('year').per.sum()
    return {**dict(shape=shape, dt=step, solver=solver, seed=seed), **{int(k): float(v) for k, v in df.items()}}


if __name__ == '__main__':
    jobs = [dict(shape=sh, step=st, solver='ode') for sh in shapes for st in ['month', 'week']]
    jobs += [dict(shape=sh, step=st, solver='abm', seed=s) for sh in shapes for st in ['month', 'week'] for s in seeds]
    df = pd.DataFrame(sc.parallelize(run, iterkwargs=jobs, ncpus=tbode4.ncpus))
    df.to_csv(outdir / 'ramp_timing_runs.csv', index=False)
    years = [1995, 1996, 1997, 1998, 2000, 2005, 2009, 2015, 2023]
    rows = []
    for (sh, st), g in df.groupby(['shape', 'dt']):
        o = g[g.solver == 'ode'].iloc[0]
        a = g[g.solver == 'abm']
        rows.append({**dict(shape=sh, dt=st), **{y: (a[y]/o[y]).mean() for y in years}})
    out = pd.DataFrame(rows)
    out.to_csv(outdir / 'ramp_timing.csv', index=False)
    with pd.option_context('display.width', 200, 'display.precision', 4):
        print(f'Flow per person, agents ÷ ODE by year ({n_agents:,} agents × {len(seeds)} seeds):'); print(out.to_string(index=False))
