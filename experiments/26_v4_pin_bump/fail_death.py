"""
Exp 26: the share of failing courses that end in TB death, agents and ODE against its exact value, at the new pins.
Exp 23's checks/fail_death.py with more agents and several seeds, since exp 26's agents read TB deaths on treatment about 1% low.
People enter treatment from Y at a constant rate, face TB death on treatment, and leave after a fixed (or exponential) course. Run from the repo root:

    PYTHONPATH=v4 v4/.venv/bin/python experiments/26_v4_pin_bump/fail_death.py
"""
import os
os.environ.setdefault('NUMBA_NUM_THREADS', '1')
import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss
import tbode4

outdir = sc.thispath(__file__)
d_rate, D = 0.3959, 0.5 # Exp 21 set 1's sym_dead; course length (years)
n_agents, seeds = 400_000, range(8)


def run(step, fixed, solver, seed=0):
    dt = dict(month=ss.months(1), week=ss.weeks(1))[step]
    Y, Tf, Z = ss.states(3)
    dur = ss.constant(ss.years(D)) if fixed else ss.years(D)
    dis = ss.Disease(name='tx', Y=Y, Tf=Tf, Z=Z, init=dict(Y=1.0),
                     start=ss.Transition(Y, Tf, rate=ss.peryear(1.0)),
                     fail=ss.Transition(Tf, Z, dur=dur),
                     tb_death=ss.death(Tf, rate=ss.peryear(d_rate)))
    kw = dict(n_agents=n_agents, rand_seed=seed) if solver == 'abm' else dict(total_pop=1)
    sim = ss.Sim(dis, solver=solver, start=2000, stop=2012, dt=dt, verbose=0, **kw).run()
    r = sim.results.tx
    died, failed = r['tb_death'].values.sum(), r['fail'].values.sum()
    expect = 1 - np.exp(-d_rate*D) if fixed else d_rate/(d_rate + 1/D)
    return dict(dt=step, fixed=fixed, solver=solver, seed=seed, died=died, failed=failed, share=died/(died + failed), exact=expect)


if __name__ == '__main__':
    jobs = []
    for step in ['month', 'week']:
        for fixed in [True, False]:
            jobs += [dict(step=step, fixed=fixed, solver='abm', seed=s) for s in seeds]
            jobs.append(dict(step=step, fixed=fixed, solver='ode'))
    df = pd.DataFrame(sc.parallelize(run, iterkwargs=jobs, ncpus=tbode4.ncpus))
    df['ratio'] = df.share/df.exact
    df.to_csv(outdir / 'fail_death_runs.csv', index=False)
    g = df.groupby(['dt', 'fixed', 'solver']).agg(ratio=('ratio', 'mean'), se=('ratio', 'sem'), reps=('ratio', 'size'), courses=('failed', 'mean')).reset_index()
    g.to_csv(outdir / 'fail_death.csv', index=False)
    print('Share of failing courses ending in TB death ÷ its exact value (agents: mean ± SE over seeds):')
    print(g.to_string(index=False))
