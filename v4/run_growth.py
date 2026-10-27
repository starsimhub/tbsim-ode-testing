"""
Growth from a seed, agents vs the ODE (the open question from finding 7 of exp 07: does the step's infection timing make
each generation lag?)

Each arm starts with 0.1% of people in symptomatic TB and everyone else susceptible, at the calibrated `beta`, and runs for
60 years. Starsim 3's TBExact had to choose between dating infections at the start of the step (right generation interval,
end-of-step latent snapshots low) and within it (snapshots right, generation interval half a step too long). Under v4, an
infection is dated at the end of its step, and people who become infectious have half a step of their state's rates at their
first draw, so that each stay counts in transmission for as long as it lasts (to first order in dt). A lag would show as
agents reaching each prevalence later than the ODE, by more at a monthly step than a weekly one.
"""
import numpy as np
import sciris as sc
import starsim as ss
import common as c
from calibrate import arms
from run_equilibrium import dts

n_agents = 100_000
reps = 8
start, stop = 1900, 1960
seed_share = 0.001
design = [('Bk', 'month'), ('Bk', 'week')]


def run_one(arm, dt, seed, solver='abm'):
    cal = c.get_calib(arm, **arms[arm])
    spec = arms[arm]
    tb = c.make_tb(spec['parset'], cal.beta, init=dict(symptomatic=seed_share), model=spec['model'])
    sim = c.make_sim(tb, solver=solver, dt=dts[dt], start=start, stop=stop, n_agents=n_agents, seed=seed)
    sim.run()
    out = c.outputs(sim)
    return sc.objdict(year=out.year, prevalence=out.prevalence, infection=out.infection)


def crossing(out, level):
    """ The first time prevalence reaches `level`, interpolated linearly """
    i = np.argmax(out.prevalence >= level)
    if i == 0:
        return np.nan
    y0, y1 = out.prevalence[i-1], out.prevalence[i]
    return out.year[i-1] + (level - y0)/(y1 - y0)*(out.year[i] - out.year[i-1])


if __name__ == '__main__':
    fn = c.folder / 'results' / 'growth.obj'
    data = sc.objdict(ode=sc.objdict(), abm=sc.objdict(), versions=c.versions())
    for arm in sorted({d[0] for d in design}):
        data.ode[arm] = run_one(arm, 'week', 0, solver='ode')
    for arm, dt in design:
        data.abm[f'{arm}_{dt}'] = sc.parallelize(run_one, iterkwargs=[dict(arm=arm, dt=dt, seed=s) for s in range(reps)], ncpus=c.ncpus)
        sc.save(fn, data)
        for level in [0.0025, 0.004]:
            ode = crossing(data.ode[arm], level)
            t = np.array([crossing(o, level) for o in data.abm[f'{arm}_{dt}']])
            print(f'{arm} {dt:5s} reaches {1e5*level:.0f}/100k: ODE {ode:.2f}, agents {np.nanmean(t):.2f} ± {np.nanstd(t, ddof=1)/np.sqrt(np.isfinite(t).sum()):.2f} (lag {np.nanmean(t) - ode:+.2f} y)', flush=True)
