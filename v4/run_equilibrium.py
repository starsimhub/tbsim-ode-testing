"""
Agents vs the ODE from the ODE's equilibrium (findings 1–3, 5, 8, 10 and 11 of exps 01–11, redone under Starsim v4)

Each arm's ODE `beta` is calibrated to 500/100k active prevalence (calibrate.py). The agents start at the ODE's equilibrium
shares and run for 100 years (1900–2000), with homogeneous mixing (v4's default with no network) or, where marked, a random
network with new contacts every step. Each run is reduced to yearly means of the outputs in `common.outputs()`.

Run with::

    python run_equilibrium.py baseline      # A0 and B0 at weekly and monthly dt, and A0 on a network (two per-contact rates)
    python run_equilibrium.py progression   # Ak (monthly), Bk and Bremote (weekly and monthly)
"""
import sys
import numpy as np
import sciris as sc
import starsim as ss
import common as c

dts = sc.objdict(day=ss.days(1), week=ss.weeks(1), month=ss.months(1))
start, stop = 1900, 2000
n_agents = 50_000
n_contacts = 10 # Mean number of contacts on the random network

designs = sc.objdict( # Abbreviated (Cliff, for time): no daily step, fewer reps than exps 01–11
    baseline = [ # (arm, dt, network, reps)
        ('A0', 'month', None, 40), ('A0', 'week', None, 20), ('A0', 'month', 'net', 20), ('A0', 'month', 'net_matched', 20),
        ('B0', 'month', None, 20), ('B0', 'week', None, 12),
    ],
    progression = [
        ('Ak', 'month', None, 16),
        ('Bk', 'month', None, 20), ('Bk', 'week', None, 12),
        ('Bremote', 'month', None, 20), ('Bremote', 'week', None, 12),
    ],
)


def yearly(out):
    """ Yearly means of each output (the first point, the starting state, is left out) """
    years = np.floor(out.year[1:] + 1e-9).astype(int)
    keep = sc.objdict(year=np.unique(years))
    for key, vals in out.items():
        if key != 'year':
            keep[key] = np.array([vals[1:][years == y].mean() for y in keep.year])
    return keep


def run_one(arm, dt, network, seed, solver='abm'):
    """
    One run. `network`: None for homogeneous mixing; 'net' for a random network of K contacts with the per-contact rate
    beta/K; 'net_matched' for the same network with the per-contact rate whose per-step probability, 1 − exp(−rate·dt), is
    the ODE's beta·dt/K (as `tbode.beta_edge_from_ode`, exp 10)
    """
    from calibrate import arms
    cal = c.get_calib(arm, **arms[arm])
    spec = arms[arm]
    K = n_contacts if network else None
    if network == 'net':
        beta = cal.beta/K
    elif network == 'net_matched':
        beta = -np.log(1 - cal.beta*dts[dt].years/K)/dts[dt].years
    else:
        beta = cal.beta
    tb = c.make_tb(spec['parset'], beta, init=cal.equilibrium, model=spec['model'], **{k:v for k,v in spec.items() if k not in ['parset', 'model']})
    sim = c.make_sim(tb, solver=solver, dt=dts[dt], start=start, stop=stop, n_agents=n_agents, seed=seed, network=K)
    T = sc.timer()
    sim.run()
    out = yearly(c.outputs(sim))
    out.time = T.toc(output=True)
    out.n_steps = len(sim.t.yearvec) - 1
    return out


if __name__ == '__main__':
    which = sys.argv[1]
    design = designs[which]
    fn = c.folder / 'results' / f'{which}.obj'
    data = sc.load(fn) if fn.exists() else sc.objdict(ode=sc.objdict(), abm=sc.objdict(), versions=c.versions())
    for arm in sorted({d[0] for d in design}): # The ODE, at a monthly step, as the reference
        if arm not in data.ode:
            data.ode[arm] = run_one(arm, 'month', None, 0, solver='ode')
    for arm, dt, network, reps in design:
        key = f'{arm}_{dt}' + (f'_{network}' if network else '')
        if key in data.abm and len(data.abm[key]) >= reps:
            continue
        T = sc.timer()
        data.abm[key] = sc.parallelize(run_one, iterkwargs=[dict(arm=arm, dt=dt, network=network, seed=s) for s in range(reps)], ncpus=c.ncpus)
        sc.save(fn, data)
        prev = np.array([o.prevalence[-50:].mean() for o in data.abm[key]])/data.ode[arm].prevalence[-50:].mean()
        print(f'{key:16s} {reps} reps in {T.toc(output=True):.0f} s; prevalence 1950–2000, agents ÷ ODE: {prev.mean():.3f} ± {prev.std(ddof=1)/np.sqrt(reps):.3f}', flush=True)
