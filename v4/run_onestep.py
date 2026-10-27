"""
One step of the natural history as agents, against the exact one-step probabilities expm(Q·dt) (findings 4–5 of exps 04–05)

Everyone starts in one state, with no transmission, births, or background deaths, and the agents take one step. The share
ending in each state (or dead of TB) is compared with the row of expm(Q·dt) for that state, which is what the ODE (and
Starsim 3's TBExact) gives. Starsim 3's tbsim.TB took transitions in a fixed order, so people moving forward could carry on
within a step and people moving back could not; v4 draws each state's exits jointly and gives people who arrive in a state
during a step half a step of its rates.
"""
import numpy as np
import sciris as sc
import starsim as ss
from scipy.linalg import expm
import common as c

n_agents = 1_000_000
reps = 2
states = ['infection', 'cleared', 'recovered', 'non_infectious', 'asymptomatic', 'symptomatic', 'dead']
starts = ['infection', 'non_infectious', 'asymptomatic', 'symptomatic']
flows = [ # (from, to, parameter)
    ('infection', 'cleared', 'inf_cle'), ('infection', 'non_infectious', 'inf_non'), ('infection', 'asymptomatic', 'inf_asy'),
    ('non_infectious', 'recovered', 'non_rec'), ('non_infectious', 'asymptomatic', 'non_asy'),
    ('asymptomatic', 'non_infectious', 'asy_non'), ('asymptomatic', 'symptomatic', 'asy_sym'),
    ('symptomatic', 'asymptomatic', 'sym_asy'), ('symptomatic', 'dead', 'sym_dead'),
]
dts = sc.objdict(month=ss.months(1), week=ss.weeks(1))


def generator(parset):
    import tbsim
    value = lambda v: v.raw.value if isinstance(v, ss.par) else v.value if isinstance(v, ss.TimePar) else v
    pars = sc.mergedicts({k: value(v) for k, v in tbsim.tb_pars().items()}, c.parsets[parset])
    Q = np.zeros((len(states), len(states)))
    for src, dst, key in flows:
        i, j = states.index(src), states.index(dst)
        Q[i, j] += pars[key]
        Q[i, i] -= pars[key]
    return Q


def run_one(parset, dt, start, seed, solver='abm'):
    tb = c.make_tb(parset, beta=0, init={start: 1.0}, mu=0)
    sim = c.make_sim(tb, solver=solver, dt=dts[dt], start=0, stop=dts[dt].years*1, n_agents=n_agents, seed=seed)
    sim.run()
    r = sim.results.tb
    n0 = sim.results.alive.values[0]
    shares = {s: r[s].values[1]/n0 for s in states if s != 'dead'}
    shares['dead'] = 1 - sum(shares.values())
    return shares


if __name__ == '__main__':
    rows = []
    for parset in ['A', 'B']:
        Q = generator(parset)
        for dt in dts:
            P = expm(Q*dts[dt].years)
            for start in starts:
                runs = sc.parallelize(run_one, iterkwargs=[dict(parset=parset, dt=dt, start=start, seed=s) for s in range(reps)], ncpus=c.ncpus)
                ode = run_one(parset, dt, start, 0, solver='ode')
                i = states.index(start)
                for s in states:
                    vals = np.array([o[s] for o in runs])
                    exact = P[i, states.index(s)]
                    se = np.sqrt(max(exact*(1 - exact), 1e-12)/(n_agents*reps))
                    rows.append(dict(parset=parset, dt=dt, start=start, state=s, exact=exact, ode=ode[s], abm=vals.mean(), z=(vals.mean() - exact)/se))
                print(f'{parset} {dt:5s} from {start}: done', flush=True)
    df = sc.dataframe(rows)
    df.to_csv(c.folder / 'results' / 'onestep.csv', index=False)
    print(df.to_string())
