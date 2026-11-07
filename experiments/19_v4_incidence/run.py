"""
Experiment 19, v4 environment: incidence as the first entry to asymptomatic TB per episode, in v4 agents and in the standalone ODE
(tbode4.fastode), against v4's every-entry count. See README.md. Run from the repo root:

    v4/.venv/bin/python experiments/19_v4_incidence/run.py           # all arms (resumes)
    v4/.venv/bin/python experiments/19_v4_incidence/run.py summary   # summary only

The standalone ODE's agreement with v4's ODE (question 3) is v4/tests/test_fastode.py.
"""
import sys
import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss
from scipy.optimize import brentq

outdir = sc.thispath(__file__)
root = outdir.parent.parent
sys.path.insert(0, str(root / 'v4'))
sys.path.insert(0, str(root / 'experiments' / '17_latent_reinfection'))
import common as c # noqa: E402
import run_equilibrium as req # noqa: E402
import tbode4 # noqa: E402
from tbode4 import fastode as fo # noqa: E402
import run_v4 as rv17 # noqa: E402

dts = sc.objdict(month=ss.months(1), week=ss.weeks(1))
reps, seed0, n_agents = 40, 19000, 50_000
start, stop = 1900, 2000
fn = outdir / 'results.obj'


def set_pars(name):
    """ (tbsim parameter set name, extra tbsim.TB parameters) for set A, or set B two-state at phi = 0.01 (exp 17's model) """
    if name == 'A':
        return 'A', {}
    return 'B', {k: v for k, v in rv17.two_state_pars(0.01).items() if k != 'fit'}


def fast_pars(name, beta):
    parset, extra = set_pars(name)
    return fo.make_pars(beta=beta, **sc.mergedicts(c.parsets[parset], extra))


def calibrate(name, target=0.005):
    """ beta for 500/100k active prevalence at the standalone ODE's equilibrium, and that equilibrium """
    prev = lambda b: (lambda x: 0 if x is None else x[fo.A] + x[fo.Y])(fo.equilibrium(fast_pars(name, b)))
    beta = brentq(lambda b: prev(b) - target, 1, 60, xtol=1e-8)
    p = fast_pars(name, beta)
    return beta, p, fo.equilibrium(p)


class Incidence(ss.Analyzer):
    """
    Incidence from transitions, two ways: every entry to asymptomatic TB from latent infection or non-infectious TB (v4's
    count, here including remote latent -> asymptomatic), and the first entry per episode. A flag per agent records whether
    they have been asymptomatic in the current episode: set on any entry to asymptomatic TB, cleared on recovery. Agents who
    start in non-infectious TB start flagged with the equilibrium's share NON1 / (NON0 + NON1).
    """
    def __init__(self, share_flagged=0.0, **kwargs):
        super().__init__(**kwargs)
        self.share_flagged = share_flagged
        return

    def init_post(self):
        super().init_post()
        npts = self.sim.t.npts
        self.flag = np.zeros(2*len(self.sim.people), bool) # Grown as births add agents (see _ensure)
        tb = self.sim.diseases.tb
        rng = np.random.default_rng(int(self.sim.pars.rand_seed) + 7)
        for st in ['asymptomatic', 'symptomatic']:
            self._set(np.asarray(getattr(tb, st).uids), True)
        non = np.asarray(tb.non_infectious.uids)
        self._set(non[rng.random(len(non)) < self.share_flagged], True)
        self.first, self.every, self.alive = np.zeros(npts), np.zeros(npts), np.zeros(npts)
        return

    def _ensure(self, u):
        """ Grow the flag array (False for new agents) to hold every uid in u """
        if len(u) and u.max() >= len(self.flag):
            self.flag = np.concatenate([self.flag, np.zeros(max(len(self.flag), u.max() + 1 - len(self.flag)), bool)])
        return

    def _set(self, u, val):
        self._ensure(u)
        self.flag[u] = val
        return

    def step(self):
        tb = self.sim.diseases.tb
        get = lambda name: np.asarray(getattr(tb, name).uids) if hasattr(tb, name) else np.array([], int)
        latent_asy = np.concatenate([get('latent_asy'), get('remote_asy')])
        worsening, remission, recovery = get('worsening'), get('remission'), get('recovery')
        self._ensure(np.concatenate([latent_asy, worsening, remission, recovery]))
        w_first = worsening[~self.flag[worsening] & ~np.isin(worsening, latent_asy)]
        ti = self.sim.ti
        self.first[ti] = len(latent_asy) + len(w_first)
        self.every[ti] = len(latent_asy) + len(worsening)
        self.alive[ti] = len(self.sim.people)
        self._set(np.concatenate([latent_asy, worsening, remission]), True)
        self._set(recovery, False)
        return


def run_one(name, dkey, seed, beta, x):
    """ One agent run from the standalone ODE's equilibrium x, reduced to yearly means of incidence per person-year """
    parset, extra = set_pars(name)
    init = dict(zip(fo.state_names, x))
    non0, non1 = init.pop('non0'), init.pop('non1')
    init['non_infectious'] = non0 + non1
    if not extra.get('k_latent'):
        init.pop('infection_remote')
    tb = c.make_tb(parset, beta, init=init, model='tb', **extra)
    sim = ss.Sim(tb, solver='abm', n_agents=n_agents, total_pop=n_agents, start=start, stop=stop, dt=dts[dkey], rand_seed=seed, verbose=0,
                 analyzers=Incidence(share_flagged=non1/(non0 + non1)))
    sim.run()
    an = sim.analyzers[0]
    r = sim.results.tb
    alive = sim.results.alive.values
    dt = sim.t.dt_year
    out = sc.objdict(year=np.array(sim.t.yearvec), first=an.first/alive/dt, every=an.every/alive/dt,
                     every_v4=sum(r[k].values for k in ['latent_asy', 'worsening'])/alive/dt, # v4/common.py's outputs(): no remote -> asymptomatic
                     prevalence=(r.asymptomatic.values + r.symptomatic.values)/alive, tb_deaths=r.tb_death.values/alive/dt)
    return req.yearly(out)


def main():
    data = sc.load(fn) if fn.exists() else sc.objdict(cal=sc.objdict(), runs=sc.objdict())
    for name in ['A', 'B2']:
        if name not in data.cal:
            beta, p, x = calibrate(name)
            data.cal[name] = sc.objdict(beta=beta, pars=p, x=x, summary=fo.summary(x, p), flows=fo.flows(x, p))
            sc.save(fn, data)
            s = data.cal[name].summary
            print(f'{name}: beta {beta:.5f}; ODE incidence per 100k: first entry {1e5*s.incidence:.2f}, every entry {1e5*s.into_asy:.2f} '
                  f'(re-entries {1 - s.incidence/s.into_asy:.3f}); prevalence ÷ incidence {s.prev_to_inc:.3f} y (every entry: {s.prevalence/s.into_asy:.3f})', flush=True)
        cal = data.cal[name]
        for dkey in dts:
            key = f'{name}_{dkey}'
            if key in data.runs and len(data.runs[key]) >= reps:
                continue
            T = sc.timer()
            data.runs[key] = sc.parallelize(run_one, iterkwargs=[dict(name=name, dkey=dkey, seed=seed0 + s, beta=cal.beta, x=cal.x) for s in range(reps)], ncpus=tbode4.ncpus)
            sc.save(fn, data)
            print(f'{key}: {reps} reps in {T.toc(output=True):.0f} s', flush=True)
    tbode4.save_versions(outdir, 'v4')
    summarise(data)
    return data


def summarise(data):
    """ Agents ÷ ODE over 1950–2000 (mean ± SE, 10–90% across reps), for each count """
    rows = []
    for key, runs in data.runs.items():
        name, dkey = key.split('_')
        s, f = data.cal[name].summary, data.cal[name].flows
        ref = dict(first=s.incidence, every=s.into_asy, every_v4=f.latent_asy + f.worsening, prevalence=s.prevalence, tb_deaths=s.mortality)
        for out, val in ref.items():
            r = np.array([o[out][(o.year >= 1950) & (o.year < 2000)].mean() for o in runs])/val
            rows.append(dict(set=name, dt=dkey, output=out, ode_per_100k=1e5*val, reps=len(r), mean=r.mean(), se=r.std(ddof=1)/np.sqrt(len(r)),
                             p10=np.percentile(r, 10), p90=np.percentile(r, 90)))
        # The re-entry share, as agents count it
        fr = np.array([o.first[(o.year >= 1950) & (o.year < 2000)].mean()/o.every[(o.year >= 1950) & (o.year < 2000)].mean() for o in runs])
        rows.append(dict(set=name, dt=dkey, output='first_over_every', ode_per_100k=s.incidence/s.into_asy, reps=len(fr), mean=fr.mean(), se=fr.std(ddof=1)/np.sqrt(len(fr)),
                         p10=np.percentile(fr, 10), p90=np.percentile(fr, 90)))
    df = pd.DataFrame(rows)
    df.to_csv(outdir / 'summary.csv', index=False)
    with pd.option_context('display.width', 200, 'display.max_rows', 200):
        print(df.round(4).to_string(index=False))
    return df


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'summary':
        summarise(sc.load(fn))
    else:
        main()
