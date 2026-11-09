import numpy as np, starsim as ss
class Track(ss.Analyzer):
    """ Per step: agents who failed and restarted in the same step """
    def init_post(self):
        super().init_post(); self.both = []; self.fails = []
    def step(self):
        d = self.sim.diseases.tx
        f, s = np.asarray(d.fail.uids), np.asarray(d.start.uids)
        self.both.append(len(np.intersect1d(f, s))); self.fails.append(len(f))

def make(fixed=True, s=0.5, theta=12.0, pulse=False):
    Y, Ts, Tf, R, L = ss.states('Y Ts Tf R L')
    D = ss.constant(ss.years(0.5)) if fixed else ss.years(0.5)
    rate = ss.peryear(theta)
    d = ss.Disease(name='tx', init=dict(Y=1.0),
        start   = ss.Transition(Y, [Ts, Tf], rate=rate, p=[s, None]),
        success = ss.Transition(Ts, [R, L], dur=D, p=[None, 0.1]),
        fail    = ss.Transition(Tf, Y, dur=D),
        relapse = ss.Transition(L, Y, dur=ss.years(1.5)))
    return d

def run(solver, fixed, s=0.5, theta=12.0, mods=None, analyzers=None, stop=2003):
    d = make(fixed, s, theta)
    sim = ss.Sim(d, solver=solver, n_agents=20000, start=2000, stop=stop, dt=ss.months(1), verbose=0, rand_seed=3, interventions=mods(d) if mods else None, analyzers=analyzers)
    sim.run(); return sim

# 1. Fixed course timing: a pulse of starts in the first month only
pulse = lambda d: ss.modify(d.start, multiply=ss.timeseries({'2000-01-01': 1.0, '2000-02-01': 0.0}))
for solver in ['abm', 'ode']:
    for fixed in [True, False]:
        sim = run(solver, fixed, mods=pulse)
        r = sim.results.tx; ex = r['fail'].values + r['success'].values
        print(f'[timing] {solver} fixed={fixed}: starts by month {np.round(r["start"].values[:3]).astype(int)}; exits months 0-9 {np.round(ex[:10]).astype(int)}')

# 2. Time-varying efficacy: s 0.2 in 2000, 0.8 from 2001 -- is it applied at start?
sts = ss.timeseries({'2000-01-01': 0.2, '2001-01-01': 0.8})
for solver in ['abm', 'ode']:
    try:
        sim = run(solver, True, s=sts, theta=2.0)
        r = sim.results.tx
        new_ts = np.diff(r['Ts'].values, prepend=0) + r['success'].values  # inflow to Ts each step
        st = r['start'].values; yr = np.array(sim.t.yearvec)
        for lo, hi in [(2000.1, 2001), (2001, 2002)]:
            m = (yr >= lo) & (yr < hi)
            print(f'[s(t) as p] {solver} starts {lo}-{hi}: share to Ts = {new_ts[m].sum()/st[m].sum():.3f}')
        # success exits in first half of 2001 come from 2000 H2 starts
        m = (yr >= 2001.05) & (yr < 2001.45); sx, fx = r['success'].values[m].sum(), r['fail'].values[m].sum()
        print(f'[s(t) as p] {solver} exits Jan-May 2001: success share {sx/(sx+fx):.3f} (0.2 if drawn at start, 0.8 if at exit)')
    except Exception as e:
        print(f'[s(t) as p] {solver}: ERROR {type(e).__name__}: {str(e)[:400]}')

# 3. Same-step failure and restart, agents
for fixed in [True, False]:
    tr = Track()
    sim = run('abm', fixed, theta=12.0, analyzers=tr)
    a = sim.analyzers[0]
    print(f'[same step] fixed={fixed}: failures {sum(a.fails)}, failed and restarted in the same step {sum(a.both)} ({sum(a.both)/max(1,sum(a.fails)):.3f}; theta*dt = {12/12:.2f}, 1-exp = {1-np.exp(-1):.3f})')
