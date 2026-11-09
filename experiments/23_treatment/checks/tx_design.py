""" Does the proposed design run in both solvers? Y -> {Ts, Tf} (outcome drawn at entry), Ts fixed 6m -> {R, L} (relapse branch), Tf fixed 6m -> Y, L -> Y at rate """
import numpy as np, starsim as ss
def make(fixed=True, s=0.5, theta=12.0, p_rel=0.1, svary=None):
    Y, Ts, Tf, R, L = ss.states('Y Ts Tf R L')
    D = ss.constant(ss.years(0.5)) if fixed else ss.years(0.5)
    tr = dict(
        start   = ss.Transition(Y, [Ts, Tf], rate=ss.peryear(theta), p=[s, None]),
        success = ss.Transition(Ts, [R, L], dur=D, p=[None, p_rel]),
        fail    = ss.Transition(Tf, Y, dur=D),
        relapse = ss.Transition(L, Y, dur=ss.years(1.5)),
    )
    d = ss.Disease(name='tx', init=dict(Y=1.0), **tr)
    return d
for solver in ['abm', 'ode']:
    for fixed in [True, False]:
        try:
            d = make(fixed)
            sim = ss.Sim(d, solver=solver, n_agents=20000, start=2000, stop=2003, dt=ss.months(1), verbose=0, rand_seed=2)
            sim.run()
            r = sim.results.tx
            keys = [k for k in r.keys()]
            tot = lambda k: r[k].values.sum()
            print(f'{solver} fixed={fixed}: starts {tot("start"):.0f}, success-exits {tot("success"):.0f}, fail-exits {tot("fail"):.0f}; R {r["R"].values[-1]:.0f}, L+relapsed share of success-exits ~ {(tot("success") and (r["L"].values[-1]+tot("relapse"))/tot("success")):.3f}')
        except Exception as e:
            print(f'{solver} fixed={fixed}: ERROR {type(e).__name__}: {str(e)[:400]}')
