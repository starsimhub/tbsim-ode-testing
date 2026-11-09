import numpy as np, starsim as ss
def run(solver, mod):
    Y, T = ss.states('Y T')
    d = ss.Disease(name='tx', init=dict(Y=1.0), start=ss.Transition(Y, T, rate=ss.peryear(1.2)))
    mods = ss.modify(d.start, multiply=ss.timeseries({'2000-01-01': 1.0, '2000-07-01': 0.0})) if mod else None
    sim = ss.Sim(d, solver=solver, n_agents=100000, start=2000, stop=2001, dt=ss.months(1), verbose=0, rand_seed=1, interventions=mods)
    sim.run()
    r = sim.results.tx
    print(f'{solver} mod={mod}: dates', [str(x)[:7] for x in sim.t.datevec[:9]])
    print('   starts', np.round(r['start'].values[:9]).astype(int), ' Y', np.round(r['Y'].values[:9]).astype(int))
for solver in ['abm', 'ode']:
    for mod in [False, True]:
        run(solver, mod)
