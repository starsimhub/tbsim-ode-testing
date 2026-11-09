import numpy as np, starsim as ss, tbsim
for solver in ['abm', 'ode']:
    for dur in [ss.years(1.5), ss.constant(ss.years(1.5))]:
        tb = tbsim.TB(init=dict(relapsing=1.0), beta=0, p_relapse=0.1, delta=2.0, sym_dead=0, sym_asy=0, mu=0, dur_relapse=dur)
        sim = ss.Sim(tb, solver=solver, n_agents=20000, start=2000, stop=2004, dt=ss.months(1), verbose=0, rand_seed=1)
        sim.run()
        r = sim.results.tb
        left = r['relapsing'].values/r['relapsing'].values[0]
        yrs = np.array(sim.t.yearvec) - 2000
        i = [np.argmin(abs(yrs - y)) for y in [0.75, 1.5, 3.0]]
        print(f'{solver} {type(dur).__name__:9s} share still relapsing at 0.75/1.5/3 y:', np.round(left[i], 3), ' exp:', np.round(np.exp(-np.array([0.75,1.5,3])/1.5), 3))
