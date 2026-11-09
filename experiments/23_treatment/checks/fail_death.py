""" Can people whose treatment will fail face TB death while on a fixed (or exponential) 6-month course, in both solvers?
People enter treatment during the run (from Y at a constant rate), and the share of entrants who die on treatment is compared with its exact value. """
import numpy as np, starsim as ss
d_rate, D = 0.3959, 0.5 # exp 21 set 1's sym_dead; course length (years)
for dt, label in [(ss.months(1), 'monthly'), (ss.weeks(1), 'weekly')]:
    for fixed in [True, False]:
        for solver in ['abm', 'ode']:
            Y, Tf, Z = ss.states('Y Tf Z')
            dur = ss.constant(ss.years(D)) if fixed else ss.years(D)
            dis = ss.Disease(name='tx', init=dict(Y=1.0),
                start    = ss.Transition(Y, Tf, rate=ss.peryear(1.0)),
                fail     = ss.Transition(Tf, Z, dur=dur),
                tb_death = ss.death(Tf, rate=ss.peryear(d_rate)))
            sim = ss.Sim(dis, solver=solver, n_agents=100000, start=2000, stop=2012, dt=dt, verbose=0, rand_seed=1).run()
            r = sim.results.tx
            died, failed = r['tb_death'].values.sum(), r['fail'].values.sum()
            expect = 1 - np.exp(-d_rate*D) if fixed else d_rate/(d_rate + 1/D)
            print(f'{label} fixed={fixed} {solver}: share of courses ending in TB death {died/(died+failed):.4f} (exact {expect:.4f}; ratio {died/(died+failed)/expect:.3f})')
