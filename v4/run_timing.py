"""
Run time of the agents and the ODE (finding 6 of exp 06), for arm A0 at 50,000 agents

One process, one run at a time, best of 3 after a run that compiles. Exp 06 measured Starsim 3's tbsim.TB and TBExact at
about 10.5 ms a step at 50,000 agents; that was another code base and another day's load on the laptop, so it is a
reference, not an A/B comparison.
"""
import sciris as sc
import starsim as ss
import common as c
from calibrate import arms
from run_equilibrium import dts

n_agents = 50_000
years = sc.objdict(day=5, week=20, month=50)
repeats = 3


def time_one(arm, dt, solver):
    cal = c.get_calib(arm, **arms[arm])
    spec = arms[arm]
    tb = c.make_tb(spec['parset'], cal.beta, init=cal.equilibrium, model=spec['model'])
    sim = c.make_sim(tb, solver=solver, dt=dts[dt], start=1900, stop=1900 + (100 if solver == 'ode' else years[dt]), n_agents=n_agents)
    T = sc.timer()
    sim.run()
    return T.toc(output=True), len(sim.t.yearvec) - 1


if __name__ == '__main__':
    rows = []
    for arm in ['A0', 'Bk']:
        for dt in dts:
            for solver in ['abm', 'ode']:
                time_one(arm, dt, solver) # Compiles
                best, n = min(time_one(arm, dt, solver) for _ in range(repeats))
                rows.append(dict(arm=arm, dt=dt, solver=solver, steps=n, seconds=best, ms_per_step=1e3*best/n))
                print(rows[-1], flush=True)
    df = sc.dataframe(rows)
    df.to_csv(c.folder / 'results' / 'timing.csv', index=False)
    print(df.to_string())
