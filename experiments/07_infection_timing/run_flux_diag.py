"""
Experiment 07 diagnostic: per-step flows in weekly TBExact runs, 'within' vs 'start' (12 seeds each, 15 years, 50k agents).
Infections from SUSCEPTIBLE / CLEARED vs the ODE expectation, new infections' end-of-step states, and natural-history flows. Output: run_flux_diag.log.
"""
import numpy as np, sciris as sc, tbode
from tbsim import TBS
beta, dt = 9.917337993130424, 7/365

def run(timing, seed):
    ode = tbode.TBODE(beta=beta, replace_tb_deaths=False, init_sym=0.0)
    flux = tbode.FluxCounts(beta_ode=beta, name='flux')
    sim = tbode.make_abm(ode.pars, seed=seed, n_agents=50_000, start=1900, stop=1915, dt=dt, analyzers=flux, tb_class=tbode.TBExact, infection_timing=timing)
    sim.init(); tbode.init_from_ode(sim, ode.equilibrium(), seed=seed)
    f = sim.analyzers.flux; f.snapshot(); sim.run()
    return dict(flows=f.flows.sum(0), stock=f.stock.sum(0), exp_S=f.exp_inf_S.sum(), exp_C=f.exp_inf_C.sum())

if __name__ == '__main__':
    seeds = range(12)
    outs = sc.parallelize(run, iterkwargs=[dict(timing=t, seed=s) for t in ['within', 'start'] for s in seeds], ncpus=4)
    for i, t in enumerate(['within', 'start']):
        oo = outs[i*len(seeds):(i+1)*len(seeds)]
        F = sum(o['flows'] for o in oo); st = sum(o['stock'] for o in oo)
        for src, ek in [(TBS.SUSCEPTIBLE, 'exp_S'), (TBS.CLEARED, 'exp_C')]:
            row = F[src]; inf = row.sum() - row[src] - row[TBS.DEAD]
            exp = sum(o[ek] for o in oo)
            ends = {s.name: int(row[s]) for s in [TBS.INFECTION, TBS.NON_INFECTIOUS, TBS.ASYMPTOMATIC, TBS.SYMPTOMATIC, TBS.CLEARED] if s != src}
            print(f'{t:6s} from {src.name:11s}: infections {inf} / expected {exp:.0f} = {inf/exp:.4f} ± {1/np.sqrt(inf):.4f}; end states {ends}; stays {int(row[src])}')
        for src in [TBS.INFECTION, TBS.NON_INFECTIOUS, TBS.ASYMPTOMATIC, TBS.SYMPTOMATIC]:
            print(f'   {t} {src.name}: stock {st[src]}  ->', {s.name: round(F[src, s]/st[src], 6) for s in [TBS.INFECTION, TBS.CLEARED, TBS.NON_INFECTIOUS, TBS.ASYMPTOMATIC, TBS.SYMPTOMATIC, TBS.DEAD] if F[src, s]})
