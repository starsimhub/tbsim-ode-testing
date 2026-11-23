"""
Reference runs of tbode4.TBTreat as it was before exp 35's options (commit 2550e63), for test_treat.py's regression test: the new defaults must
reproduce them step for step. Run once, before changing TBTreat:

    PYTHONPATH=v4 v4/.venv/bin/python v4/tests/make_tbtreat_reference.py
"""
import numpy as np
import sciris as sc
import starsim as ss
import tbode4

here = sc.thispath(__file__)
pi, phi, inf_non, inf_asy = 0.3518, 0.01194, 0.3153, 0.3735
pars = dict(beta=22.10, trans_asymp=0.7488, rr_reinfection_cleared=pi, rr_reinfection_rec=pi, rr_reinfection_treat=pi, rr_reinfection_inf=pi,
            inf_cle=0.6638, inf_non=inf_non, inf_asy=inf_asy, k_latent=7.517, inf_non_remote=phi*inf_non, inf_asy_remote=phi*inf_asy,
            non_rec=0.2207, non_asy=0.2776, asy_non=1.612, asy_sym=0.7560, sym_asy=0.7593, mu=1/70)
init = dict(infection=0.012, infection_remote=0.25, cleared=0.2, recovered=0.03, non_infectious=0.003, asymptomatic=0.002, symptomatic=0.001)
cases = sc.objdict(
    ode_month_fixed = dict(solver='ode', dt=ss.months(1), course='fixed'),
    ode_month_rate  = dict(solver='ode', dt=ss.months(1), course='rate'),
    ode_week_fixed  = dict(solver='ode', dt=ss.weeks(1), course='fixed'),
    abm_month_fixed = dict(solver='abm', dt=ss.months(1), course='fixed'),
)


def make(solver, dt, course, start=1990, stop=2010):
    """ TBTreat with tb_LAI_TPT's treatment scale-up from 1995 and a falling sym_dead """
    tx = tbode4.lai_tpt_treatment(step=dt.years, end=stop)
    sym_dead = ss.peryear(ss.timeseries({1995.0: 0.3959, 2024.0: 0.1426}, interpolation='linear'))
    tb = tbode4.TBTreat(init=init, course=course, tx_rate=tx.tx_rate, tx_success=tx.tx_success, tx_p_relapse=tx.p_relapse, sym_dead=sym_dead, **pars)
    kw = dict(n_agents=20_000, rand_seed=35) if solver == 'abm' else dict(total_pop=1)
    return ss.Sim(tb, solver=solver, start=start, stop=stop, dt=dt, verbose=0, **kw)


def results(sim):
    """ Every TB result (states and flows) as arrays, and alive """
    r = sim.results.tb
    out = {f'tb.{k}': np.asarray(r[k].values, dtype=float) for k in r.keys() if hasattr(r[k], 'values')}
    out['alive'] = np.asarray(sim.results.alive.values, dtype=float)
    return out


if __name__ == '__main__':
    ref = {}
    for name, c in cases.items():
        sim = make(**c).run()
        ref.update({f'{name}/{k}': v for k, v in results(sim).items()})
        print(name, len(results(sim)), 'results', flush=True)
    np.savez_compressed(here / 'data' / 'tbtreat_reference.npz', **ref)
