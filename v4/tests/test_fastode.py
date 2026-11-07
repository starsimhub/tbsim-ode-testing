"""
The standalone ODE (tbode4.fastode) reproduces v4's TB ODE (tbsim.TB run by Starsim's compartmental solver): at equilibrium, and
over time from a seed. Rerun after bumping the v4 pins: a change to tbsim's declared model would show here.

    v4/.venv/bin/python -m pytest v4/tests/test_fastode.py
"""
import numpy as np
import sciris as sc
import starsim as ss
import tbsim
import tbode4
from tbode4 import fastode as fo

# Set A (tbsim's defaults) and set B two-state (tb_LAI_TPT's best fit with remote progression, exp 17's phi = 0.01 fit)
set_b = dict(inf_cle=0.5365965386150265, inf_non=0.35069644069791794, inf_asy=0.23431100544176175, inf_non_remote=0.0035069644069791795,
             inf_asy_remote=0.002343110054417618, k_latent=4.631151620621326, non_rec=0.1706711906227213, non_asy=0.2426759489829796,
             asy_non=1.738594655537684, asy_sym=0.7435451323482238, sym_asy=0.6630323475395543, trans_asymp=0.9407760328175468,
             sym_dead=0.2626559974664169, rr_reinfection_cleared=0.2415635040650483, rr_reinfection_rec=0.2415635040650483,
             rr_reinfection_treat=0.2415635040650483, rr_reinfection_inf=0.2415635040650483)
cases = dict(A=(dict(), 9.926475), B=(set_b, 21.14343))


def v4_sim(pars, beta, init, dt, years):
    tb = tbsim.TB(beta=ss.peryear(beta), init=init, **{k: (ss.peryear(v) if not k.startswith(('rr_', 'trans')) else v) for k, v in pars.items()})
    sim = ss.Sim(tb, solver='ode', dt=dt, start=0, stop=years, total_pop=1, verbose=0)
    sim.run()
    return sim


def shares(sim, i=-1):
    r = sim.results.tb
    alive = sim.results.alive.values[i]
    return {s.name: float(r[s.name].values[i]/alive) for s in sim.diseases.tb.graph.people_states}


def test_equilibrium():
    """ The root-found equilibrium is a fixed point of v4's ODE: run from it for 50 years, v4's shares stay within 1e-9 """
    for name, (pars, beta) in cases.items():
        p = fo.make_pars(beta=beta, **pars)
        x = fo.equilibrium(p)
        start = dict(zip(fo.state_names, x))
        start['non_infectious'] = start.pop('non0') + start.pop('non1')
        if not pars.get('k_latent'):
            start.pop('infection_remote')
        sim = v4_sim(pars, beta, start, ss.years(1), 50)
        end = shares(sim)
        diff = max(abs(end[k] - start.get(k, 0.0)) for k in end)
        assert diff < 1e-9, f'Set {name}: v4 drifts {diff:.2e} from the standalone equilibrium'
        assert abs(x[fo.A] + x[fo.Y] - 0.005) < 1e-4, f'Set {name}: prevalence {x[fo.A] + x[fo.Y]:.5f}, not about 500/100k'


def test_trajectory():
    """ From a 1% symptomatic seed, 100 years at weekly steps: prevalence and every state within 1e-8 of v4's """
    for name, (pars, beta) in cases.items():
        p = fo.make_pars(beta=beta, **pars)
        sim = v4_sim(pars, beta, dict(symptomatic=0.01), ss.weeks(1), 100)
        t = np.asarray(sim.t.yearvec, dtype=float)
        _, X = fo.run(p, dict(symptomatic=0.01), years=100, t_eval=t)
        r = sim.results.tb
        alive = np.asarray(sim.results.alive.values, dtype=float)
        for s in ['susceptible', 'infection', 'cleared', 'recovered', 'asymptomatic', 'symptomatic']:
            diff = np.max(np.abs(X[:, fo.state_names.index(s)] - np.asarray(r[s].values, dtype=float)/alive))
            assert diff < 1e-8, f'Set {name}, {s}: max difference {diff:.2e}'
        non = X[:, fo.N0] + X[:, fo.N1]
        assert np.max(np.abs(non - np.asarray(r.non_infectious.values, dtype=float)/alive)) < 1e-8


def test_incidence_split():
    """ NON0 + NON1 changes nothing: the split's total follows v4's non_infectious, and first entries never exceed every entry """
    for name, (pars, beta) in cases.items():
        p = fo.make_pars(beta=beta, **pars)
        x = fo.equilibrium(p)
        f = fo.flows(x, p)
        assert 0 < f.incidence <= f.into_asy
        assert np.isclose(f.into_asy - f.incidence, p['non_asy']*x[fo.N1])


if __name__ == '__main__':
    test_equilibrium()
    test_trajectory()
    test_incidence_split()
    print('Passed')
