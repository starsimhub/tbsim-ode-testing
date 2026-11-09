"""
Treatment (tbode4.TBTreat, exp 23) against the standalone ODE (tbode4.fastode): v4's ODE stays at fastode's treated equilibrium, with
a course that is a rate or Erlang stages; first-entry incidence rebuilt from v4's ODE (ode_incidence) equals fastode's; and with
treatment off, TBTreat is tbsim.TB. Rerun after bumping the v4 pins.

    PYTHONPATH=v4 v4/.venv/bin/python -m pytest v4/tests/test_treat.py
"""
import numpy as np
import starsim as ss
import tbsim
import tbode4
from tbode4 import fastode as fo

# Exp 21's set 1 (the starter script's parameters), with tb_LAI_TPT's 2024 treatment
pi, phi, inf_non, inf_asy = 0.3518, 0.01194, 0.3153, 0.3735
pars = dict(beta=22.10, trans_asymp=0.7488, rr_reinfection_cleared=pi, rr_reinfection_rec=pi, rr_reinfection_treat=pi, rr_reinfection_inf=pi,
            inf_cle=0.6638, inf_non=inf_non, inf_asy=inf_asy, k_latent=7.517, inf_non_remote=phi*inf_non, inf_asy_remote=phi*inf_asy,
            non_rec=0.2207, non_asy=0.2776, asy_non=1.612, asy_sym=0.7560, sym_asy=0.7593, sym_dead=0.3959, mu=1/70)
tx = dict(tx_rate=float(tbode4.theta_from_coverage(0.677)), tx_success=0.53, tx_p_relapse=0.034)


def start_dict(x, p):
    """ fastode's state vector as TBTreat's init= (NON0 + NON1 as non_infectious; Erlang stages summed into each treatment state) """
    d = dict(zip(fo.names(p), x))
    out = {k: d[k] for k in fo.state_names if k not in ['non0', 'non1']}
    out['non_infectious'] = d['non0'] + d['non1']
    out['treated'], out['relapsing'] = d['treated'], d['relapsing']
    k = p['tx_stages']
    out['treatment_success'] = sum(d[f'treatment_success_{i}'] for i in range(k))
    out['treatment_fail'] = sum(d[f'treatment_fail_{i}'] for i in range(k))
    return out


def run_v4(init, years=50, **kwargs):
    tb = tbode4.TBTreat(init=init, **tx, **pars, **kwargs)
    return ss.Sim(tb, solver='ode', start=0, stop=years, dt=ss.months(1), total_pop=1, verbose=0).run()


def test_equilibrium_rate():
    """ A course at a rate: v4's ODE stays at fastode's treated equilibrium (one Erlang stage) """
    p = fo.make_pars(**pars, **tx, tx_stages=1)
    x = fo.equilibrium(p)
    sim = run_v4(start_dict(x, p), course='rate')
    r, alive = sim.results.tb, sim.results.alive.values[-1]
    for state, value in start_dict(x, p).items():
        assert abs(r[state].values[-1]/alive - value) < 1e-8, f'{state} drifts from the equilibrium'


def test_equilibrium_stages():
    """ A course of 4 Erlang stages in v4's ODE (course='fixed', stages=4) equals fastode's with tx_stages=4, from the start of a course """
    p = fo.make_pars(**pars, **tx, tx_stages=4)
    x = fo.equilibrium(p)
    x0 = x.copy()
    for first in [fo.TX, fo.TX + 4]: # v4 starts each treatment state's people at its first stage, so do the same in fastode
        x0[first] = x[first:first + 4].sum(); x0[first + 1:first + 4] = 0
    t = np.arange(0, 50 + 1e-9, 1/12)
    _, X = fo.run(p, x0, years=50, t_eval=t)
    sim = run_v4(start_dict(x, p), course='fixed', stages=4)
    r, alive = sim.results.tb, np.asarray(sim.results.alive.values, dtype=float)
    for state in ['symptomatic', 'asymptomatic', 'treated', 'relapsing']:
        diff = np.max(np.abs(np.asarray(r[state].values, dtype=float)/alive - X[:, fo.names(p).index(state)]))
        assert diff < 1e-7, f'{state}: max difference {diff:.2e}'


def test_ode_incidence():
    """ Incidence rebuilt from v4's ODE (first entries, plus relapses) equals fastode's at the treated equilibrium """
    p = fo.make_pars(**pars, **tx, tx_stages=1)
    x = fo.equilibrium(p)
    sim = run_v4(start_dict(x, p), course='rate', years=20)
    out = tbode4.tx_outputs(sim, tbode4.ode_incidence(sim, non0=x[fo.N0]))
    f = fo.summary(x, p)
    for key in ['incidence', 'notifications', 'mortality', 'relapse']:
        rel = abs(np.mean(out[key][1:])/f[key] - 1)
        assert rel < 1e-6, f'{key}: v4 {np.mean(out[key][1:]):.6g} against fastode {f[key]:.6g}'


def test_off_is_tbsim():
    """ With no treatment, TBTreat runs exactly as tbsim.TB """
    init = dict(infection=0.05, symptomatic=0.005)
    runs = []
    for model in [tbsim.TB(init=init, **pars), tbode4.TBTreat(init=init, **pars)]:
        runs.append(ss.Sim(model, solver='ode', start=0, stop=30, dt=ss.months(1), total_pop=1, verbose=0).run())
    for state in ['symptomatic', 'infection', 'non_infectious']:
        assert np.allclose(runs[0].results.tb[state].values, runs[1].results.tb[state].values, rtol=0, atol=1e-14)


if __name__ == '__main__':
    test_equilibrium_rate()
    test_equilibrium_stages()
    test_ode_incidence()
    test_off_is_tbsim()
