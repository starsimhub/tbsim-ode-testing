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


#%% Exp 35: death and loss to follow-up on treatment

import sys
import sciris as sc
from tbode4 import fastode_age as fa
sys.path.insert(0, str(sc.thispath(__file__)))
import make_tbtreat_reference as ref_runs # noqa: E402

no_tb = dict(beta=0, trans_asymp=0, inf_cle=0, inf_non=0, inf_asy=0, k_latent=0, inf_non_remote=0, inf_asy_remote=0, non_rec=0, non_asy=0, asy_non=0,
             asy_sym=0, sym_asy=0, mu=0, sym_dead=0.0) # Only treatment moves people


def test_defaults_reproduce_reference():
    """ With exp 35's options at their defaults, TBTreat reproduces its runs from before them (make_tbtreat_reference.py), step for step """
    ref = np.load(sc.thispath(__file__) / 'data' / 'tbtreat_reference.npz')
    for name, c in ref_runs.cases.items():
        got = ref_runs.results(ref_runs.make(**c).run())
        keys = [k.split('/', 1)[1] for k in ref.files if k.startswith(name + '/')]
        assert set(keys) == set(got), (name, set(keys) ^ set(got))
        for k in keys:
            assert np.max(np.abs(got[k] - ref[f'{name}/{k}'])) < 1e-12, (name, k)


def cohort(solver, dt, n=200_000, years=3):
    """ A cohort with nothing but treatment, from symptomatic TB, asymptomatic TB (its own death rates), latent infection and susceptible, with
    every exit in use """
    widths = tbode4.phase_lengths(0.5, 1/6, dt.years)
    tb = tbode4.TBTreat(init=dict(symptomatic=0.4, asymptomatic=0.2, infection=0.2), tx_from={'symptomatic': 4.0, 'asymptomatic': 2.0, 'infection': 2.0, 'susceptible': 0.5},
                        tx_success=0.6, tx_clear_latent=0.9, tx_p_relapse=0.1, dur_relapse_tx=100, tx_phase_dur=widths,
                        tx_death_fail={'symptomatic': 0.5, 'asymptomatic': 0.2}, tx_death_cure={'symptomatic': [0.4, 0, 0]}, tx_ltfu=0.6, tx_ltfu_cure=lambda f: f, **no_tb)
    kw = dict(n_agents=n, rand_seed=3, analyzers=[tbode4.TxChecks()]) if solver == 'abm' else dict(total_pop=n)
    return ss.Sim(tb, solver=solver, start=2000, stop=2000 + years, dt=dt, verbose=0, **kw).run()


def test_cohort_agents_match_ode():
    """ Every treatment transition's total over 3 years, agents against v4's ODE, at monthly and weekly steps with phases of whole steps (within
    4 SE); and agents' checks: nobody treated without TB disease is cured or dies of TB, and those leaving treatment are in their starting state """
    for dt in [ss.months(1), ss.weeks(1)]:
        ode, abm = cohort('ode', dt), cohort('abm', dt)
        tb = ode.diseases.tb
        names = [k for k in tb.transitions if k.startswith(('tx_', 'tb_death_t')) and k in ode.results.tb]
        for k in names:
            o, a = np.asarray(ode.results.tb[k].values).sum(), np.asarray(abm.results.tb[k].values).sum()
            if o > 100:
                assert abs(a/o - 1) < 4/np.sqrt(o) + 0.01, (str(dt), k, a/o, o)
        c = abm.analyzers[0].c
        assert c.notb_cured == 0 and c.tb_death_notb == 0 and c.relapse_not_from_success == 0, dict(c)
        assert c.latent_to_relapsing == 0 and c.tb_death_latent == 0 and c.latent_cleared > 0, dict(c)
        assert not any(k.startswith('tb_death_treatment_success_asymptomatic') for k in names) # Asymptomatic TB: deaths on the fail track only, as given
        assert c.left_elsewhere == 0, dict(c)


def test_outcome_shares_match_ode():
    """ v4's ODE: the shares of a cohort's courses cured, dying of TB, lost and failing, against outcome_shares() (exact in continuous time);
    daily steps, everyone starting treatment within the first week and nobody after (so failures and losses don't restart) """
    widths = tbode4.phase_lengths(0.5, 1/6, 1/12)
    mu_f, mu_c, h, q = 0.3, [0.2, 0, 0], 0.4, [1/6, 0.5, 5/6]
    tb = tbode4.TBTreat(init=dict(symptomatic=1.0), tx_rate={2000.0: 1000.0, 2000.02: 1000.0, 2000.0201: 0.0}, tx_success=0.6, tx_phase_dur=widths, tx_death_fail=mu_f, tx_death_cure=mu_c,
                        tx_ltfu=h, tx_ltfu_cure=q, dur_relapse_tx=1e6, **no_tb)
    sim = ss.Sim(tb, solver='ode', start=2000, stop=2002, dt=ss.days(1), total_pop=1, verbose=0).run()
    info, r = sim.diseases.tb.tx_info, sim.results.tb
    tot = lambda names: sum(float(np.asarray(r[k].values).sum()) for k in names)
    starts = tot(info.start)
    got = dict(cured=(tot(info.tx_cure) + tot(info.ltfu_cure))/starts, died_tb=tot(info.deaths)/starts, lost=(tot(info.ltfu_cure) + tot(info.ltfu_back))/starts,
               failed=tot(info.tx_fail)/starts)
    want = tbode4.outcome_shares(0.6, mu_f, mu_c, h, q, widths)
    for k, v in got.items():
        assert abs(v - want[k]) < 2e-3, (k, v, want[k]) # Restarts within the first week (0.6%) and phases of whole days


def test_carve_out_round_trip():
    """ carve_out() finds rates that reproduce the shares it was given; and caps c at 1 when overall cure is out of reach """
    widths = tbode4.phase_lengths(0.5, 1/6, 1/12)
    for om in [0.0, 0.5]:
        r = tbode4.carve_out(0.5, 0.04, 0.06, widths, q=[1/6, 0.5, 5/6], omega=om, bg=0.01)
        assert not r.capped and r.residual < 1e-10, r
    r = tbode4.carve_out(0.97, 0.04, 0.06, widths, q=0.0, bg=0.01)
    assert r.capped and r.c == 1.0 and r.c_raw > 1


def test_fastode_age_default_options():
    """ fastode_age's run_cohort with exp 35's options at their defaults (one phase, sym_dead on failing courses, no loss) equals its exp 23 treatment """
    from test_fastode_age import pars, init
    p = pars()
    a = fa.run_cohort(p, 1/12, 1950, 2000, table=init)
    b = fa.run_cohort(p, 1/12, 1950, 2000, table=init, tx=fa.tx_options(1/12, tx_dur=p['tx_dur']))
    for k in ['prevalence', 'incidence', 'mortality', 'notifications', 'on_treatment', 'alive']:
        assert np.max(np.abs(a[k] - b[k])) < 1e-12, k
