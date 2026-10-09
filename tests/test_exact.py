"""
Tests for tbode.TBExact: one natural-history step reproduces the exact transition matrix expm(Q*dt).
"""

import numpy as np
import pytest
import starsim as ss
import tbsim
from scipy.linalg import expm
from scipy.stats import chisquare
from tbsim import TBS
import tbode

DT = 1 / 12
N_AGENTS = 200_000
START_STATES = [TBS.INFECTION, TBS.NON_INFECTIOUS, TBS.ASYMPTOMATIC, TBS.SYMPTOMATIC]


def exact_step_matrix(p, dt):
    """
    Exact one-step matrix built independently of tbode.exact, over
    [INFECTION, NON_INFECTIOUS, ASYMPTOMATIC, SYMPTOMATIC, CLEARED via INFECTION, CLEARED via NON_INFECTIOUS, DEAD].
    """
    Q = np.zeros((7, 7))
    Q[0, 4], Q[0, 1], Q[0, 2] = p.inf_cle, p.inf_non, p.inf_asy
    Q[1, 5], Q[1, 2] = p.non_rec, p.non_asy
    Q[2, 1], Q[2, 3] = p.asy_non, p.asy_sym
    Q[3, 2], Q[3, 6] = p.sym_asy, p.sym_dead
    Q[np.diag_indices_from(Q)] = -Q.sum(axis=1)
    return expm(Q * dt)


def one_step(tb_class, seed=1, **tb_kwargs):
    """Spread agents evenly over the four transient states, take one step with no transmission or demography."""
    p = tbode.TBODE().pars
    rates = {k: ss.peryear(p[k]) for k in ['inf_cle', 'inf_non', 'inf_asy', 'non_rec', 'non_asy', 'asy_non', 'asy_sym', 'sym_asy', 'sym_dead']}
    tb = tb_class(beta=ss.peryear(0), init_prev=ss.bernoulli(0), **rates, **tb_kwargs)
    sim = tbsim.Sim(sim_pars=dict(start=ss.years(0), stop=ss.years(1), dt=ss.years(DT), rand_seed=seed),
                    n_agents=N_AGENTS, tb_model=tb, networks=[], demographics=[], verbose=0)
    sim.init()
    tb = sim.get_tb()
    uids = sim.people.auids
    start = np.array(START_STATES)[np.arange(len(uids)) % len(START_STATES)]
    tb.state[uids] = start
    tb.step_bookkeeping()
    sim.run_one_step()
    return p, start, np.asarray(tb.state[uids]), np.asarray(tb.rr_reinfection[uids])


def observed_categories(end, rr, p):
    """Map end states to the 7 categories of exact_step_matrix, using rr_reinfection to tell CLEARED pathways apart."""
    cat = np.full(len(end), -1)
    for i, s in enumerate(START_STATES):
        cat[end == s] = i
    cleared = end == TBS.CLEARED
    cat[cleared & np.isclose(rr, p.rr_reinfection_cleared)] = 4
    cat[cleared & np.isclose(rr, p.rr_reinfection_rec)] = 5
    cat[end == TBS.DEAD] = 6
    return cat


def test_one_step_matches_expm():
    """End-of-step state frequencies from each start state match expm(Q*dt), including the CLEARED pathway split."""
    p, start, end, rr = one_step(tbode.TBExact)
    P = exact_step_matrix(p, DT)
    cat = observed_categories(end, rr, p)
    assert np.all(cat >= 0), 'every agent should end in a recognised state'
    for i, s in enumerate(START_STATES):
        m = start == s
        obs = np.bincount(cat[m], minlength=7)
        exp = P[i] * m.sum()
        keep = exp >= 5  # pool rare categories so the chi-square approximation holds
        obs_k = np.append(obs[keep], obs[~keep].sum())
        exp_k = np.append(exp[keep], exp[~keep].sum())
        if exp_k[-1] < 5:
            obs_k, exp_k = obs_k[:-1], exp_k[:-1]
            exp_k *= obs_k.sum() / exp_k.sum()
        pval = chisquare(obs_k, exp_k).pvalue
        assert pval > 1e-3, f'{s.name}: end-state frequencies differ from expm(Q*dt) (p={pval:.2g})'


def test_tbsim_tb_differs():
    """Control: standard tbsim.TB fails the same test (forward cascade), so the test has power to detect the bias."""
    p, start, end, rr = one_step(tbsim.TB)
    P = exact_step_matrix(p, DT)
    cat = observed_categories(end, rr, p)
    m = start == TBS.NON_INFECTIOUS
    obs = np.bincount(cat[m], minlength=7)
    # NON_INFECTIOUS -> ASYMPTOMATIC is under-represented at the end of the step (exp 04: ratio ~0.90 at monthly dt)
    assert obs[2] / (P[1, 2] * m.sum()) < 0.95


@pytest.mark.parametrize('kwargs', [dict(k_asy=1.0), dict(k_non=1.0)])
def test_rejects_time_varying_hazards(kwargs):
    """TBExact refuses to run with time-since-infection progression, which it does not support."""
    with pytest.raises(NotImplementedError):
        one_step(tbode.TBExact, **kwargs)
