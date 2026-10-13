"""
Tests for tbode.TBExact: one natural-history step reproduces the exact transition matrix expm(Q*dt).
"""

import numpy as np
import pytest
import sciris as sc
import starsim as ss
import tbsim
from scipy.linalg import expm
from scipy.stats import chisquare
from tbsim import TBS
import tbode

DT = 1 / 12
N_AGENTS = 200_000
START_STATES = [TBS.INFECTION, TBS.NON_INFECTIOUS, TBS.ASYMPTOMATIC, TBS.SYMPTOMATIC]


def generator(p):
    """
    Generator (per year) built independently of tbode.exact, over
    [INFECTION, NON_INFECTIOUS, ASYMPTOMATIC, SYMPTOMATIC, CLEARED via INFECTION, CLEARED via NON_INFECTIOUS, DEAD].
    """
    Q = np.zeros((7, 7))
    Q[0, 4], Q[0, 1], Q[0, 2] = p.inf_cle, p.inf_non, p.inf_asy
    Q[1, 5], Q[1, 2] = p.non_rec, p.non_asy
    Q[2, 1], Q[2, 3] = p.asy_non, p.asy_sym
    Q[3, 2], Q[3, 6] = p.sym_asy, p.sym_dead
    Q[np.diag_indices_from(Q)] = -Q.sum(axis=1)
    return Q


def exact_step_matrix(p, dt):
    """Exact one-step matrix expm(Q·dt)."""
    return expm(generator(p) * dt)


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
    one_step.last_tb, one_step.last_uids = tb, uids
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
        pval = chisq_pvalue(obs, P[i] * m.sum())
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


class SeededExact(tbode.TBExact):
    """TBExact whose step infects ``seed_uids`` via set_prognoses in place of transmission, then steps as usual."""
    seed_uids = ss.uids()

    def step(self):
        self._new_inf = ss.uids()
        self.set_prognoses(self.seed_uids)
        self.step_transitions()
        self.step_bookkeeping()
        return


def new_infection_step(infection_timing, seed=1):
    """Infect every agent at the first step (all start SUSCEPTIBLE) and return their end-of-step categories."""
    p = tbode.TBODE().pars
    rates = {k: ss.peryear(p[k]) for k in ['inf_cle', 'inf_non', 'inf_asy', 'non_rec', 'non_asy', 'asy_non', 'asy_sym', 'sym_asy', 'sym_dead']}
    tb = SeededExact(beta=ss.peryear(0), init_prev=ss.bernoulli(0), infection_timing=infection_timing, **rates)
    sim = tbsim.Sim(sim_pars=dict(start=ss.years(0), stop=ss.years(1), dt=ss.years(DT), rand_seed=seed),
                    n_agents=N_AGENTS, tb_model=tb, networks=[], demographics=[], verbose=0)
    sim.init()
    tb = sim.get_tb()
    uids = sim.people.auids
    tb.seed_uids = uids
    sim.run_one_step()
    return p, observed_categories(np.asarray(tb.state[uids]), np.asarray(tb.rr_reinfection[uids]), p)


def uniform_time_matrix(p, dt):
    """∫₀¹ expm(Q·s·dt) ds, the one-step matrix averaged over a uniform infection time: the top-right block of expm([[Q·dt, I], [0, 0]])."""
    Q = generator(p)
    n = len(Q)
    A = np.zeros((2 * n, 2 * n))
    A[:n, :n] = Q * dt
    A[:n, n:] = np.eye(n)
    return expm(A)[:n, n:]


def chisq_pvalue(obs, exp):
    """Chi-square p-value, pooling categories with expected count < 5."""
    keep = exp >= 5
    obs_k = np.append(obs[keep], obs[~keep].sum())
    exp_k = np.append(exp[keep], exp[~keep].sum())
    if exp_k[-1] < 5:
        obs_k, exp_k = obs_k[:-1], exp_k[:-1]
        exp_k *= obs_k.sum() / exp_k.sum()
    return chisquare(obs_k, exp_k).pvalue


def test_within_step_infection_timing():
    """With infection_timing='within', agents infected this step end in each state with probability ∫₀¹ expm(Q·s·dt) ds (row INFECTION)."""
    p, cat = new_infection_step('within')
    assert np.all(cat >= 0)
    obs = np.bincount(cat, minlength=7)
    pval = chisq_pvalue(obs, uniform_time_matrix(p, DT)[0] * len(cat))
    assert pval > 1e-3, f'new-infection end states differ from the uniform-time average (p={pval:.2g})'


def test_start_timing_differs():
    """Control: with the default 'start' timing, new infections follow expm(Q·dt) instead, which the test above rejects."""
    p, cat = new_infection_step('start')
    obs = np.bincount(cat, minlength=7)
    assert chisq_pvalue(obs, exact_step_matrix(p, DT)[0] * len(cat)) > 1e-3
    assert chisq_pvalue(obs, uniform_time_matrix(p, DT)[0] * len(cat)) < 1e-6


K = 5.85  # tb_LAI_TPT's k_asy = k_non


def test_one_step_latent_classes():
    """With k > 0, INFECTION agents end the step with probability sum_m Poisson(m; r0/k) * expm(Q_m*dt)[INFECTION]."""
    p, start, end, rr = one_step(tbode.TBExact, k_asy=K, k_non=K)
    weights, rate_non, rate_asy = tbode.latent_classes(p.inf_non, p.inf_asy, K)
    P = np.zeros(7)
    for w, rn, ra in zip(weights, rate_non, rate_asy):
        P += w * exact_step_matrix(sc.mergedicts(p, dict(inf_non=rn, inf_asy=ra), _copy=True), DT)[0]
    cat = observed_categories(end, rr, p)
    m = start == TBS.INFECTION
    pval = chisq_pvalue(np.bincount(cat[m], minlength=7), P * m.sum())
    assert pval > 1e-3, f'INFECTION end states differ from the class mixture (p={pval:.2g})'


def test_rejects_unequal_k():
    """TBExact's latent classes need one decline rate for both progression routes."""
    with pytest.raises(NotImplementedError):
        one_step(tbode.TBExact, k_asy=K, k_non=0.0)


@pytest.mark.parametrize('infection_prob', ['exp', 'linear'])
def test_mass_action_infection_probability(infection_prob):
    """TBExactMassAction infects each susceptible with 1 - exp(-x) or x, x = beta*w*rel_sus*dt, w = (trans_asymp*ASY + SYM)/N."""
    p = tbode.TBODE().pars
    beta = 600.0  # Large, so x is big enough to tell the two forms apart
    tb = tbode.TBExactMassAction(beta=ss.peryear(0), init_prev=ss.bernoulli(0), beta_ode=beta, infection_prob=infection_prob)
    sim = tbsim.Sim(sim_pars=dict(start=ss.years(0), stop=ss.years(1), dt=ss.years(DT), rand_seed=3),
                    n_agents=N_AGENTS, tb_model=tb, networks=[], demographics=[], verbose=0)
    sim.init()
    tb = sim.get_tb()
    uids = sim.people.auids
    st = np.full(len(uids), TBS.SUSCEPTIBLE)
    st[:2000], st[2000:3000] = TBS.ASYMPTOMATIC, TBS.SYMPTOMATIC
    st[3000:103_000] = TBS.CLEARED
    tb.state[uids] = st
    tb.rr_reinfection[uids[3000:103_000]] = 0.5
    tb.step_bookkeeping()
    w = (p.trans_asymp * 2000 + 1000) / len(uids)
    x_s, x_c = beta * w * DT, beta * w * 0.5 * DT
    prob = (lambda x: 1 - np.exp(-x)) if infection_prob == 'exp' else (lambda x: x)
    new, _, _ = tb.infect()
    is_new = np.isin(uids, new)
    for mask, x in [(st == TBS.SUSCEPTIBLE, x_s), (st == TBS.CLEARED, x_c)]:
        n, k = mask.sum(), is_new[mask].sum()
        se = np.sqrt(n * prob(x) * (1 - prob(x)))
        assert abs(k - n * prob(x)) < 4 * se, f'{infection_prob}: {k} infections vs {n * prob(x):.0f} expected'
    assert not is_new[st == TBS.ASYMPTOMATIC].any()


def visit_probabilities(p, dt):
    """P(path visits ASYMPTOMATIC or SYMPTOMATIC within one step) from INFECTION and NON_INFECTIOUS, built independently:
    make ASY and SYM absorbing and read off the probability of having been absorbed into either."""
    Q = generator(p)
    Q[2:4, :] = 0  # ASYMPTOMATIC and SYMPTOMATIC absorbing: reaching them at any point within the step = visited
    P = expm(Q * dt)
    return P[0, 2] + P[0, 3], P[1, 2] + P[1, 3]


def test_event_incidence_counts_visits():
    """TBExact counts a new episode for INFECTION / NON_INFECTIOUS starters (episode flag off) with P(visit ASY/SYM within the step)."""
    p, start, end, rr = one_step(tbode.TBExact)
    sim_tb = one_step.last_tb
    flags = np.asarray(sim_tb.new_episode[one_step.last_uids])
    pv_inf, pv_non = visit_probabilities(p, DT)
    for s, pv in [(TBS.INFECTION, pv_inf), (TBS.NON_INFECTIOUS, pv_non)]:
        m = start == s
        n, k = m.sum(), flags[m].sum()
        se = np.sqrt(n * pv * (1 - pv))
        assert abs(k - n * pv) < 4 * se, f'{s.name}: {k} episodes counted vs {n * pv:.0f} expected'
        # Every agent that ends the step active was counted; some that end inactive were too (the snapshot misses these)
        ends_active = np.isin(end[m], [TBS.ASYMPTOMATIC, TBS.SYMPTOMATIC])
        assert flags[m][ends_active].all()
        assert flags[m][~ends_active].any()
    assert not flags[np.isin(start, [TBS.ASYMPTOMATIC, TBS.SYMPTOMATIC])].any()
    assert sim_tb.results.new_episodes[0] == flags.sum()
