"""
Progression by time since infection (k_asy = k_non = k > 0): the Poisson-mixture latent classes in TBODE and
TBExact reproduce the exact cohort progression with hazards inf_non*exp(-k*tau), inf_asy*exp(-k*tau).
"""

import numpy as np
import sciris as sc
import tbode

K = 5.85  # tb_LAI_TPT's k_asy = k_non


def ode_cohort_from_tbode(pars, years=(1, 10), dt=1e-3):
    """
    Cohort progression computed with TBODE itself: everyone starts latent (split over the classes by their weights),
    no transmission, no background or TB death. Without reinfection, incidence (first ASY entry per episode) is the
    rate of first entry into ASYMPTOMATIC, so its integral is the cumulative fraction ever ASYMPTOMATIC.
    """
    ode = tbode.TBODE(**sc.mergedicts(pars, dict(beta=0, mu=0, sym_dead=0, k_asy=K, k_non=K)))
    y0 = np.zeros(ode.n_states + 1)
    y0[1:1 + ode.n_latent] = ode.latent_weights
    res = ode.run(start=0, stop=max(years), dt=dt, y0=y0)
    flow = res.incidence_py * res.N
    B = np.concatenate([[0], np.cumsum((flow[1:] + flow[:-1]) / 2 * dt)])
    return np.interp(years, res.time, B)


def test_latent_classes_converge():
    """The class weights sum to 1 and the truncation leaves only a negligible tail."""
    p = tbode.TBODE().pars
    w, rn, ra = tbode.latent_classes(p.inf_non, p.inf_asy, K)
    assert np.isclose(w.sum(), 1) and 3 <= len(w) <= 10
    assert np.allclose(rn + ra, np.arange(len(w)) * K)


def test_tbode_cohort_matches_exact():
    """TBODE's latent classes give the exact cohort progression (time-varying hazards integrated directly)."""
    p = tbode.TBODE().pars
    exact = tbode.ode_cohort_progression(p, k_asy=K, k_non=K).B
    assert np.allclose(ode_cohort_from_tbode(p), exact, rtol=1e-4)


def test_tbode_k0_unchanged():
    """k = 0 keeps a single INFECTION compartment, so the layout (and TB_ODE regression) is unchanged."""
    ode = tbode.TBODE()
    assert ode.states == tbode.states and ode.n_latent == 1


def test_equilibrium_with_classes():
    """The equilibrium with latent classes is stationary, and its INFECTION total is the sum of the classes."""
    ode = tbode.TBODE(beta=40, replace_tb_deaths=False, k_asy=K, k_non=K)
    x = ode.equilibrium()
    res = ode.run(start=0, stop=50, dt=1, y0=x * 1e5)
    assert res.prev_active[0] > 1e-4
    assert np.allclose(res.prev_active, res.prev_active[0], rtol=1e-6)
    assert np.isclose(res.INFECTION[0], x[1:1 + ode.n_latent].sum() * 1e5)


def test_tbexact_cohort_matches_exact():
    """TBExact's cohort progression (weekly dt, 200k agents) matches the exact curve within sampling error."""
    p = sc.mergedicts(tbode.TBODE().pars, dict(inf_non=0.6, inf_asy=0.6), _copy=True)  # Faster progression, for more events
    exact = tbode.ode_cohort_progression(p, k_asy=K, k_non=K).B
    abm = tbode.abm_cohort_progression(p, dt=7/365, n_agents=200_000, tb_class=tbode.TBExact, k_asy=K, k_non=K)
    B = abm.B[[0, 9]]
    se = np.sqrt(B * (1 - B) / 200_000)
    z = (B - exact) / se
    assert np.all(np.abs(z) < 4), f'TBExact cohort progression {B} vs exact {exact} (z = {z})'
