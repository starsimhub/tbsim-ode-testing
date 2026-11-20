"""
Starting TB states by age and sex (tbode4.start; CONVENTIONS.md 10, exp 30): agents' drawn counts keep rare states in small groups, which Starsim's
rounding loses; and a table from settled_start() starts the ODE exactly at it. Rerun after bumping the v4 pins.

    PYTHONPATH=v4 v4/.venv/bin/python -m pytest v4/tests/test_start.py
"""
import numpy as np
import pandas as pd
import starsim as ss
import tbsim
import tbode4

ages = tbode4.age_edges(1.0)
pars = dict(beta=22.10, trans_asymp=0.7488, inf_cle=0.6638, inf_non=0.3153, inf_asy=0.3735, k_latent=7.517, inf_non_remote=0.0038, inf_asy_remote=0.0045,
            non_rec=0.2207, non_asy=0.2776, asy_non=1.612, asy_sym=0.7560, sym_asy=0.7593, sym_dead=0.3959, mu=0)


def infected_50plus(seed, n=20_000):
    """ Agents infected at 50+ at the start, with 0.3% infected at every age and agents stratified by single years """
    sim = ss.Sim(diseases=ss.SIS(beta=0, init=dict(I=ss.by(pd.DataFrame(dict(age=ages, value=0.003))))), people=ss.population('India', year=1950, n_agents=n),
                 total_pop=n, start=1950, stop=1951, solver='abm', verbose=0, rand_seed=seed)
    sim.stratify(age=ages)
    sim.init()
    u = np.asarray(sim.people.auids)
    expected = 0.003*np.sum(np.asarray(sim.people.age.raw[u], dtype=float) >= 50)
    infected = np.asarray(sim.diseases.sis.I.uids)
    return int(np.sum(np.asarray(sim.people.age.raw[infected], dtype=float) >= 50)), expected


def test_drawn_start():
    """ Rounded, a rare state all but vanishes in small old-age groups, in every seed; drawn, it's right on average """
    seeds = range(1, 11)
    try:
        tbode4.use_drawn_start(False)
        rounded = [infected_50plus(s)[0] for s in seeds]
        tbode4.use_drawn_start()
        tbode4.use_drawn_start() # Idempotent
        drawn = [infected_50plus(s) for s in seeds]
    finally:
        tbode4.use_drawn_start(False)
    n = np.array([d[0] for d in drawn])
    expected = np.mean([d[1] for d in drawn])
    assert np.mean(rounded) < 0.25*expected, (rounded, expected)
    assert abs(n.mean() - expected) < 3*np.sqrt(expected/len(n)), (n, expected)


def ode(init, width, years):
    sim = ss.Sim(tbsim.TB(init=init, **pars), solver='ode', people=ss.population('India', year=1950), total_pop=1.0, start=1950, stop=1950 + years,
                 dt=ss.months(1), verbose=0)
    e = tbode4.age_edges(width)
    sim.stratify(age=e, sex=tbode4.joint_sex(e, year=1950))
    return sim.run()


def test_settled_start():
    """ The table from a quarter-year ODE, averaged to single years, starts a single-year ODE exactly at it """
    run = ode(dict(infection=0.05, infection_remote=0.3, cleared=0.2, non_infectious=0.002, asymptomatic=0.002, symptomatic=0.002), 0.25, 5)
    st = tbode4.settled_start(run, drift_years=2)
    assert set(st.init) <= set(tbode4.tb_states) - {'susceptible'}
    assert st.drift['infection'] > 0 # Not settled after 5 years: the drift shows it
    sim = ode(st.init, 1.0, 1)
    sh, lows, _ = tbode4.shares_by_group(sim, 0)
    assert np.array_equal(lows, st.edges)
    for k, v in st.table.items():
        assert np.max(np.abs(sh[k] - v)) < 1e-12, k
