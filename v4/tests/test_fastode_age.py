"""
tbode4.fastode_age, the standalone ODE with ages and demography (exp 31): its copy of v4's ODE (run_rate) reproduces v4's quarter-year ODE with
TB, treatment and India's WPP demography; and its exact aging by cohort (run_cohort) at monthly steps is within 0.2% of half-month steps.
Rerun after bumping the v4 pins (the first test runs v4's ODE itself).

    PYTHONPATH=v4 v4/.venv/bin/python -m pytest v4/tests/test_fastode_age.py
"""
import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss
import tbode4
from tbode4 import fastode as fo, fastode_age as fa

start, stop = 1950, 2041
pi = 0.3518
nat = dict(beta=22.10, trans_asymp=0.7488, rr_reinfection_cleared=pi, rr_reinfection_rec=pi, rr_reinfection_treat=pi, rr_reinfection_inf=pi, inf_cle=0.6638,
           inf_non=0.3153, inf_asy=0.3735, k_latent=7.517, inf_non_remote=0.0038, inf_asy_remote=0.0045, non_rec=0.2207, non_asy=0.2776, asy_non=1.612,
           asy_sym=0.7560, sym_asy=0.7593) # Exp 21's set 1
trend = sc.objdict(sym_dead_1995=0.3959, sym_dead_2024=0.1426, tx_covg_2010=0.5479, tx_covg_2024=0.6772, efficacy_1995=0.4525, efficacy_2024=0.5300,
                   p_relapse=0.0339, over_detection=0.1831) # tb_LAI_TPT's best-fit treatment
init = dict(infection=0.012, infection_remote=0.25, cleared=0.2, recovered=0.03, non_infectious=0.003, asymptomatic=0.002, symptomatic=0.001) # The same at every age
ss.solvers.Engine.group_bytes = 4e9 # Quarter-year groups keep their results at every step (Starsim keeps them only every few steps above 100 MB per run)
keys = ['alive', 'prevalence', 'mortality', 'notifications'] + [f'prev_{g}' for g in fa.groups5] + [f'share_{g}' for g in fa.groups5]


def v4_ode():
    """ v4's ODE as exp 30 runs it: TBTreat with exp 24's schedules, WPP demography, quarter-year groups split jointly by sex, monthly steps """
    tx = tbode4.lai_tpt_treatment(row=dict(tx_covg_2010=trend.tx_covg_2010, tx_covg_2024=trend.tx_covg_2024, efficacy_1995=trend.efficacy_1995,
                                           efficacy_2024=trend.efficacy_2024, p_relapse=trend.p_relapse), step=1/12, end=stop)
    sym_dead = ss.peryear(ss.timeseries({1995.0: trend.sym_dead_1995, 2024.0: trend.sym_dead_2024}, interpolation='linear'))
    tb = tbode4.TBTreat(init=init, course='fixed', tx_rate=tx.tx_rate, tx_success=tx.tx_success, tx_p_relapse=trend.p_relapse, dur_relapse_tx=1.5,
                        sym_dead=sym_dead, mu=0, **nat)
    sim = ss.Sim(tb, demographics=tbode4.ode_demographics(tbode4.WPP()), people=ss.population('India', year=start), total_pop=1.0, start=start, stop=stop,
                 dt=ss.months(1), solver='ode', verbose=0)
    e = tbode4.age_edges(0.25)
    sim.stratify(age=e, sex=tbode4.joint_sex(e, year=start))
    sim.run()
    out = tbode4.tx_outputs(sim)
    r = sim.results.tb
    def by_age(res):
        b = res.by('age')
        return np.asarray(b.values, dtype=float), np.array(ss.data.age_edges(b.columns), dtype=float)
    pop = sum(by_age(r[k])[0] for k in tbode4.tb_states if k in r)
    act = by_age(r.asymptomatic)[0] + by_age(r.symptomatic)[0]
    lows = by_age(r.susceptible)[1]
    g = np.clip(np.digitize(lows + 1e-9, [lo for lo, _ in fa.groups5.values()]) - 1, 0, len(fa.groups5) - 1)
    a5 = np.stack([pop[:, g == k].sum(axis=1) for k in range(len(fa.groups5))], axis=1)
    p5 = np.stack([act[:, g == k].sum(axis=1) for k in range(len(fa.groups5))], axis=1)
    alive = np.asarray(sim.results.alive.values, dtype=float)
    out.alive = alive/alive[0]
    for k, name in enumerate(fa.groups5):
        out[f'prev_{name}'] = p5[:, k]/a5[:, k]
        out[f'share_{name}'] = a5[:, k]/a5.sum(axis=1)
    return yearly(out, 1/12)


def yearly(out, dt):
    """ Yearly means of `keys` as fa.yearly() takes them """
    year = np.floor(np.asarray(out.year)[1:] - dt/2 + 1e-9).astype(int)
    df = pd.DataFrame({k: (out[k][1:] if k in fa.stock_keys else out[k][:-1]) for k in keys})
    df['year'] = year
    return df.groupby('year').mean()


def pars():
    return fo.make_pars(**nat, mu=0, **fa.schedules(trend))


def test_copy_matches_v4_ode():
    """ run_rate (continuous time, 12 Erlang stages for the course) against v4's ODE (RK4, monthly steps, its exact fixed course): within 0.1% every
    year before treatment, and 0.2% from 1995, where the Erlang stages stand in for the fixed course (exp 31) """
    ref = v4_ode()
    out = fa.run_rate(pars(), 0.25, start, stop, 1/12, table=init, rtol=1e-6, atol=1e-12)
    got = yearly(out, 1/12)
    yrs = got.index[(got.index >= 1951) & (got.index < stop)]
    for k in keys:
        r = got.loc[yrs, k]/ref.loc[yrs, k]
        if k == 'notifications':
            r = r[r.index >= 1996]
        before, after = r[r.index < 1995], r[r.index >= 1995]
        if len(before):
            assert np.nanmax(np.abs(before - 1)) < 1e-3, (k, before.describe())
        assert np.nanmax(np.abs(after - 1)) < 2e-3, (k, after.describe())


def test_cohort_converges():
    """ run_cohort at monthly steps against half-month steps (mid-step force of infection): within 0.2% every year, and nearer than quarter-year steps """
    p = pars()
    runs = {h: yearly(fa.run_cohort(p, h, start, stop, table=init), h) for h in [0.25, 1/12, 1/24]}
    yrs = runs[1/24].index[(runs[1/24].index >= 1951) & (runs[1/24].index < stop)]
    for k in keys:
        dev = {h: np.nanmax(np.abs(runs[h].loc[yrs, k]/runs[1/24].loc[yrs, k] - 1)) for h in [0.25, 1/12]}
        assert dev[1/12] < 2e-3, (k, dev)
        assert dev[1/12] <= dev[0.25] + 1e-6, (k, dev)


def test_settle_round_trip():
    """ start_table() of a settled run, put back on the same grid by initial(), keeps each single year's shares """
    p = pars()
    st = fa.settle(p, 0.25, init, 0.5, years=30, drift_years=5)
    table, flagged = fa.start_table(st)
    X = fa.initial(st.layout, st.lows, table, flagged, pop=st.X[:st.layout.ns].sum(axis=0))
    back, flagged2 = fa.start_table(sc.objdict(X=X, lows=st.lows, layout=st.layout))
    for k in table:
        assert np.allclose(back[k], table[k], atol=1e-12), k
    assert abs(flagged2 - flagged) < 1e-12
