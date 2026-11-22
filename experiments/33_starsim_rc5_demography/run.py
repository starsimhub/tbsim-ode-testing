"""
Experiment 33: moving the v4 pins to Starsim 4.0.0rc5 (0089ac0f) and tbsim a2ae583. See README.md. Run from the repo root.

Stage 1 (the bump):

    PYTHONPATH=v4 experiments/33_starsim_rc5_demography/.venv/bin/python experiments/33_starsim_rc5_demography/run.py tests     # v4/tests/ at the new pins
    PYTHONPATH=v4 v4/.venv/bin/python experiments/33_starsim_rc5_demography/run.py pins_ode                                     # exp 30's ODE, current pins
    PYTHONPATH=v4 experiments/33_starsim_rc5_demography/.venv/bin/python experiments/33_starsim_rc5_demography/run.py pins_ode  # the same, new pins
    PYTHONPATH=v4 experiments/33_starsim_rc5_demography/.venv/bin/python experiments/33_starsim_rc5_demography/run.py pins_compare

`pins_ode` runs v4's ODE as exp 30 does (TBTreat with exp 24's schedules, India's WPP demography, quarter-year groups split jointly by sex,
monthly steps, 1950–2041) with births' female share fixed at the 1950 population's overall share, which the current pins use whatever
`ss.birth`'s p_female says, so that both pins run the same model. It keeps every step's states and flows, so a change in how flows are dated
would show (equilibrium checks can't see one).

Stage 2 (demography alone, exp 29's setting), at the new pins, on 4 cores:

    PYTHONPATH=v4 experiments/33_starsim_rc5_demography/.venv/bin/python experiments/33_starsim_rc5_demography/run.py demog
    PYTHONPATH=v4 experiments/33_starsim_rc5_demography/.venv/bin/python experiments/33_starsim_rc5_demography/run.py demog_summary

Stage 3 (TB with demography, a reduced rerun of exp 30), at the new pins (EXP33_NCPUS for the agents' workers):

    PYTHONPATH=v4 experiments/33_starsim_rc5_demography/.venv/bin/python experiments/33_starsim_rc5_demography/run.py tb
    PYTHONPATH=v4 experiments/33_starsim_rc5_demography/.venv/bin/python experiments/33_starsim_rc5_demography/run.py tb_summary

Stage 4 (births' female share in fastode_age, and its effect on exp 32's targets), at the new pins:

    PYTHONPATH=v4 experiments/33_starsim_rc5_demography/.venv/bin/python experiments/33_starsim_rc5_demography/run.py fastode_age
"""
import os
os.environ.setdefault('NUMBA_NUM_THREADS', '1')
import sys
import subprocess
import warnings
import importlib.util
import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss
import tbode4
from tbode4 import demog

here = sc.thispath(__file__)
root = here.parents[1]
start, stop = 1950, 2041
p_female_1950 = demog.female_share('India', start)[2] # 0.48635: the current pins' births, whatever ss.birth's p_female
nat = dict(beta=22.10, trans_asymp=0.7488, rr_reinfection_cleared=0.3518, rr_reinfection_rec=0.3518, rr_reinfection_treat=0.3518, rr_reinfection_inf=0.3518,
           inf_cle=0.6638, inf_non=0.3153, inf_asy=0.3735, k_latent=7.517, inf_non_remote=0.01194*0.3153, inf_asy_remote=0.01194*0.3735, non_rec=0.2207,
           non_asy=0.2776, asy_non=1.612, asy_sym=0.7560, sym_asy=0.7593) # Exp 24's re-identification truth (exp 21's set 1), as v4/tests/test_fastode_age.py
trend = sc.objdict(sym_dead_1995=0.3959, sym_dead_2024=0.1426, tx_covg_2010=0.5479, tx_covg_2024=0.6772, efficacy_1995=0.4525, efficacy_2024=0.5300, p_relapse=0.0339)
init = dict(infection=0.012, infection_remote=0.25, cleared=0.2, recovered=0.03, non_infectious=0.003, asymptomatic=0.002, symptomatic=0.001) # The same at every age


def env_label():
    return f'starsim_{ss.__version__}'


def ode_births(wpp, p_female):
    """ tbode4.ode_births with an explicit p_female (honoured by the new pins' ODE; the current pins' ODE takes its starting share instead) """
    zero = ss.peryear(0.0)
    fert = {a: demog.wpp_rate(wpp, wpp.asfr[:, j]) for j, a in enumerate(wpp.fert_ages)}
    rows = [dict(age=a, sex='female', value=fert.get(a, zero)) for a in demog.wpp_edges] + [dict(age=a, sex='male', value=zero) for a in demog.wpp_edges]
    return ss.birth(rate=ss.by(pd.DataFrame(rows)), p_female=p_female)


def pins_ode():
    """ Exp 30's ODE at whichever pins this environment has; every step's states (per person) and flows saved """
    ss.solvers.Engine.group_bytes = 4e9
    wpp = tbode4.WPP()
    tx = tbode4.lai_tpt_treatment(row=dict(tx_covg_2010=trend.tx_covg_2010, tx_covg_2024=trend.tx_covg_2024, efficacy_1995=trend.efficacy_1995,
                                           efficacy_2024=trend.efficacy_2024, p_relapse=trend.p_relapse), step=1/12, end=stop)
    sym_dead = ss.peryear(ss.timeseries({1995.0: trend.sym_dead_1995, 2024.0: trend.sym_dead_2024}, interpolation='linear'))
    tb = tbode4.TBTreat(init=init, course='fixed', tx_rate=tx.tx_rate, tx_success=tx.tx_success, tx_p_relapse=trend.p_relapse, dur_relapse_tx=1.5,
                        sym_dead=sym_dead, mu=0, **nat)
    sim = ss.Sim(tb, demographics=[ode_births(wpp, p_female_1950), demog.death_rate(wpp)], people=ss.population('India', year=start), total_pop=1.0,
                 start=start, stop=stop, dt=ss.months(1), solver='ode', verbose=0)
    e = tbode4.age_edges(0.25)
    sim.stratify(age=e, sex=tbode4.joint_sex(e, year=start))
    T = sc.timer()
    sim.run()
    r = sim.results
    alive = np.asarray(r.alive.values, dtype=float)
    out = sc.objdict(year=np.asarray(sim.t.yearvec, dtype=float), alive=alive, secs=T.toc(output=True), env=env_label())
    for k in tbode4.tb_states:
        if k in r.tb:
            out[f'state_{k}'] = np.asarray(r.tb[k].values, dtype=float)/alive
    tx_out = tbode4.tx_outputs(sim)
    for k in ['notifications', 'mortality', 'relapse', 'on_treatment']:
        if k in tx_out:
            out[f'out_{k}'] = np.asarray(tx_out[k], dtype=float)
    for k in ['births', 'died']: # Population flows
        if k in r:
            out[f'flow_{k}'] = np.asarray(r[k].values, dtype=float)
    sc.save(here / f'pins_ode_{env_label()}.obj', out)
    out.versions = tbode4.versions('v4') # Each pin's versions with its run (versions.json holds the last environment's)
    tbode4.save_versions(here, 'v4')
    print(f'{env_label()}: {out.secs:.1f} s; {len([k for k in out if k.startswith(("state_", "out_", "flow_"))])} series', flush=True)
    return out


def pins_compare():
    """ Step-by-step comparison of pins_ode's runs: the largest relative difference in each series """
    files = sorted(here.glob('pins_ode_starsim_*.obj'))
    assert len(files) == 2, f'Expected a run at each pin, found {[f.name for f in files]}'
    a, b = [sc.load(f) for f in files]
    lines = [f'{a.env} against {b.env} (exp 30\'s ODE, monthly, 1950–2041); run time {a.secs:.1f} s against {b.secs:.1f} s',
             f'Steps: {len(a.year)} against {len(b.year)}; year grids equal: {np.allclose(a.year, b.year) if len(a.year) == len(b.year) else False}']
    rows = []
    for k in [k for k in a if k.startswith(('state_', 'out_', 'flow_')) or k == 'alive']:
        if k not in b or np.shape(a[k]) != np.shape(b[k]):
            rows.append(dict(series=k, note='missing or a different shape'))
            continue
        x, y = np.asarray(a[k]), np.asarray(b[k])
        scale = np.maximum(np.abs(x), 1e-300)
        rel = np.abs(y - x)/scale
        ok = np.abs(x) > 1e-12*np.nanmax(np.abs(x)) # Ignore exact zeros (e.g. treatment before 1995)
        row = dict(series=k, max_rel=float(np.nanmax(rel[ok])) if ok.any() else 0.0)
        if len(x) > 2: # A shift by one step would match the neighbouring index better
            row['max_rel_shifted_by_one'] = float(np.nanmax((np.abs(y[1:] - x[:-1])/scale[:-1])[ok[:-1]])) if ok[:-1].any() else 0.0
        rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(here / 'pins_compare.csv', index=False)
    with pd.option_context('display.width', 200):
        lines.append(df.to_string(index=False, float_format=lambda v: f'{v:.2e}'))
    text = '\n'.join(lines)
    print(text)
    sc.savetext(here / 'pins_compare.log', text)
    return df


def tests():
    """ v4/tests/ in this environment (the trial environment at the new pins), with the log saved """
    cmd = [sys.executable, '-m', 'pytest', '-q', '-p', 'no:cacheprovider', str(root / 'v4' / 'tests')]
    res = subprocess.run(cmd, capture_output=True, text=True, env={**os.environ, 'PYTHONPATH': str(root / 'v4')}, cwd=root)
    sc.savetext(here / f'tests_{env_label()}.log', res.stdout + res.stderr)
    print(res.stdout[-3000:], flush=True)
    return res.returncode


#%% Stage 2: demography alone (exp 29's setting), Starsim's defaults against exp 29's corrected setup

def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


r28 = load_module('run28', root / 'experiments' / '28_age_demography' / 'run.py') # Exp 28's WPP rates, deaths, annual ASFR table and helpers for results by age
wpp = r28.wpp # Exp 28's WPP (annual ASFR table, realized fertility)
wpp4 = tbode4.WPP() # tbode4's, for its births and the shifted table
d_start, d_stop = 1950, 2051 # Exp 29's years
n_agents, reps = 100_000, 10
ncpus = int(os.environ.get('EXP33_NCPUS', tbode4.ncpus))
arms = sc.objdict(Q='ode', E='ode', Qj='ode', D='abm', W='abm')
widths = dict(Q=0.25, E=0.125, Qj=0.25)
report_years = [1950, 1951, 1955, 1960, 1970, 1980, 1990, 2000, 2010, 2020, 2030, 2040, 2050]
decades_fert = [1950, 1960, 1980, 2000, 2020, 2040]
fn_demog = here / 'demog.obj'


def starting_people(arm):
    """ D and the ODEs: ss.population (female by its fraction at each age). W: exp 29's P2b start, WPP's 1950 ages with one female share (the overall) """
    if arm != 'W':
        return ss.population('India', year=d_start, n_agents=n_agents if arms[arm] == 'abm' else None)
    pyr = ss.data.population('India', year=d_start)
    ages = np.array(ss.data.age_edges(pyr.index), dtype=float)
    return ss.People(n_agents=n_agents, age_data=np.column_stack([ages, (pyr.female + pyr.male).values]), p_female=p_female_1950)


def make_demog_sim(arm, rand_seed=0):
    abm = arms[arm] == 'abm'
    if arm in ['Q', 'E']: # Starsim's defaults: the ODE split by sex as its population is at each age, births' female share its default
        demog_mods = [ode_births(wpp4, None), r28.death_rate()]
    elif arm == 'Qj': # Exp 29's R4x: the joint split by sex, births at the 1950 population's overall female share
        demog_mods = [ode_births(wpp4, p_female_1950), r28.death_rate()]
    elif arm == 'D': # Starsim's defaults: WPP's rates as given (by age at delivery), the default start of pregnancies and births' female share
        demog_mods = [ss.Pregnancy(fertility_rate=r28.annual_asfr(), rate_per=1000, min_age=10, max_age=55, name='pregnancy'), r28.death_rate()]
    else: # W, exp 29's P2b: the shifted table read at the age at conception (so not shifted twice); Starsim's start of pregnancies is now P2b's (h·d)
        tab, lo, hi = demog.shifted_fertility(wpp4)
        demog_mods = [ss.Pregnancy(fertility_rate=tab, rate_per=1000, min_age=lo, max_age=hi, fertility_age='conception', p_female=ss.bernoulli(p_female_1950),
                                   name='pregnancy'), r28.death_rate()]
    kw = dict(rand_seed=rand_seed) if abm else {}
    sim = ss.Sim(diseases=ss.SIS(beta=0, init_prev=0), demographics=demog_mods, people=starting_people(arm), total_pop=n_agents if abm else 1.0,
                 start=d_start, stop=d_stop, dt=ss.months(1), solver=arms[arm], verbose=0, **kw)
    if abm:
        sim.stratify(age=r28.edges, sex=r28.sexes) # Agents age exactly; their groups only report results
    else:
        e = tbode4.age_edges(widths[arm])
        sim.stratify(age=e, sex=tbode4.joint_sex(e, year=d_start) if arm == 'Qj' else r28.sexes)
    return sim


def result(r, *names):
    """ The first of `names` among the sim's results (Starsim renamed several: died -> death, births -> birth) """
    for path in names:
        node = r
        try:
            for part in path.split('.'):
                node = node[part]
            return np.asarray(node.values, dtype=float)
        except (KeyError, AttributeError):
            continue
    raise KeyError(f'None of {names} in the results')


def demog_outputs(sim):
    """ Exp 28's demog_outputs(), reading the new results' names: population by 5-year group and sex (1 January), crude rates, realized fertility """
    r = sim.results
    yv = np.asarray(sim.t.yearvec, dtype=float)
    years = np.arange(d_start, int(np.floor(yv[-1])) + 1)
    idx = np.array([r28.year_index(sim, y) for y in years])
    pyr = {s: r28.regroup(*r28.by_age_sex(r.sis.S, s), r28.pyr_edges)[idx] for s in r28.sexes}
    alive = np.asarray(r.alive.values, dtype=float)
    births = result(r, 'birth', 'sis.births', 'pregnancy.birth')
    deaths = result(r, 'death', 'died')
    rows = []
    for y in years[:-1]:
        m = (yv >= y - 1e-9) & (yv < y + 1 - 1e-9) # Flows indexed by their step's start
        rows.append(dict(year=y, cbr=1000*births[m].sum()/alive[m].mean(), cdr=1000*deaths[m].sum()/alive[m].mean()))
    total = pyr['female'].sum(axis=1) + pyr['male'].sum(axis=1)
    out = sc.objdict(years=years, total=total, pyr=pyr, alive=alive[idx], rates=pd.DataFrame(rows), asfr_real=None)
    if 'pregnancy' in r:
        fbins = [0] + list(wpp.fert_ages) + [55]
        dv = r28.regroup(*r28.by_age_sex(r.pregnancy.delivery, 'female'), fbins)[:, 1:-1]
        women = r28.regroup(*r28.by_age_sex(r.sis.S, 'female'), fbins)[:, 1:-1]
        real = []
        for y in years[:-1]:
            m = (yv >= y - 1e-9) & (yv < y + 1 - 1e-9)
            real.append(dv[m].sum(axis=0)/women[m].mean(axis=0))
        out.asfr_real = pd.DataFrame(real, index=years[:-1], columns=wpp.fert_ages)
    return out, births, deaths, alive, yv


def run_demog_one(arm, rep=0):
    T = sc.timer()
    if arm == 'E':
        ss.solvers.Engine.group_bytes = 2e9 # Keep results by group at every step
    sim = make_demog_sim(arm, rand_seed=33000 + 100*list(arms).index(arm) + rep)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        sim.run()
    out, births, deaths, alive, yv = demog_outputs(sim)
    r = sim.results
    out.warnings = sorted({str(w.message)[:160] for w in caught})
    gap = alive[1:] - alive[:-1] - (births[:-1] - deaths[:-1])
    out.balance = dict(max_step=float(np.max(np.abs(gap))), total=float(gap.sum()), births=float(births.sum()))
    fem0, all0 = r28.by_age_sex(r.sis.S, 'female')[0][1:, 0], r28.by_age_sex(r.sis.S)[0][1:, 0]
    out.birth_female = float(fem0.sum()/all0.sum()) # ODE: the female share of the youngest group (the newest births; the birth result isn't split by sex). Agents: replaced below
    out.start_female = float(r28.by_age_sex(r.sis.S, 'female')[0][0].sum()/r28.by_age_sex(r.sis.S)[0][0].sum())
    out.losses = None
    if 'pregnancy' in r: # Losses at the mother's death by balance, against ss.Pregnancy's own loss results
        rp = r.pregnancy
        P, conc, deliv = (np.asarray(rp[k].values, dtype=float) for k in ['P', 'conception', 'delivery'])
        rows = []
        for y in out.years[:-1]:
            m = (yv >= y - 1e-9) & (yv < y + 1 - 1e-9)
            i0, i1 = np.flatnonzero(m)[0], np.flatnonzero(m)[-1] + 1
            rec = sum(np.asarray(rp[k].values, dtype=float)[m].sum() for k in ['miscarriage', 'stillbirth', 'abortion', 'antepartum'] if k in rp)
            rows.append(dict(year=y, losses=float(P[i0] + conc[m].sum() - deliv[m].sum() - P[i1]), births=float(births[m].sum()), deliveries=float(deliv[m].sum()), recorded=float(rec)))
        out.losses = pd.DataFrame(rows)
        out.p_start = float(P[0])
        ppl = sim.people
        alive_uids = ppl.auids
        age = ppl.age.raw[alive_uids]
        newborn = ppl.parent.raw[alive_uids] != ppl.parent.nan
        gap_age = np.ceil(age[newborn]) - age[newborn]
        out.ages = dict(dtype=str(age.dtype), born_in_run=int(newborn.sum()), just_below_whole_year=int(np.sum((gap_age > 0) & (gap_age < 1e-6)))) # Exp 29: ±1.7% bands from these
        b_now = np.histogram(age[newborn], bins=list(r28.pyr_edges) + [200])[0] # Those born in the run, whose ages sit on the grid of steps
        b_off = np.histogram(age[newborn] + 1e-3, bins=list(r28.pyr_edges) + [200])[0]
        out.ages['groups_moved_by_1e-3'] = int(np.sum(b_now != b_off))
        kids = alive_uids[newborn]
        out.birth_female = float(ppl.female.raw[kids].mean()) # Female share of the newborns still alive at the end (births' share, less any difference in deaths)
    out.update(arm=arm, rep=rep, secs=T.toc(output=True), env=env_label())
    return out


def run_demog():
    T = sc.timer()
    data = sc.objdict()
    print(f'Mean gestation {tbode4.gestation():.4f} years', flush=True)
    for arm in [a for a, s in arms.items() if s == 'ode']:
        data[arm] = [run_demog_one(arm)]
        o = data[arm][0]
        print(f'{arm} done ({o.secs:.0f} s; female share at the start {o.start_female:.5f}, of births {o.birth_female:.5f})', flush=True)
    jobs = [dict(arm=a, rep=r) for a, s in arms.items() if s == 'abm' for r in range(reps)]
    res = sc.parallelize(run_demog_one, iterkwargs=jobs, ncpus=min(ncpus, len(jobs)))
    for arm in [a for a, s in arms.items() if s == 'abm']:
        data[arm] = [o for o in res if o.arm == arm]
        print(f'{arm}: {np.mean([o.secs for o in data[arm]]):.0f} s per run; female share at the start {np.mean([o.start_female for o in data[arm]]):.5f}, '
              f'of newborns alive at the end {np.mean([o.birth_female for o in data[arm]]):.5f}', flush=True)
    sc.save(fn_demog, data)
    tbode4.save_versions(here, 'v4')
    print(f'Stage 2 runs done ({T.toc(output=True):.0f} s)', flush=True)
    demog_summary()
    return data


def by_year(outs, key, years):
    """ [rep, year] of a yearly output: 'pop' (÷ 1950), 'cbr', 'cdr', or 'share_<group>' (as exp 29's) """
    rows = []
    for o in outs:
        yi = {int(y): i for i, y in enumerate(o.years)}
        if key == 'pop':
            v = o.total/o.total[0]
            rows.append([v[yi[y]] for y in years])
        elif key.startswith('share_'):
            sh = r28.shares5(o.pyr)
            rows.append([sh[yi[y], list(r28.groups5).index(key[6:])] for y in years])
        else:
            rt = o.rates.set_index('year')[key]
            rows.append([rt.get(y, np.nan) for y in years])
    return np.array(rows)


def demog_summary():
    """ Ratios between arms as mean ± SE over reps: D against the ODE extrapolated to exact aging (2E − Q), W against exp 29's P2b, Qj against exp 29's R4x """
    data = sc.load(fn_demog)
    e29 = sc.load(root / 'experiments' / '29_pregnancy_timing' / 'runs.obj')
    years = [y for y in report_years if y < d_stop]
    yrs_rates = [y for y in report_years if y < d_stop - 1]
    fmt = lambda x: [f'{m:.4f} ± {s:.4f}' if len(x) > 1 else f'{m:.6f}' for m, s in zip(x.mean(0), x.std(0, ddof=1)/np.sqrt(len(x)) if len(x) > 1 else np.zeros(x.shape[1]))]
    def ratios(key, yrs):
        v = {a: by_year(data[a], key, yrs) for a in arms}
        old = {a: by_year(e29[a], key, yrs) for a in ['R4x', 'R8x', 'P2b']}
        exact = 2*v['E'] - v['Q']
        out = {'D ÷ exact': fmt(v['D']/exact), 'D ÷ E': fmt(v['D']/v['E']), 'W ÷ exact (Qj-style births)': fmt(v['W']/(2*old['R8x'] - old['R4x'])),
               'W ÷ P2b (exp 29)': fmt(v['W']/old['P2b'].mean(0)), 'Qj ÷ R4x (exp 29)': fmt(v['Qj']/old['R4x']), 'Q ÷ Qj': fmt(v['Q']/v['Qj']),
               'D ÷ W': fmt(v['D']/v['W'].mean(0)), 'P2b ÷ exact (exp 29)': fmt(old['P2b']/(2*old['R8x'] - old['R4x']))}
        return pd.DataFrame(out, index=yrs), v
    lines = [f'Stage 2: demography alone, {data.D[0].env}. Agents: {reps} reps of {n_agents:,} at the start, monthly; mean ± SE over reps (ODEs deterministic). '
             '"exact": the ODE extrapolated to exact aging (2 × eighth-year − quarter-year)']
    pop, vpop = ratios('pop', years)
    lines += ['\nPopulation ÷ its 1950 value, ratios between arms:', pop.to_string()]
    lines += ['\nPopulation ÷ 1950 by arm (mean):', pd.DataFrame({a: vpop[a].mean(0) for a in arms}, index=years).round(4).to_string()]
    sh = []
    for g in r28.groups5:
        t, _ = ratios(f'share_{g}', years)
        t.insert(0, 'group', g)
        sh.append(t)
    shares = pd.concat(sh)
    lines += ['\nAge shares, ratios between arms:', shares[['group', 'D ÷ exact', 'W ÷ P2b (exp 29)', 'Qj ÷ R4x (exp 29)', 'Q ÷ Qj']].to_string()]
    cbr, _ = ratios('cbr', yrs_rates)
    lines += ['\nCrude birth rate, ratios between arms:', cbr[['D ÷ exact', 'D ÷ E', 'W ÷ P2b (exp 29)', 'Qj ÷ R4x (exp 29)', 'Q ÷ Qj']].to_string()]
    frows = []
    mids = np.array(wpp.fert_ages, dtype=float) + 2.5
    for arm in ['D', 'W']:
        for d in decades_fert:
            yrs = [y for y in range(d, d + 10) if y < d_stop - 1]
            real = np.mean([o.asfr_real.loc[yrs].values for o in data[arm]], axis=0)
            inp = np.array([wpp.asfr_at(y + 0.5) for y in yrs])
            row = dict(arm=arm, decade=d, **{f'{a}': float(real[:, j].sum()/inp[:, j].sum()) if inp[:, j].sum() > 0 else np.nan for j, a in enumerate(wpp.fert_ages)})
            row.update(TFR=float(real.sum()/inp.sum()), mac_real=float((real.mean(0)*mids).sum()/real.mean(0).sum()), mac_wpp=float((inp.mean(0)*mids).sum()/inp.mean(0).sum()))
            frows.append(row)
    ft = pd.DataFrame(frows)
    ft.to_csv(here / 'demog_fertility.csv', index=False)
    with pd.option_context('display.precision', 3, 'display.width', 200):
        lines += ['\nRealized fertility by age at delivery ÷ WPP\'s rate, by decade (mean over reps); TFR ratio; mean age at childbearing (realized, WPP):', ft.to_string(index=False)]
    sexes = pd.DataFrame([dict(arm=a, start_female=np.mean([o.start_female for o in data[a]]), birth_female=np.nanmean([o.birth_female for o in data[a]])) for a in arms])
    lines += ['\nFemale share at the start and of births (agents: of newborns alive at the end):', sexes.to_string(index=False, float_format=lambda v: f'{v:.5f}')]
    L = pd.concat([o.losses.assign(arm=a) for a in ['D', 'W'] for o in data[a]])
    lt = L.groupby('arm')[['losses', 'recorded', 'births']].sum()
    lines += ['\nPregnancies lost (by balance) against ss.Pregnancy\'s own loss results, summed over reps and years:', lt.to_string()]
    lines += ['\nAges of agents alive at the end (exp 29: 32-bit ages put whole-year ages a hair below, giving ±1.7% bands by group):',
              pd.DataFrame([dict(arm=a, rep=o.rep, **o.ages) for a in ['D', 'W'] for o in data[a]]).groupby('arm').agg('max').to_string()]
    warns = sorted({w for a in arms for o in data[a] for w in o.warnings})
    lines += ['\nWarnings raised during the runs:'] + [f'  {w}' for w in warns] if warns else ['\nWarnings raised during the runs: none']
    bal = pd.DataFrame([dict(arm=o.arm, rep=o.rep, **o.balance) for a in arms for o in data[a]])
    lines += ['\nBalance (population change − (births − deaths)), largest per step and total over the run:', bal.groupby('arm')[['max_step', 'total', 'births']].agg('max').to_string()]
    lines += ['\nRun time (s): ' + ', '.join(f'{a} {np.mean([o.secs for o in data[a]]):.0f}' for a in arms)]
    pd.concat([pop.add_prefix('pop '), cbr.add_prefix('cbr ')], axis=1).to_csv(here / 'demog_ratios.csv')
    shares.to_csv(here / 'demog_shares.csv')
    text = '\n'.join(lines)
    sc.savetext(here / 'demog_summary.log', text)
    print(text, flush=True)
    return


#%% Stage 3: TB with demography, a reduced rerun of exp 30 with Starsim's default demography

exp30 = root / 'experiments' / '30_tb_corrected_demography'
tb_sets, tb_dts, tb_widths, tb_reps = [0, 6], ['month', 'week'], ['quarter', 'eighth'], int(os.environ.get('EXP33_TB_REPS', 20))
fn_tb = sc.objdict(ode=here / 'tb_ode.obj', runs=here / 'tb_runs.obj')
_r30 = None


def r30():
    """ Exp 30's module, loaded on first use: its TB model, starting tables (start.obj), agents' drawn start (applied on import), AgeCounts and the ODE's
    first-entry incidence. Its demography (shifted table, PregnancyStart, joint_sex()) is not used here """
    global _r30
    if _r30 is None:
        _r30 = load_module('run30', exp30 / 'run.py')
    return _r30


def make_tb_sim(i, dt, solver, width='quarter', st=None, rand_seed=0, years=None):
    """ Exp 30's sim for set i with Starsim's default demography: agents with ss.Pregnancy on WPP's table as given (and every default); the ODE split by
    sex as its population is at each age, with births at their default female share; both from exp 30's starting table """
    m = r30()
    s = m.sets[i]
    abm = solver == 'abm'
    tb = m.tb_model(s, dt, st.init)
    kw = dict(people=ss.population('India', year=m.start, n_agents=m.n_agents if abm else None), total_pop=m.n_agents if abm else 1.0)
    if abm:
        kw.update(rand_seed=rand_seed, analyzers=[tbode4.Incidence(share_flagged=st.share_flagged), tbode4.TxChecks(), m.AgeCounts()])
        demog_mods = [ss.Pregnancy(fertility_rate=r28.annual_asfr(), rate_per=1000, min_age=10, max_age=55, name='pregnancy'), demog.death_rate(wpp4)]
    else:
        demog_mods = [ode_births(wpp4, None), demog.death_rate(wpp4)]
    t0, t1 = (m.start, m.stop) if years is None else years
    sim = ss.Sim(tb, demographics=demog_mods, start=t0, stop=t1, dt=m.dts[dt], solver=solver, verbose=0, **kw)
    if abm:
        sim.stratify(age=m.agent_edges, sex=tbode4.sexes)
    else:
        sim.stratify(age=tbode4.age_edges(m.widths[width]), sex=tbode4.sexes)
    return sim


def tb_seed(i, dt, rep):
    return 33500 + 100*(2*tb_sets.index(i) + tb_dts.index(dt)) + rep


def run_tb_agents(i, dt, rep, st):
    m = r30()
    T = sc.timer()
    sim = make_tb_sim(i, dt, 'abm', st=st, rand_seed=tb_seed(i, dt, rep))
    sim.run()
    inc, chk, ac = sim.analyzers[0], sim.analyzers[1], sim.analyzers[2]
    out = m.r28.tb_per_step(sim, inc)
    for k, g in enumerate(m.groups5):
        out[f'prev_{g}'] = ac.active[:, k]/ac.alive[:, k]
        out[f'share_{g}'] = ac.alive[:, k]/ac.alive.sum(axis=1)
    return sc.objdict(set=i, dt=dt, rep=rep, seed=tb_seed(i, dt, rep), yearly=m.r28.yearly(out, sim.t.dt_year), checks=dict(chk.c), pop_check=out.pop_check,
                      secs=T.toc(output=True))


def run_tb_ode(i, dt, w, st):
    m = r30()
    T = sc.timer()
    sim = make_tb_sim(i, dt, 'ode', width=w, st=st)
    sim.run()
    inc, _, _ = m.ode_incidence(sim, st.non0[w])
    out = m.r28.tb_per_step(sim, inc)
    return sc.objdict(set=i, dt=dt, width=w, yearly=m.r28.yearly(out, sim.t.dt_year), pop_check=out.pop_check, secs=T.toc(output=True))


def run_tb():
    T = sc.timer()
    m = r30()
    starts = sc.load(m.fn.start)
    ode = sc.objdict()
    for i in tb_sets:
        for dt in tb_dts:
            for w in tb_widths:
                ode[f'set{i}_{dt}_{w}'] = run_tb_ode(i, dt, w, starts[f'set{i}'])
                print(f'ODE set {i} {dt} {w}: {ode[f"set{i}_{dt}_{w}"].secs:.0f} s', flush=True)
    sc.save(fn_tb.ode, ode)
    jobs = [dict(i=i, dt=dt, rep=r, st=starts[f'set{i}']) for i in tb_sets for dt in tb_dts for r in range(tb_reps)]
    runs = sc.parallelize(run_tb_agents, iterkwargs=jobs, ncpus=min(ncpus, len(jobs)))
    sc.save(fn_tb.runs, runs)
    tbode4.save_versions(here, 'v4')
    print(f'Stage 3 runs done in {T.toc(output=True):.0f} s; agent runs {np.mean([o.secs for o in runs]):.0f} s each', flush=True)
    tb_summary()
    return runs


def tb_summary():
    """ Agents ÷ the ODE extrapolated to exact aging (2 × eighth − quarter), by window, side by side with exp 30's (its by_window.csv) """
    m = r30()
    runs, ode = sc.load(fn_tb.runs), sc.load(fn_tb.ode)
    exact = [sc.objdict(set=i, dt=dt, yearly=2*ode[f'set{i}_{dt}_eighth'].yearly - ode[f'set{i}_{dt}_quarter'].yearly) for i in tb_sets for dt in tb_dts]
    bw = m.r28.by_window(runs, exact)
    bw.to_csv(here / 'tb_by_window.csv', index=False)
    old = pd.read_csv(exp30 / 'by_window.csv')
    main = ['prevalence', 'prop_asymptomatic', 'incidence_per_case', 'notifications_per_case', 'on_treatment_per_case', 'mortality_per_case',
            'recent', 'remote', 'cleared', 'alive'] + [f'prev_{g}' for g in m.groups5] + [f'share_{g}' for g in m.groups5]
    chk = pd.DataFrame([dict(set=r.set, dt=r.dt, **r.checks) for r in runs]).groupby(['set', 'dt']).sum()
    chk.to_csv(here / 'tb_pitfall_checks.csv')
    nrep = ', '.join(f'set {i} {dt}: {sum(1 for r in runs if r.set == i and r.dt == dt)}' for i in tb_sets for dt in tb_dts)
    lines = [f'Stage 3: agents ÷ the ODE with exact aging (2 × eighth − quarter), {env_label()}, Starsim\'s default demography in both solvers. '
             f'Mean ± SE over reps ({nrep}) of window means of yearly values; per_case = each ÷ prevalence. Windows: ' + ', '.join(f'{k} {a}–{b - 1}' for k, (a, b) in m.r28.windows.items())]
    with pd.option_context('display.width', 300, 'display.max_columns', 40):
        lines += ['\nThis experiment (new pins, no demography workarounds):', m.table(bw, main).to_string()]
        lines += ['\nExp 30 (current pins, with the workarounds), for comparison:', m.table(old, main).to_string()]
        lines += ['\nPitfall checks (summed over reps):', chk.to_string()]
    lines += ['\nODE run times (s): ' + ', '.join(f'{k} {o.secs:.0f}' for k, o in ode.items()),
              f'Agent run times: {np.mean([r.secs for r in runs if r.dt == "month"]):.0f} s monthly, {np.mean([r.secs for r in runs if r.dt == "week"]):.0f} s weekly']
    text = '\n'.join(lines)
    sc.savetext(here / 'tb_summary.log', text)
    print(text, flush=True)
    return bw


#%% Stage 4: births' female share in fastode_age, and its effect on exp 32's targets

def fastode_age_stage(n_draws=20, seed=33):
    """ (a) fastode_age with births at the 1950 share at ages 0–4 against v4's ODE with Starsim's default births (v4/tests' copy check);
    (b) exp 32's targets with births at that share ÷ at the overall share, at exp 32's best point and n_draws of its weighted posterior draws """
    from tbode4 import fastode_age as fa
    lines = []
    t = load_module('test_fastode_age', root / 'v4' / 'tests' / 'test_fastode_age.py')
    fa.Demography.birth_female_default = 'young'
    try:
        t.test_copy_matches_v4_ode()
        lines.append('(a) fastode_age (births at the ages 0–4 share) against v4\'s ODE with Starsim\'s default births: test_copy_matches_v4_ode passes')
    except AssertionError as e:
        lines.append(f'(a) test_copy_matches_v4_ode FAILS: {str(e)[:400]}')
    finally:
        fa.Demography.birth_female_default = 'start'
    _, share, overall = demog.female_share('India', 1950)
    lines.append(f'Births\' female share: overall 1950 population {overall:.5f}; ages 0–4 {share[0]:.5f}')
    c32 = load_module('calib32', root / 'experiments' / '32_demography_calibration' / 'calib.py')
    draws = pd.read_csv(root / 'experiments' / '32_demography_calibration' / 'real_draws.csv')
    rng = np.random.default_rng(seed)
    best = sc.load(root / 'experiments' / '32_demography_calibration' / 'real.obj')
    i_best = int(np.argmax(best.ll))
    thetas = [('best', best.theta[i_best])] + [(f'draw {int(r["index"])}', r[c32.names].values.astype(float)) for _, r in draws.iloc[rng.choice(len(draws), n_draws, replace=False)].iterrows()]
    keys = list(c32.targets) + ['prevalence_all_2020', 'population_2020']
    def run_all(bf):
        fa.Demography.birth_female_default = bf
        try:
            return [c32.simulate(th) for _, th in thetas] # Serially: about a minute, and the class-wide setting can't be copied stale into workers
        finally:
            fa.Demography.birth_female_default = 'start'
    old, new = run_all('start'), run_all('young')
    rows = [dict(point=name, **{k: new[i][k]/old[i][k] for k in keys}, ll_start=c32.loglik(old[i]), ll_young=c32.loglik(new[i])) for i, (name, _) in enumerate(thetas)]
    df = pd.DataFrame(rows)
    df.to_csv(here / 'birth_female_targets.csv', index=False)
    ratio_cols = [k for k in keys]
    summ = pd.DataFrame({'min': df[ratio_cols].min(), 'median': df[ratio_cols].median(), 'max': df[ratio_cols].max()})
    lines += [f'(b) Exp 32\'s targets, births at the ages 0–4 share ÷ at the overall share, at its best point and {n_draws} posterior draws (seed 1):',
              summ.to_string(float_format=lambda v: f'{v:.4f}'),
              f'Log-likelihood change (young − start): median {np.median(df.ll_young - df.ll_start):.3f}, range {np.min(df.ll_young - df.ll_start):.3f} to {np.max(df.ll_young - df.ll_start):.3f}']
    text = '\n'.join(lines)
    sc.savetext(here / 'fastode_age.log', text)
    print(text, flush=True)
    return df


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'tests'
    dict(tests=tests, pins_ode=pins_ode, pins_compare=pins_compare, demog=run_demog, demog_summary=demog_summary, fastode_age=fastode_age_stage, tb=run_tb, tb_summary=tb_summary)[which]()
