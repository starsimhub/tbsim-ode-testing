"""
Experiment 29, v4 environment: do agents' births match the ODE's when ss.Pregnancy uses fertility by age at delivery? See README.md.
Run from the repo root:

    PYTHONPATH=v4 v4/.venv/bin/python experiments/29_pregnancy_timing/run.py run      # every arm, then the summaries
    PYTHONPATH=v4 v4/.venv/bin/python experiments/29_pregnancy_timing/run.py more     # only the arms added later (P2, P3, R8), merged into the saved runs
    PYTHONPATH=v4 v4/.venv/bin/python experiments/29_pregnancy_timing/run.py summary  # summaries and figures from the saved runs

Environment variables: EXP29_NCPUS (cores for the agent runs; default tbode4.ncpus), EXP29_QUICK=1 (a small, short test, written to quick/).
"""
import os
os.environ.setdefault('NUMBA_NUM_THREADS', '1') # One thread per process; before importing starsim
import sys
import importlib.util
import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

here = sc.thispath(__file__)
root = here.parent.parent
sys.path.insert(0, str(root / 'v4'))
import tbode4 # noqa: E402

spec = importlib.util.spec_from_file_location('run28', root / 'experiments' / '28_age_demography' / 'run.py')
r28 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(r28) # Exp 28's WPP rates, deaths, ODE births and output helpers
wpp = r28.wpp

# Settings (README)
quick = os.environ.get('EXP29_QUICK') == '1'
ncpus = int(os.environ.get('EXP29_NCPUS', tbode4.ncpus))
outdir = here / 'quick' if quick else here
outdir.mkdir(exist_ok=True)
start, stop = 1950, 2051
n_agents, reps = 100_000, 10
report_years = [1950, 1951, 1955, 1960, 1970, 1980, 1990, 2000, 2010, 2020, 2030, 2040, 2050]
decades_fert = [1950, 1960, 1980, 2000, 2020, 2040]
if quick:
    stop, n_agents, reps = 1961, 5_000, 2
    report_years = [1950, 1951, 1955, 1960]
    decades_fert = [1950]
grid = 0.25 # Age grid of the shifted fertility table, in years
arms = sc.objdict(R='ode', P0='abm', P1='abm', P1w='abm', O1='ode', R2='ode', R4='ode', P2='abm', P3='abm', R8='ode', P3s='abm', P2s='abm', R4x='ode', R8x='ode', P2b='abm') # New arms at the end, so the others keep their seeds
later = ['R4x', 'R8x', 'P2b'] # The arms the 'more' step runs: the latest added (README)
dts = {a: ss.weeks(1) if a == 'P1w' else ss.months(1) for a in arms}
age_edges = {k: [float(a) for a in np.round(np.arange(0, 100.001, w), 4)] for k, w in dict(R2=0.5, R4=0.25, R8=0.125, R4x=0.25, R8x=0.125).items()} # Finer ODE age groups (the last, 100+, open): is R's aging converged?
fn = here / ('quick' if quick else '.') / 'runs.obj'


def gestation():
    """ ss.Pregnancy's mean duration of pregnancy, in years (0.769 at the default dur_pregnancy) """
    sim = ss.Sim(diseases=ss.SIS(beta=0, init_prev=0), demographics=[ss.Pregnancy(fertility_rate=30)], n_agents=100, dt=ss.months(1), start=2000, stop=2001, verbose=0)
    sim.init()
    return float(sim.demographics.pregnancy.mean_dur_pregnancy)


g = gestation()


def asfr_by_age(values, ages):
    """ WPP's rates (one decade's or year's, by 5-year group of age at delivery) at each age in `ages`, as a step function; 0 outside 10–54 """
    fa = np.array(wpp.fert_ages, dtype=float)
    j = np.clip(np.digitize(ages, fa) - 1, 0, len(fa) - 1)
    return np.where((ages >= 10) & (ages < 55), values[j], 0.0)


def shifted_table(offset_losses=False):
    """ P1's table for ss.Pregnancy: conception at age x takes WPP's rate for delivery age x + g, on a grid of `grid` years, every year (per 1,000 women)

    With offset_losses (P3), each rate is raised by exp(m·g), with m women's death rate (WPP, that year) at age x + g/2, mid-pregnancy: the inverse of
    the chance of surviving the pregnancy, so that births per woman-year, after pregnancies lost to the mother's death, are WPP's rate
    """
    ages = np.round(np.arange(10 - g, 55 - g, grid), 6)
    fem = r28.sexes.index('female')
    mid = np.clip(np.digitize(ages + g/2, r28.edges) - 1, 0, len(r28.edges) - 1) # WPP's mortality group of each age mid-pregnancy
    rows = []
    for y in np.arange(1950, 2101):
        vals = asfr_by_age(wpp.asfr_at(y), ages + g)
        if offset_losses:
            vals = vals*np.exp(wpp.mx_at(y)[fem][mid]*g)
        rows += [dict(Time=y, AgeGrp=a, ASFR=1000*v) for a, v in zip(ages, vals)]
    return pd.DataFrame(rows), float(ages[0]), float(ages[-1] + grid)


def wpp_female_share(year=start):
    """ WPP's population (as bundled) in `year`: the female share in each of its age groups (their lower edges), and overall """
    pyr = ss.data.population('India', year=year)
    lows = np.array(ss.data.age_edges(pyr.index), dtype=float)
    return lows, pyr.female.values/(pyr.female.values + pyr.male.values), float(pyr.female.sum()/(pyr.female.sum() + pyr.male.sum()))


def joint_sex(edges):
    """ R4x, R8x: the ODE's starting split by sex within each of its age groups, from WPP's 1950 population (each group takes the share of WPP's group holding it) """
    lows, share, _ = wpp_female_share()
    f = share[np.clip(np.digitize(edges, lows) - 1, 0, len(lows) - 1)]
    return ss.levels(['female', 'male'], fractions=ss.by(pd.DataFrame({'female': f, 'male': 1 - f}, index=edges), dims=['age', 'sex']))


class PregnancyStart(ss.Pregnancy):
    """ ss.Pregnancy with the steady-state share pregnant at the start (P2, P3): each woman pregnant with probability h·d (her conception rate × the
    mean gestation), instead of burnin_plan()'s 1 − exp(−h·steps·dt), with steps·dt the mean gestation rounded to whole steps """
    def burnin_plan(self, engine):
        plan = super().burnin_plan(engine)
        if plan is None:
            return plan
        orig, steps = plan.share, plan.steps
        factor = ss.solvers.convert(ss.years(self.mean_dur_pregnancy), engine.sim.unit_base)/engine.dt/steps # d ÷ (steps·dt)
        plan.share = lambda: np.minimum(-np.log1p(-np.minimum(orig(), 1 - 1e-12))*factor, 1.0) # orig = 1 − exp(−h·steps·dt) as agents (no losses before the start), so this is h·d
        return plan


def ode_births_older():
    """ O1's births: each single year of age [a, a+1) takes the mean, over the year of age, of WPP's rate at age − g (by decade, linear between) """
    sub = np.linspace(0, 1, 41)[:-1] + 1/80 # Midpoints of 40 parts of each year of age
    zero = ss.peryear(0.0)
    rows = []
    for a in r28.edges1:
        vals = np.array([asfr_by_age(wpp.asfr[d], a + sub - g).mean() for d in range(len(r28.decades))])
        rows.append(dict(age=a, sex='female', value=r28.rate(vals, False) if vals.any() else zero))
        rows.append(dict(age=a, sex='male', value=zero))
    return ss.birth(rate=ss.by(pd.DataFrame(rows)))


def make_sim(arm, rand_seed=0):
    solver = arms[arm]
    abm = solver == 'abm'
    kw = dict(people=ss.population('India', year=start, n_agents=n_agents if abm else None), total_pop=n_agents if abm else 1.0)
    if arm in ['P2s', 'P3s']: # The compartmental solver's start: WPP's 1950 ages, with half of each age women (Starsim's ODE takes a population's ages, not its sexes)
        pyr = ss.data.population('India', year=start)
        kw['people'] = ss.People(n_agents=n_agents, age_data=np.column_stack([np.array(ss.data.age_edges(pyr.index), dtype=float), (pyr.female + pyr.male).values]), p_female=0.5)
    if abm:
        kw['rand_seed'] = rand_seed
    if arm in ['R', 'R2', 'R4', 'R8', 'R4x', 'R8x']:
        demog = r28.ode_demographics()
    elif arm == 'O1':
        demog = [ode_births_older(), r28.death_rate()]
    elif arm == 'P0':
        demog = [ss.Pregnancy(fertility_rate=r28.annual_asfr(), min_age=10, max_age=55), r28.death_rate()]
    else:
        tab, lo, hi = shifted_table(offset_losses=arm in ['P3', 'P3s'])
        cls = PregnancyStart if arm in ['P2', 'P3', 'P2s', 'P3s', 'P2b'] else ss.Pregnancy
        extra = dict(p_female=ss.bernoulli(wpp_female_share()[2])) if arm == 'P2b' else {} # P2b: births as female as the ODE's with WPP's start (the starting share overall; Starsim's ODE ignores ss.birth's p_female)
        demog = [cls(fertility_rate=tab, min_age=lo, max_age=hi, name='pregnancy', **extra), r28.death_rate()]
    sim = ss.Sim(diseases=ss.SIS(beta=0, init_prev=0), demographics=demog, start=start, stop=stop, dt=dts[arm], solver=solver, verbose=0, **kw)
    sex = joint_sex(age_edges[arm]) if arm in ['R4x', 'R8x'] else r28.sexes
    sim.stratify(age=r28.edges if abm else age_edges.get(arm, r28.edges1), sex=sex) # Agents age exactly; their groups only report results
    return sim


def run_one(arm, rep=0):
    T = sc.timer()
    if arm in ['R8', 'R8x']:
        ss.solvers.Engine.group_bytes = 2e9 # Keep results by group at every step (above 100 MB per run, Starsim keeps them only every few steps)
    sim = make_sim(arm, rand_seed=29000 + 100*list(arms).index(arm) + rep)
    sim.run()
    out = r28.demog_outputs(sim)
    r = sim.results
    alive = np.asarray(r.alive.values, dtype=float)
    births = np.asarray(r.sis.births.values if 'births' in r.sis else r.pregnancy.birth.values, dtype=float)
    deaths = np.asarray(r.died.values, dtype=float)
    gap = alive[1:] - alive[:-1] - (births[:-1] - deaths[:-1]) # Flows at index i are the step from t_i to t_(i+1)
    bs = r.sis.births if 'births' in r.sis else r.pregnancy.birth
    out.birth_female = float(np.asarray(bs.select(sex='female').values, dtype=float).sum()/max(np.asarray(bs.values, dtype=float).sum(), 1e-12))
    out.start_female = float(r28.by_age_sex(r.sis.S, 'female')[0][0].sum()/r28.by_age_sex(r.sis.S)[0][0].sum())
    out.balance = dict(max_step=float(np.max(np.abs(gap))), total=float(gap.sum()), births=float(births.sum()))
    out.losses = None
    if 'pregnancy' in r: # Pregnancies lost to the mother's death, by balance: pregnant at the start + conceptions − deliveries − pregnant at the end (no abortions or spontaneous losses here); ss.Pregnancy's miscarriage and stillbirth results don't count these losses
        yv = np.asarray(sim.t.yearvec, dtype=float)
        rp = r.pregnancy
        P, conc, deliv = (np.asarray(rp[k].values, dtype=float) for k in ['P', 'conception', 'delivery'])
        rows = []
        for y in out.years[:-1]:
            m = (yv >= y - 1e-9) & (yv < y + 1 - 1e-9) # Flows indexed by their step's start
            i0, i1 = np.flatnonzero(m)[0], np.flatnonzero(m)[-1] + 1
            rows.append(dict(year=y, losses=float(P[i0] + conc[m].sum() - deliv[m].sum() - P[i1]), births=float(births[m].sum()), deliveries=float(deliv[m].sum()),
                             recorded=float(sum(np.asarray(rp[k].values, dtype=float)[m].sum() for k in ['miscarriage', 'stillbirth', 'abortion', 'antepartum']))))
        out.losses = pd.DataFrame(rows)
        out.p_start = float(P[0])
    out.update(arm=arm, rep=rep, secs=T.toc(output=True))
    return out


def run_all():
    T = sc.timer()
    print(f'Mean gestation g = {g:.4f} years', flush=True)
    data = sc.objdict()
    for arm in [a for a, s in arms.items() if s == 'ode']:
        data[arm] = [run_one(arm)]
        print(f'{arm} done ({data[arm][0].secs:.0f} s)', flush=True)
    jobs = [dict(arm=a, rep=r) for a, s in arms.items() if s == 'abm' for r in range(reps)]
    res = sc.parallelize(run_one, iterkwargs=jobs, ncpus=min(ncpus, len(jobs)))
    for arm in [a for a, s in arms.items() if s == 'abm']:
        data[arm] = [o for o in res if o.arm == arm]
        print(f'{arm}: {np.mean([o.secs for o in data[arm]]):.0f} s per run', flush=True)
    data.g = g
    sc.save(fn, data)
    tbode4.save_versions(outdir, 'v4')
    print(f'Runs done ({T.toc(output=True):.0f} s)', flush=True)
    summarize()
    return data


def run_more():
    """ Only the arms added later, merged into the saved runs (the others' seeds and results are unchanged) """
    T = sc.timer()
    data = sc.load(fn)
    for arm in [a for a in later if arms[a] == 'ode']:
        data[arm] = [run_one(arm)]
        print(f'{arm} done ({data[arm][0].secs:.0f} s; female share at the start {data[arm][0].start_female:.4f}, of births {data[arm][0].birth_female:.4f})', flush=True)
    jobs = [dict(arm=a, rep=r) for a in later if arms[a] == 'abm' for r in range(reps)]
    res = sc.parallelize(run_one, iterkwargs=jobs, ncpus=min(ncpus, len(jobs)))
    for arm in [a for a in later if arms[a] == 'abm']:
        data[arm] = [o for o in res if o.arm == arm]
        print(f'{arm}: female share at the start {np.mean([o.start_female for o in data[arm]]):.4f}, of births {np.mean([o.birth_female for o in data[arm]]):.4f}; {np.mean([o.secs for o in data[arm]]):.0f} s per run; pregnant at the start {np.mean([o.p_start for o in data[arm]]):.0f} (P1: {np.mean([o.p_start for o in data.P1]) if "p_start" in data.P1[0] else float("nan"):.0f})', flush=True)
    sc.save(fn, data)
    tbode4.save_versions(outdir, 'v4')
    print(f'Added arms done ({T.toc(output=True):.0f} s)', flush=True)
    summarize()
    return data


#%% Summaries

def by_year(outs, key, years):
    """ [rep, year] of a yearly output: 'pop' (÷ 1950), 'cbr', 'cdr', or 'share_<group>' """
    rows = []
    for o in outs:
        yi = {int(y): i for i, y in enumerate(o.years)}
        if key == 'pop':
            v = o.total/o.total[0]
            rows.append([v[yi[y]] for y in years])
        elif key.startswith('share_'):
            sh = r28.shares5(o.pyr)
            k = list(r28.groups5).index(key[6:])
            rows.append([sh[yi[y], k] for y in years])
        else:
            rt = o.rates.set_index('year')[key]
            rows.append([rt.get(y, np.nan) for y in years])
    return np.array(rows)


def ratio_table(data, key, years):
    """ Mean ± SE of each comparison, per year: P0 ÷ R, P1 ÷ R, P1w ÷ P1, O1 ÷ P0 (SE from the agents' reps) """
    v = {a: by_year(data[a], key, years) for a in arms if a in data}
    out = {}
    def fmt(x):
        return [f'{m:.4f} ± {s:.4f}' for m, s in zip(x.mean(0), x.std(0, ddof=1)/np.sqrt(len(x)) if len(x) > 1 else np.zeros(x.shape[1]))]
    out['P0 ÷ R'] = fmt(v['P0']/v['R'])
    out['P1 ÷ R'] = fmt(v['P1']/v['R'])
    out['P1w ÷ P1'] = fmt(v['P1w']/v['P1'].mean(0))
    out['O1 ÷ P0'] = fmt(v['O1']/v['P0'])
    out['P1w ÷ R'] = fmt(v['P1w']/v['R'])
    out['R2 ÷ R'] = fmt(v['R2']/v['R'])
    out['R4 ÷ R'] = fmt(v['R4']/v['R'])
    out['P1 ÷ R4'] = fmt(v['P1']/v['R4'])
    if 'P2' in v:
        out['R8 ÷ R'] = fmt(v['R8']/v['R'])
        out['P2 ÷ P1'] = fmt(v['P2']/v['P1'].mean(0))
        out['P3 ÷ P2'] = fmt(v['P3']/v['P2'].mean(0))
        out['P2 ÷ R4'] = fmt(v['P2']/v['R4'])
        out['P3 ÷ R4'] = fmt(v['P3']/v['R4'])
        out['P3 ÷ R8'] = fmt(v['P3']/v['R8'])
        out['P3 ÷ R2'] = fmt(v['P3']/v['R2'])
    if 'P3s' in v:
        out['P3s ÷ R4'] = fmt(v['P3s']/v['R4'])
        out['P3s ÷ R8'] = fmt(v['P3s']/v['R8'])
        out['P3s ÷ P3'] = fmt(v['P3s']/v['P3'].mean(0))
    if 'P2b' in v:
        out['R8x ÷ R8'] = fmt(v['R8x']/v['R8'])
        out['R4x ÷ R4'] = fmt(v['R4x']/v['R4'])
        out['P2b ÷ R4x'] = fmt(v['P2b']/v['R4x'])
        out['P2b ÷ R8x'] = fmt(v['P2b']/v['R8x'])
        out['P2b ÷ P2'] = fmt(v['P2b']/v['P2'].mean(0))
    if 'P2s' in v:
        out['P2s ÷ R4'] = fmt(v['P2s']/v['R4'])
        out['P2s ÷ R8'] = fmt(v['P2s']/v['R8'])
    return pd.DataFrame(out, index=years), v


def fertility_table(data):
    """ Realized fertility by age at delivery ÷ WPP's rate (mid-year), mean over reps and over each decade; total fertility; mean age at childbearing """
    rows = []
    mids = np.array(wpp.fert_ages, dtype=float) + 2.5
    for arm in [a for a in ['P0', 'P1', 'P1w', 'P2', 'P3', 'P2s', 'P3s', 'P2b'] if a in data]:
        for d in decades_fert:
            yrs = [y for y in range(d, d + 10) if y < stop - 1]
            real = np.mean([o.asfr_real.loc[yrs].values for o in data[arm]], axis=0) # [year, group]
            inp = np.array([wpp.asfr_at(y + 0.5) for y in yrs])
            row = dict(arm=arm, decade=d)
            row.update({f'{a}': float(real[:, j].sum()/inp[:, j].sum()) if inp[:, j].sum() > 0 else np.nan for j, a in enumerate(wpp.fert_ages)})
            row['TFR'] = float(real.sum()/inp.sum())
            row['mac_real'] = float((real.mean(0)*mids).sum()/real.mean(0).sum())
            row['mac_wpp'] = float((inp.mean(0)*mids).sum()/inp.mean(0).sum())
            rows.append(row)
    return pd.DataFrame(rows)


def summarize():
    data = sc.load(fn)
    years = [y for y in report_years if y < stop]
    yrs_rates = [y for y in report_years if y < stop - 1]
    txt = f'Mean gestation g = {data.g:.4f} years. Agents: {reps} reps of {n_agents:,} at the start; mean ± SE over reps (R and O1 are deterministic)\n\n'
    pop, vpop = ratio_table(data, 'pop', years)
    txt += 'Population ÷ its 1950 value, ratios between arms:\n' + pop.to_string() + '\n\n'
    txt += 'Population ÷ 1950 by arm (mean):\n' + pd.DataFrame({a: vpop[a].mean(0) for a in arms}, index=years).round(4).to_string() + '\n\n'
    sh_rows = []
    for gname in r28.groups5:
        t, v = ratio_table(data, f'share_{gname}', years)
        t.insert(0, 'group', gname)
        sh_rows.append(t)
        if gname == list(r28.groups5)[-1]:
            pass
    shares = pd.concat(sh_rows)
    txt += 'Age shares, ratios between arms:\n' + shares.to_string() + '\n\n'
    cbr, vcbr = ratio_table(data, 'cbr', yrs_rates)
    cdr, vcdr = ratio_table(data, 'cdr', yrs_rates)
    txt += 'Crude birth rate, ratios between arms (calendar years):\n' + cbr.to_string() + '\n\n'
    txt += 'Crude death rate, ratios between arms (calendar years):\n' + cdr.to_string() + '\n\n'
    ft = fertility_table(data)
    ft.to_csv(outdir / 'fertility.csv', index=False)
    with pd.option_context('display.precision', 3):
        txt += 'Realized fertility by mother\'s age at delivery ÷ WPP\'s rate, by decade (mean over reps); TFR ratio; mean age at childbearing (realized, WPP):\n' + ft.to_string(index=False) + '\n\n'
    lrows = []
    for arm in [a for a in ['P0', 'P1', 'P1w', 'P2', 'P3', 'P2s', 'P3s', 'P2b'] if a in data]:
        for o in data[arm]:
            if o.losses is not None:
                l = o.losses.copy(); l['arm'] = arm; l['rep'] = o.rep; lrows.append(l)
    if lrows:
        L = pd.concat(lrows)
        L['decade'] = (L.year//10)*10
        lt = L.groupby(['arm', 'decade'])[['losses', 'births', 'recorded']].sum()
        lt['per_1000_births'] = 1000*lt.losses/lt.births
        txt += 'Pregnancies lost to the mother\'s death (by balance), per 1,000 births, summed over reps:\n' + lt.per_1000_births.unstack(0).round(2).to_string() + '\n'
        txt += f'ss.Pregnancy\'s own miscarriage + stillbirth + abortion + antepartum results over the same runs: {lt.recorded.sum():.0f}, against {lt.losses.sum():.0f} by balance\n\n'
    bal = pd.DataFrame([dict(arm=o.arm, rep=o.rep, **o.balance) for a in arms for o in data[a]])
    bal.to_csv(outdir / 'balance.csv', index=False)
    txt += 'Balance check (population change − (births − deaths), per step; total over the run, against total births):\n' + bal.groupby('arm')[['max_step', 'total', 'births']].agg(['min', 'max']).to_string() + '\n\n'
    txt += 'Run times (s): ' + ', '.join(f'{a} {np.mean([o.secs for o in data[a]]):.0f}' for a in arms) + '\n'
    pd.concat([pop.add_prefix('pop '), cbr.add_prefix('cbr ')], axis=1).to_csv(outdir / 'ratios.csv')
    shares.to_csv(outdir / 'shares.csv')
    with open(outdir / 'summary.log', 'w') as fh:
        fh.write(txt)
    print(txt)
    plot(data, ft)
    plot_ages(data)
    return


def plot(data, ft):
    allyears = list(range(start, stop))
    colors = dict(P0='tab:orange', P1='tab:blue', P1w='tab:green', O1='tab:purple', P2='tab:red', P3='tab:brown', P2s='tab:olive', P3s='tab:pink', P2b='tab:cyan', R4='0.4', R8x='k')
    labels = dict(P0='P0 ÷ R (agents, unshifted)', P1='P1 ÷ R (agents, shifted)', P1w='P1w ÷ R (shifted, weekly)', O1='O1 ÷ R (ODE, shifted older)',
                  P2='P2 ÷ R (shifted, start fixed)', P3='P3 ÷ R (shifted, start fixed, losses offset)', P2s='P2s ÷ R (as P2, half women at the start)', P3s='P3s ÷ R (as P3, half women at the start)', P2b='P2b ÷ R (as P2, births 48.6% girls)', R8x='R8x ÷ R (ODE, eighth-year ages, WPP sex split)', R4='R4 ÷ R (ODE, quarter-year ages)')
    fig, axs = plt.subplots(2, 3, figsize=(17, 8.5))
    for ax, key, title in [(axs[0, 0], 'pop', 'Population, ÷ R'), (axs[0, 1], 'cbr', 'Crude birth rate, ÷ R'), (axs[0, 2], 'cdr', 'Crude death rate, ÷ R')]:
        yrs = allyears if key == 'pop' else allyears[:-1]
        ref = by_year(data.R, key, yrs)[0]
        for a in [a for a in ['P0', 'P1', 'P1w', 'O1', 'P2', 'P3', 'P2s', 'P3s', 'P2b', 'R4', 'R8x'] if a in data]:
            v = by_year(data[a], key, yrs)/ref
            m = v.mean(0)
            ax.plot(yrs, m, color=colors[a], label=labels[a])
            if len(v) > 1:
                se = v.std(0, ddof=1)/np.sqrt(len(v))
                ax.fill_between(yrs, m - 2*se, m + 2*se, color=colors[a], alpha=0.2, lw=0)
        ax.axhline(1, color='k', lw=0.6); ax.axhspan(0.99, 1.01, color='0.9', zorder=0); ax.set_title(title, fontsize=10)
    axs[0, 0].legend(fontsize=8)
    for k, gname in enumerate(['0-14', '50-64', '65+']):
        ax = axs[1, k]
        ref = by_year(data.R, f'share_{gname}', allyears)[0]
        for a in [a for a in ['P0', 'P1', 'P1w', 'O1', 'P2', 'P3', 'P2s', 'P3s', 'P2b', 'R4', 'R8x'] if a in data]:
            v = by_year(data[a], f'share_{gname}', allyears)/ref
            ax.plot(allyears, v.mean(0), color=colors[a])
        ax.axhline(1, color='k', lw=0.6); ax.axhspan(0.99, 1.01, color='0.9', zorder=0); ax.set_title(f'Share aged {gname}, ÷ R', fontsize=10)
    fig.suptitle(f'Agents (mean of {reps} reps, ± 2 SE) and the ODE with fertility shifted older, against the ODE with single years of age (R)')
    fig.tight_layout(); fig.savefig(outdir / 'ratios.png', dpi=110); plt.close(fig)
    fig, axs = plt.subplots(1, len(decades_fert), figsize=(4*len(decades_fert), 3.8), squeeze=False)
    for ax, d in zip(axs[0], decades_fert):
        for a in [a for a in ['P0', 'P1', 'P1w', 'P2', 'P3', 'P2s', 'P3s', 'P2b'] if a in data]:
            row = ft[(ft.arm == a) & (ft.decade == d)].iloc[0]
            ax.plot(wpp.fert_ages, [row[str(x)] for x in wpp.fert_ages], 'o-', ms=3, color=colors[a], label=a)
        ax.axhline(1, color='k', lw=0.6); ax.set_ylim(0, 2); ax.set_xlabel("Mother's age at delivery (5-year group)"); ax.set_title(f'{d}s', fontsize=10)
    axs[0, 0].set_ylabel("Realized fertility ÷ WPP's rate"); axs[0, 0].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(outdir / 'fertility_by_age.png', dpi=110); plt.close(fig)
    return


def plot_ages(data):
    """ ages.png: agents (P2b) against the ODE with WPP's start (R8x, and exact aging from R4x and R8x) and WPP: population, five age shares, pyramids,
    and agents ÷ exact aging by 5-year group and year """
    if 'P2b' not in data:
        return
    yrs = list(range(start, stop))
    def pyr5(o): # Shares by 5-year group and sex, [year, group], for each sex
        tot = o.pyr['female'] + o.pyr['male']
        idx = [list(o.years).index(y) for y in yrs]
        return {s: o.pyr[s][idx]/tot[idx].sum(1, keepdims=True) for s in r28.sexes}
    ag = [pyr5(o) for o in data.P2b]
    q8, q4 = pyr5(data.R8x[0]), pyr5(data.R4x[0])
    ex = {s: 2*q8[s] - q4[s] for s in r28.sexes} # Exact aging, extrapolated
    w = r28.wpp_population(yrs)
    wsh = {s: w.pyr[s]/(w.pyr['female'] + w.pyr['male']).sum(1, keepdims=True) for s in r28.sexes}
    both = lambda d: d['female'] + d['male']
    lows = np.array(r28.pyr_edges, dtype=float)
    g5 = lambda m: r28.regroup(m, lows, [gg[0] for gg in r28.groups5.values()])
    fig = plt.figure(figsize=(22, 9.5))
    gs = fig.add_gridspec(2, 6, height_ratios=[1, 1.15])
    ax = fig.add_subplot(gs[0, 0])
    pop = r29_pop = by_year(data.P2b, 'pop', yrs)
    ax.fill_between(yrs, pop.min(0), pop.max(0), color='tab:blue', alpha=0.25, lw=0)
    ax.plot(yrs, pop.mean(0), color='tab:blue', label='Agents (P2b), mean and range of 10 reps')
    ax.plot(yrs, by_year(data.R8x, 'pop', yrs)[0], 'k--', lw=1.2, label='ODE, eighth-year ages (R8x)')
    pex = 2*by_year(data.R8x, 'pop', yrs)[0] - by_year(data.R4x, 'pop', yrs)[0]
    ax.plot(yrs, pex, color='tab:red', ls=':', lw=1.6, label='ODE, exact aging (extrapolated)')
    ax.plot(yrs, w.total/w.total[0], 'o', color='0.5', ms=2.5, label='WPP (includes migration)')
    ax.set_title('Population ÷ 1950', fontsize=10); ax.legend(fontsize=7)
    for k, gname in enumerate(r28.groups5):
        ax = fig.add_subplot(gs[0, k + 1])
        a = np.array([g5(both(d))[:, k] for d in ag])
        ax.fill_between(yrs, a.min(0), a.max(0), color='tab:blue', alpha=0.25, lw=0)
        ax.plot(yrs, a.mean(0), color='tab:blue')
        ax.plot(yrs, g5(both(q8))[:, k], 'k--', lw=1.2)
        ax.plot(yrs, g5(both(ex))[:, k], color='tab:red', ls=':', lw=1.6)
        ax.plot(yrs, g5(both(wsh))[:, k], 'o', color='0.5', ms=2.5)
        ax.set_title(f'Share aged {gname}', fontsize=10)
    show = [y for y in [1960, 1990, 2020, 2050] if y in yrs]
    keep = lows < 100
    for c, y in enumerate(show):
        ax = fig.add_subplot(gs[1, c])
        i = yrs.index(y)
        mf = np.mean([d['female'][i] for d in ag], 0); mm = np.mean([d['male'][i] for d in ag], 0)
        ax.barh(lows[keep], -mm[keep], height=4.5, color='tab:blue', alpha=0.5, label='Agents (P2b), male / female')
        ax.barh(lows[keep], mf[keep], height=4.5, color='tab:pink', alpha=0.5)
        ax.plot(-ex['male'][i][keep], lows[keep], 'r-', lw=1.2, label='ODE, exact aging')
        ax.plot(ex['female'][i][keep], lows[keep], 'r-', lw=1.2)
        ax.axvline(0, color='k', lw=0.5); ax.set_title(f'{y}: share by 5-year group and sex', fontsize=10); ax.set_ylabel('Age' if c == 0 else '')
        if c == 0:
            ax.legend(fontsize=7, loc='upper right')
    ax = fig.add_subplot(gs[1, len(show):])
    ratio = np.mean([both(d) for d in ag], 0)/both(ex) # [year, group]
    top = lows < 85
    im = ax.pcolormesh(yrs, lows[top] + 2.5, ratio[:, top].T, cmap='RdBu_r', vmin=0.95, vmax=1.05, shading='nearest')
    fig.colorbar(im, ax=ax, label='Agents ÷ ODE (exact aging)')
    ax.set_title('Agents ÷ ODE, share in each 5-year group (to 80–84), by year', fontsize=10); ax.set_ylabel('Age group'); ax.set_xlabel('Year')
    fig.suptitle("Ages with both pregnancy fixes and WPP's start in both models: agents (P2b) against the ODE (R8x, and exact aging from R4x and R8x) and WPP", fontsize=11)
    fig.tight_layout(); fig.savefig(outdir / 'ages.png', dpi=110); plt.close(fig)
    return


if __name__ == '__main__':
    what = sys.argv[1] if len(sys.argv) > 1 else 'summary'
    T = sc.timer()
    if what == 'run':
        run_all()
    elif what == 'more':
        run_more()
    elif what == 'summary':
        summarize()
    else:
        raise ValueError(f'Unknown step "{what}"')
    print(f'{what} finished in {T.toc(output=True):.0f} s')
