"""
Experiment 28, v4 environment: do v4's agents match its ODE with India's WPP fertility and mortality, ages and a growing population?
See README.md. Run from the repo root:

    PYTHONPATH=v4 v4/.venv/bin/python experiments/28_age_demography/run.py inputs   # the bundled rates against tb_LAI_TPT's annual files
    PYTHONPATH=v4 v4/.venv/bin/python experiments/28_age_demography/run.py demog    # stage 1: demography alone (D1–D5)
    PYTHONPATH=v4 v4/.venv/bin/python experiments/28_age_demography/run.py burnin   # stage 2, ODE only (B1–B3, and the main arm's ODE)
    PYTHONPATH=v4 v4/.venv/bin/python experiments/28_age_demography/run.py tb       # stage 2: TB with demography, agents and the ODE
    PYTHONPATH=v4 v4/.venv/bin/python experiments/28_age_demography/run.py single   # stage 2, arm S: the ODE with single years of age (after tb)
    PYTHONPATH=v4 v4/.venv/bin/python experiments/28_age_demography/run.py summary  # summaries and figures from the saved runs

Environment variables: EXP28_NCPUS (cores for the agent runs; default tbode4.ncpus), EXP28_QUICK=1 (a small, short test of every
step, written to quick/ instead of this folder).
"""
import os
os.environ.setdefault('NUMBA_NUM_THREADS', '1') # One thread per process (Starsim's kernels otherwise use every core); before importing starsim
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


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


r27 = load_module('run27', root / 'experiments' / '27_calibrated_treatment' / 'run.py') # Exp 27's sets, starting state and per-step outputs

# Settings (README)
quick = os.environ.get('EXP28_QUICK') == '1'
ncpus = int(os.environ.get('EXP28_NCPUS', tbode4.ncpus))
outdir = here / 'quick' if quick else here
outdir.mkdir(exist_ok=True)
lai_data = root.parent / 'tb_LAI_TPT' / 'data'
edges = [0, 1] + list(range(5, 101, 5)) # WPP's mortality groups: 0, 1–4, 5–9, ..., 95–99, 100+
edges1 = list(range(0, 101)) # D2: single years of age
pyr_edges = list(range(0, 101, 5)) # WPP's population groups, for the pyramids
groups5 = sc.objdict({'0-14': (0, 15), '15-29': (15, 30), '30-49': (30, 50), '50-64': (50, 65), '65+': (65, 999)})
sexes = ['female', 'male']
decades = np.arange(1950, 2101, 10)
start = 1950
dts = sc.objdict(month=ss.months(1), week=ss.weeks(1))
sets_used = [0, 6] # Exp 27's set 0 (exp 24's first IMIS seed) and set 6 (its second)
windows = sc.objdict(pre=(1950, 1995), scaleup=(1995, 2025), late=(2025, 2041)) # 1950–1994, 1995–2024, 2025–2040
n_agents, reps1, reps2 = 100_000, 10, 20
stop1, stop2, stop_b2 = 2051, 2041, 2151
report_years = [1950, 1960, 1970, 1980, 1990, 2000, 2010, 2020, 2030, 2040, 2050]
if quick:
    n_agents, reps1, reps2 = 5_000, 2, 2
    stop1, stop2, stop_b2 = 1961, 1998, 1980
    report_years = [1950, 1955, 1960]
    windows = sc.objdict(pre=(1950, 1995), scaleup=(1995, 1998), late=(1997, 1998))
arms1 = sc.objdict(D1='ode', D2='ode', D3='abm', D4='abm', D5='abm')
states_tb = ['susceptible', 'infection', 'infection_remote', 'cleared', 'recovered', 'non_infectious', 'asymptomatic', 'symptomatic', 'treated',
             'treatment_success', 'treatment_fail', 'relapsing']
outputs = ['prevalence', 'prop_asymptomatic', 'incidence', 'notifications', 'on_treatment', 'mortality', 'mortality_tx', 'relapse', 'recent', 'remote', 'cleared']
outputs += [f'prev_{g}' for g in groups5] + [f'share_{g}' for g in groups5] + ['alive']
stocks = ['prevalence', 'prop_asymptomatic', 'on_treatment', 'recent', 'remote', 'cleared', 'alive'] + [f'prev_{g}' for g in groups5] + [f'share_{g}' for g in groups5]
unscaled = ['prop_asymptomatic', 'alive'] + [f'share_{g}' for g in groups5] # Not per 100,000
per_case = ['incidence', 'notifications', 'on_treatment', 'mortality', 'mortality_tx']
compared = ['prevalence', 'prop_asymptomatic', 'incidence', 'notifications', 'on_treatment', 'mortality', 'mortality_tx', 'recent', 'remote', 'cleared'] + \
           [f'prev_{g}' for g in groups5] + [f'share_{g}' for g in groups5] + ['alive']
burnin_outputs = ['prevalence'] + [f'prev_{g}' for g in groups5] + ['recent', 'remote', 'incidence', 'mortality'] + [f'share_{g}' for g in groups5]
fn = sc.objdict(inputs=outdir / 'inputs.csv', demog=outdir / 'demog.obj', burnin=outdir / 'burnin.obj', ode=outdir / 'ode.obj', runs=outdir / 'runs.obj',
                ode_single=outdir / 'ode_single.obj')


#%% Demography inputs: WPP 2024 as bundled with starsim_data, interpolated linearly by year between its decades

class WPP:
    """ The bundled rates as arrays by decade, and their linear interpolation to any year """
    def __init__(self):
        mx = ss.data.rates('India', 'deaths')
        fr = ss.data.rates('India', 'fertility')
        self.mx = {s: mx[mx.Sex == s.capitalize()].pivot(index='Time', columns='AgeGrpStart', values='mx').loc[decades, edges].values for s in sexes} # [decade, age]
        self.fert_ages = sorted(fr.AgeGrpStart.unique()) # 10, 15, ..., 50
        self.asfr = fr.pivot(index='Time', columns='AgeGrpStart', values='ASFR').loc[decades, self.fert_ages].values/1000 # Per woman per year
        return

    @staticmethod
    def at(table, year):
        return np.array([np.interp(year, decades, table[:, j]) for j in range(table.shape[1])])

    def mx_at(self, year): # [sex, age group of `edges`]
        return np.array([self.at(self.mx[s], year) for s in sexes])

    def asfr_at(self, year):
        return self.at(self.asfr, year)


wpp = WPP()


def rate(values, frozen):
    """ A rate per year over time: a timeseries on WPP's decades (linear between them), or frozen at 1950's value """
    if frozen:
        return ss.peryear(float(values[0]))
    return ss.peryear(ss.timeseries({float(y): float(v) for y, v in zip(decades, values)}, interpolation='linear'))


def death_rate(frozen=False):
    """ ss.death of everyone, by age group (WPP's) and sex """
    rows = [dict(age=a, sex=s, value=rate(wpp.mx[s][:, j], frozen)) for s in sexes for j, a in enumerate(edges)]
    return ss.death(rate=ss.by(pd.DataFrame(rows)))


def ode_demographics(frozen=False):
    """ The ODE's births (each woman's age-specific fertility, by her age group; none for men) and deaths """
    zero = ss.peryear(0.0)
    fert = {a: rate(wpp.asfr[:, j], frozen) for j, a in enumerate(wpp.fert_ages)}
    rows = [dict(age=a, sex='female', value=fert.get(a, zero)) for a in edges] + [dict(age=a, sex='male', value=zero) for a in edges]
    return [ss.birth(rate=ss.by(pd.DataFrame(rows))), death_rate(frozen)]


def annual_asfr():
    """ The bundled fertility rates interpolated to every year, as a table for ss.Pregnancy (per 1,000 women) """
    years = np.arange(1950, 2101)
    return pd.DataFrame([dict(Time=y, AgeGrp=a, ASFR=1000*v) for y in years for a, v in zip(wpp.fert_ages, wpp.asfr_at(y))])


def lai_files():
    return sc.objdict(fert=pd.read_csv(lai_data / 'fertility_rates_India.csv', index_col=0), deaths=pd.read_csv(lai_data / 'death_rates_India.csv', index_col=0),
                      cbr=pd.read_csv(lai_data / 'crude_birth_rates_India.csv', index_col=0), cdr=pd.read_csv(lai_data / 'crude_death_rates_India.csv', index_col=0))


def agent_demographics(kind):
    """ D3 ('matched'): ss.Pregnancy set to match the ODE, and the ODE's deaths; D4 ('builtin'): Starsim's modules with the bundled data, as
    tb_LAI_TPT would call them; D5 ('lai'): the same modules with tb_LAI_TPT's annual files """
    if kind == 'matched':
        return [ss.Pregnancy(fertility_rate=annual_asfr(), min_age=10, max_age=55), death_rate()]
    if kind == 'builtin':
        return [ss.Pregnancy(fertility_rate='India'), ss.Deaths(death_rate='India')]
    if kind == 'lai':
        f = lai_files()
        return [ss.Pregnancy(fertility_rate=f.fert[['Time', 'AgeGrp', 'ASFR']]), ss.Deaths(death_rate=f.deaths[['Time', 'Sex', 'AgeGrpStart', 'mx']])]
    raise ValueError(kind)


#%% Helpers for results by age

def by_age_sex(res, sex=None):
    """ A stratified result split by age group (optionally one sex): values [time, group] and the groups' lower edges """
    x = res.select(sex=sex) if sex is not None else res
    b = x.by('age')
    return np.asarray(b.values, dtype=float), np.array(ss.data.age_edges(b.columns), dtype=float)


def regroup(values, lows, bins):
    """ Sum the columns of `values` (groups with lower edges `lows`) into the groups given by `bins` (lower edges, the last open) """
    idx = np.clip(np.digitize(lows, bins) - 1, 0, len(bins) - 1)
    out = np.zeros((values.shape[0], len(bins)))
    for j, k in enumerate(idx):
        out[:, k] += values[:, j]
    return out


def to_groups5(values, lows):
    return regroup(values, lows, [g[0] for g in groups5.values()])


def year_index(sim, year):
    yv = np.asarray(sim.t.yearvec, dtype=float)
    return int(np.argmin(np.abs(yv - year)))


def life_expectancy(mx_row):
    """ Life expectancy at birth from death rates in `edges`' groups, constant within each, with an open last group """
    widths = np.diff(edges)
    S, e = 1.0, 0.0
    for m, w in zip(mx_row[:-1], widths):
        e += S*(1 - np.exp(-m*w))/m
        S *= np.exp(-m*w)
    return e + S/mx_row[-1]


#%% Input check

def run_inputs():
    """ The bundled rates, interpolated by year, against tb_LAI_TPT's annual files: total fertility, life expectancy at birth, and the largest relative difference by age """
    f = lai_files()
    fert, deaths = f.fert, f.deaths
    rows = []
    for y in range(1950, 2051):
        lf = fert[fert.Time == y].set_index('AgeGrp').ASFR/1000
        lai5 = np.array([lf[(lf.index >= a) & (lf.index < a + 5)].mean() if ((lf.index >= a) & (lf.index < a + 5)).any() else 0.0 for a in wpp.fert_ages])
        b5 = wpp.asfr_at(y)
        both = (lai5 > 0) & (b5 > 0)
        row = dict(year=y, tfr_bundled=5*b5.sum(), tfr_lai=lf.sum(), tfr_bundled_15_49=5*b5[1:8].sum(),
                   asfr_max_rel_diff=float(np.max(np.abs(b5[both]/lai5[both] - 1))), asfr_age_of_max=int(np.array(wpp.fert_ages)[both][np.argmax(np.abs(b5[both]/lai5[both] - 1))]))
        mxb = wpp.mx_at(y)
        for k, s in enumerate(sexes):
            d = deaths[(deaths.Time == y) & (deaths.Sex == s.capitalize())].set_index('AgeGrpStart').mx.loc[edges].values
            row[f'e0_{s}_bundled'] = life_expectancy(mxb[k])
            row[f'e0_{s}_lai'] = life_expectancy(d)
            rel = np.abs(mxb[k]/d - 1)
            row[f'mx_{s}_max_rel_diff'] = float(rel.max())
            row[f'mx_{s}_age_of_max'] = int(np.array(edges)[np.argmax(rel)])
        rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(fn.inputs, index=False)
    with open(outdir / 'inputs.log', 'w') as fh, pd.option_context('display.width', 250, 'display.max_columns', 30, 'display.precision', 3):
        txt = ('Bundled WPP rates (decades, linear between them) against tb_LAI_TPT annual files, every 5 years. TFR from 5-year groups (bundled, 10-54) '
               'and single years (tb_LAI_TPT, 15-49); tfr_bundled_15_49 leaves out 10-14 and 50-54. e0 from the death rates (constant within each group).\n')
        txt += df[df.year % 5 == 0].to_string(index=False) + '\n'
        txt += f'\nLargest relative difference over 1950-2050: ASFR {df.asfr_max_rel_diff.max():.3f}, mx female {df.mx_female_max_rel_diff.max():.3f}, mx male {df.mx_male_max_rel_diff.max():.3f}\n'
        txt += f'Largest |TFR difference| (15-49): {np.max(np.abs(df.tfr_bundled_15_49 - df.tfr_lai)):.3f}; |e0 difference|: female {np.max(np.abs(df.e0_female_bundled - df.e0_female_lai)):.2f}, male {np.max(np.abs(df.e0_male_bundled - df.e0_male_lai)):.2f} years\n'
        fh.write(txt)
        print(txt)
    fig, axs = plt.subplots(1, 3, figsize=(15, 4))
    axs[0].plot(df.year, df.tfr_bundled, label='Bundled, 10–54'); axs[0].plot(df.year, df.tfr_bundled_15_49, label='Bundled, 15–49'); axs[0].plot(df.year, df.tfr_lai, '--', label="tb_LAI_TPT's file, 15–49")
    axs[0].set_title('Total fertility'); axs[0].legend()
    for s in sexes:
        axs[1].plot(df.year, df[f'e0_{s}_bundled'], label=f'{s}, bundled'); axs[1].plot(df.year, df[f'e0_{s}_lai'], '--', label=f"{s}, tb_LAI_TPT's file")
    axs[1].set_title('Life expectancy at birth (years)'); axs[1].legend()
    axs[2].plot(df.year, df.asfr_max_rel_diff, label='Fertility'); axs[2].plot(df.year, df.mx_female_max_rel_diff, label='Mortality, female'); axs[2].plot(df.year, df.mx_male_max_rel_diff, label='Mortality, male')
    axs[2].set_title('Largest relative difference in any age group'); axs[2].legend()
    fig.tight_layout(); fig.savefig(outdir / 'inputs.png', dpi=120); plt.close(fig)
    return df


#%% Stage 1: demography alone

def make_demog_sim(arm, rand_seed=0):
    solver = arms1[arm]
    kw = dict(people=ss.population('India', year=start, n_agents=n_agents if solver == 'abm' else None), total_pop=n_agents if solver == 'abm' else 1.0)
    if solver == 'abm':
        kw['rand_seed'] = rand_seed
    demog = ode_demographics() if solver == 'ode' else agent_demographics(dict(D3='matched', D4='builtin', D5='lai')[arm])
    sim = ss.Sim(diseases=ss.SIS(beta=0, init_prev=0), demographics=demog, start=start, stop=stop1, dt=ss.months(1), solver=solver, verbose=0, **kw)
    sim.stratify(age=edges1 if arm == 'D2' else edges, sex=sexes)
    return sim


def demog_outputs(sim):
    """ Yearly (1 January): the population, by WPP's 5-year groups and sex; and crude birth and death rates per 1,000 over each calendar year """
    r = sim.results
    yv = np.asarray(sim.t.yearvec, dtype=float)
    years = np.arange(start, int(np.floor(yv[-1])) + 1)
    idx = np.array([year_index(sim, y) for y in years])
    pyr = {}
    for s in sexes:
        v, lows = by_age_sex(r.sis.S, s)
        pyr[s] = regroup(v, lows, pyr_edges)[idx]
    alive = np.asarray(r.alive.values, dtype=float)
    births = np.asarray(r.sis.births.values if 'births' in r.sis else r.pregnancy.birth.values, dtype=float)
    deaths = np.asarray(r.died.values, dtype=float)
    rows = []
    for y in years[:-1]:
        m = (yv >= y - 1e-9) & (yv < y + 1 - 1e-9) # Flows indexed by their step's start
        rows.append(dict(year=y, cbr=1000*births[m].sum()/alive[m].mean(), cdr=1000*deaths[m].sum()/alive[m].mean()))
    total = pyr['female'].sum(axis=1) + pyr['male'].sum(axis=1)
    out = sc.objdict(years=years, total=total, pyr=pyr, alive=alive[idx], rates=pd.DataFrame(rows), asfr_real=None)
    if 'pregnancy' in r: # Agents: realized fertility by mother's age group (deliveries in each calendar year ÷ mean women), for the inputs' groups
        try:
            fbins = list(wpp.fert_ages) + [55]
            dv = regroup(*by_age_sex(r.pregnancy.delivery, 'female'), fbins)[:, :-1]
            women = regroup(*by_age_sex(r.sis.S, 'female'), fbins)[:, :-1]
            real = []
            for y in years[:-1]:
                m = (yv >= y - 1e-9) & (yv < y + 1 - 1e-9)
                real.append(dv[m].sum(axis=0)/women[m].mean(axis=0))
            out.asfr_real = pd.DataFrame(real, index=years[:-1], columns=wpp.fert_ages)
        except Exception as E:
            print(f'Realized fertility not available: {E}', flush=True)
    return out


def run_demog_one(arm, rep=0):
    T = sc.timer()
    sim = make_demog_sim(arm, rand_seed=28000 + 100*list(arms1).index(arm) + rep)
    sim.run()
    out = demog_outputs(sim)
    out.update(arm=arm, rep=rep, secs=T.toc(output=True))
    return out


def run_demog():
    T = sc.timer()
    data = sc.objdict()
    for arm in [a for a, s in arms1.items() if s == 'ode']:
        data[arm] = [run_demog_one(arm)]
        print(f'{arm} done ({T.toc(output=True):.0f} s)', flush=True)
    jobs = [dict(arm=a, rep=r) for a, s in arms1.items() if s == 'abm' for r in range(reps1)]
    res = sc.parallelize(run_demog_one, iterkwargs=jobs, ncpus=min(ncpus, len(jobs)))
    for arm in [a for a, s in arms1.items() if s == 'abm']:
        data[arm] = [o for o in res if o.arm == arm]
    sc.save(fn.demog, data)
    tbode4.save_versions(outdir, 'v4')
    print(f'Stage 1 done ({T.toc(output=True):.0f} s; agent runs {np.mean([o.secs for o in res]):.0f} s each on average)', flush=True)
    summarize_demog()
    return data


def wpp_population(years):
    """ WPP's population (interpolated between decades by ss.data.population): totals, and shares by 5-year group and sex """
    out = sc.objdict(total=[], pyr={s: [] for s in sexes})
    for y in years:
        p = ss.data.population('India', year=int(y))
        lows = np.array(ss.data.age_edges(p.index), dtype=float)
        for s in sexes:
            out.pyr[s].append(regroup(p[s].values[None, :], lows, pyr_edges)[0])
        out.total.append(p.female.sum() + p.male.sum())
    out.total = np.array(out.total)
    out.pyr = {s: np.array(v) for s, v in out.pyr.items()}
    return out


def shares5(pyr):
    """ Shares of the population in the five groups, from a pyramid by WPP's 5-year groups and sex [time, group] """
    both = pyr['female'] + pyr['male']
    g = regroup(both, np.array(pyr_edges, dtype=float), [g[0] for g in groups5.values()])
    return g/both.sum(axis=1, keepdims=True)


def cohort_spread(years_after):
    """ The ODE's aging alone: the distribution of model age (groups' midpoints) for people born at time 0, by the exponential
    exits from each of `edges`' groups at 1 ÷ its width (no deaths); mean and SD """
    widths = np.diff(edges).astype(float)
    n = len(edges)
    rates = np.append(1/widths, 0.0)
    mids = np.append(np.array(edges[:-1]) + widths/2, 102.5)
    out = []
    p = np.zeros(n); p[0] = 1.0
    h = 0.01
    t = 0.0
    for target in years_after:
        while t < target - 1e-9:
            flow = p*rates*h
            p = p - flow
            p[1:] += flow[:-1]
            t += h
        mean = (p*mids).sum()
        out.append(dict(years=target, mean_age=mean, sd_age=np.sqrt((p*(mids - mean)**2).sum())))
    return pd.DataFrame(out)


def summarize_demog():
    data = sc.load(fn.demog)
    years = data.D1[0].years
    yrs = [y for y in report_years if y in years]
    w = wpp_population(yrs)
    w_share = shares5(w.pyr)
    rows = []
    for arm, runs in data.items():
        for j, y in enumerate(yrs):
            i = list(years).index(y)
            tot = np.array([o.total[i]/o.total[0] for o in runs])
            sh = np.array([shares5({s: o.pyr[s][[i]] for s in sexes})[0] for o in runs])
            row = dict(arm=arm, year=y, pop_ratio=tot.mean(), pop_se=tot.std(ddof=1)/np.sqrt(len(tot)) if len(tot) > 1 else 0.0, wpp_pop_ratio=w.total[j]/w.total[0])
            for k, g in enumerate(groups5):
                row[f'share_{g}'] = sh[:, k].mean()
                row[f'wpp_share_{g}'] = w_share[j, k]
            rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(outdir / 'demog.csv', index=False)
    piv = df.pivot(index='year', columns='arm', values='pop_ratio')
    comp = pd.DataFrame({'D1 vs WPP': piv.D1/df[df.arm == 'D1'].set_index('year').wpp_pop_ratio, 'D2/D1': piv.D2/piv.D1, 'D3/D1': piv.D3/piv.D1, 'D4/D3': piv.D4/piv.D3, 'D5/D4': piv.D5/piv.D4})
    se = df.pivot(index='year', columns='arm', values='pop_se')
    shp = {g: df.pivot(index='year', columns='arm', values=f'share_{g}') for g in groups5}
    with open(outdir / 'demog.log', 'w') as fh, pd.option_context('display.width', 250, 'display.max_columns', 40, 'display.precision', 4):
        txt = f'Stage 1: population ÷ 1950 (mean over reps; agents {reps1} reps of {n_agents:,} at the start), and its ratios between arms\n'
        txt += piv.join(se.add_suffix(' SE')).to_string() + '\n\n' + comp.to_string() + '\n\n'
        txt += 'Age shares (mean over reps), by arm, with WPP:\n'
        for g in groups5:
            t = shp[g].copy(); t['WPP'] = df[df.arm == 'D1'].set_index('year')[f'wpp_share_{g}']
            txt += f'{g}:\n' + t.to_string() + '\n'
        rates = {arm: pd.concat([o.rates.set_index('year') for o in runs]).groupby(level=0).mean() for arm, runs in data.items()}
        f = lai_files()
        cbr_w = f.cbr.set_index(f.cbr.columns[0]) if 'Time' not in f.cbr else f.cbr.set_index('Time')
        txt += '\nCrude birth and death rates per 1,000, every 10 years (calendar-year totals ÷ mean population):\n'
        tab = pd.DataFrame({f'{arm} {k}': rates[arm][k] for arm in rates for k in ['cbr', 'cdr']})
        txt += tab.loc[[y for y in yrs if y in tab.index]].to_string() + '\n'
        txt += '\nAgents\' realized fertility (deliveries ÷ women, by mother\'s age group) ÷ the bundled rates interpolated to mid-year, mean over reps:\n'
        for arm, runs in data.items():
            if runs[0].get('asfr_real') is None:
                continue
            real = sum(o.asfr_real for o in runs)/len(runs)
            ry = [y for y in yrs if y in real.index]
            inp = pd.DataFrame([wpp.asfr_at(y + 0.5) for y in ry], index=ry, columns=wpp.fert_ages)
            ratio = real.loc[ry]/inp
            ratio['TFR'] = real.loc[ry].sum(axis=1)/inp.sum(axis=1)
            txt += f'{arm}:\n' + ratio.to_string() + '\n'
        txt += '\nThe ODE\'s aging alone (exponential exits at 1 ÷ width, no deaths): model age of people born at time 0\n' + cohort_spread([10, 25, 50, 75]).to_string(index=False) + '\n'
        txt += '\nWPP populations are as given by ss.data.population (decadal, interpolated); the ODE and agents leave out migration.\n'
        fh.write(txt)
        print(txt)
    plot_demog(data, years)
    return df


def lai_crude():
    """ WPP's crude birth and death rates by year, from tb_LAI_TPT's files """
    f = lai_files()
    pick = lambda d, col: d.set_index([c for c in d.columns if c.lower() in ['time', 'year']][0])[col]
    vcol = lambda d: [c for c in d.columns if c.lower() not in ['time', 'year', 'location', 'iso3_code', 'loc_id']][-1]
    return pick(f.cbr, vcol(f.cbr)), pick(f.cdr, vcol(f.cdr))


def plot_demog(data, years):
    colors = dict(D1='k', D2='tab:gray', D3='tab:blue', D4='tab:orange', D5='tab:green')
    labels = dict(D1='D1 ODE', D2='D2 ODE, single years', D3='D3 agents, matched', D4='D4 agents, built-in', D5="D5 agents, tb_LAI_TPT's files")
    wy = [y for y in years if y % 10 == 0]
    w = wpp_population(wy)
    w_share = shares5(w.pyr)
    fig, axs = plt.subplots(2, 3, figsize=(16, 8))
    ax = axs.flat[0]
    for arm, runs in data.items():
        tot = np.array([o.total/o.total[0] for o in runs])
        ax.plot(years, tot.mean(axis=0), color=colors[arm], label=labels[arm], lw=1.5 if arm in ['D1', 'D3'] else 1)
    ax.plot(wy, w.total/w.total[0], 'ro', label='WPP')
    ax.set_title('Population ÷ 1950'); ax.legend(fontsize=8)
    for k, g in enumerate(groups5):
        ax = axs.flat[k + 1]
        for arm, runs in data.items():
            sh = np.array([shares5(o.pyr)[:, k] for o in runs]).mean(axis=0)
            ax.plot(years, sh, color=colors[arm], lw=1.5 if arm in ['D1', 'D3'] else 1)
        ax.plot(wy, w_share[:, k], 'ro')
        ax.set_title(f'Share aged {g}')
    fig.tight_layout(); fig.savefig(outdir / 'demog_population.png', dpi=120); plt.close(fig)

    pyr_years = [y for y in [1950, 1980, 2000, 2020, 2050] if y in years] or [years[0], years[-1]]
    w = wpp_population(pyr_years)
    fig, axs = plt.subplots(1, len(pyr_years), figsize=(4*len(pyr_years), 5), sharey=True)
    axs = np.atleast_1d(axs)
    mids = np.array(pyr_edges) + 2.5
    for j, y in enumerate(pyr_years):
        ax = axs[j]
        i = list(years).index(y)
        for arm, runs in data.items():
            tot = lambda o: o.pyr['female'][i].sum() + o.pyr['male'][i].sum()
            f = np.mean([o.pyr['female'][i]/tot(o) for o in runs], axis=0)
            m = np.mean([o.pyr['male'][i]/tot(o) for o in runs], axis=0)
            ax.plot(-m, mids, color=colors[arm], lw=1.2); ax.plot(f, mids, color=colors[arm], lw=1.2, label=labels[arm] if j == 0 else None)
        wt = w.pyr['female'][j].sum() + w.pyr['male'][j].sum()
        ax.plot(-w.pyr['male'][j]/wt, mids, 'ro', ms=3); ax.plot(w.pyr['female'][j]/wt, mids, 'ro', ms=3, label='WPP' if j == 0 else None)
        ax.axvline(0, color='k', lw=0.5); ax.set_title(f'{y}: men (left), women (right)'); ax.set_xlabel('Share of population')
    axs[0].set_ylabel('Age (5-year groups)'); axs[0].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(outdir / 'demog_pyramids.png', dpi=120); plt.close(fig)

    cbr_w, cdr_w = lai_crude()
    fig, axs = plt.subplots(1, 2, figsize=(12, 4))
    for arm, runs in data.items():
        rt = pd.concat([o.rates.set_index('year') for o in runs]).groupby(level=0).mean()
        axs[0].plot(rt.index, rt.cbr, color=colors[arm], label=labels[arm]); axs[1].plot(rt.index, rt.cdr, color=colors[arm])
    lo, hi = years[0], years[-1]
    axs[0].plot(cbr_w.loc[lo:hi].index, cbr_w.loc[lo:hi].values, 'r--', label='WPP (tb_LAI_TPT file)'); axs[1].plot(cdr_w.loc[lo:hi].index, cdr_w.loc[lo:hi].values, 'r--')
    axs[0].set_title('Crude birth rate per 1,000'); axs[1].set_title('Crude death rate per 1,000'); axs[0].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(outdir / 'demog_rates.png', dpi=120); plt.close(fig)
    return


#%% Stage 2: TB with demography

def tb_model(s, dt, scenario):
    """ Exp 27's TBTreat for set s, with mu = 0 (demography from outside). B1 and B2: no treatment, and sym_dead held at its 1995 value """
    init, non0, share_flagged = r27.starting_state(s)
    t = s.trend
    nat = {k: v for k, v in s.nat.items() if k not in ['sym_dead', 'mu']}
    if scenario == 'main':
        tx = tbode4.lai_tpt_treatment(row=dict(tx_covg_2010=t.tx_covg_2010, tx_covg_2024=t.tx_covg_2024, efficacy_1995=t.efficacy_1995,
                                               efficacy_2024=t.efficacy_2024, p_relapse=t.p_relapse), step=dts[dt].years, end=stop2)
        sym_dead = ss.peryear(ss.timeseries({1995.0: t.sym_dead_1995, 2024.0: t.sym_dead_2024}, interpolation='linear'))
        tb = tbode4.TBTreat(init=init, course='fixed', tx_rate=tx.tx_rate, tx_success=tx.tx_success, tx_p_relapse=t.p_relapse, dur_relapse_tx=1.5,
                            sym_dead=sym_dead, mu=0, **nat)
    else:
        tb = tbode4.TBTreat(init=init, course='fixed', sym_dead=float(t.sym_dead_1995), mu=0, **nat)
    return tb, non0, share_flagged


def make_tb_sim(i, dt, solver, scenario='main', rand_seed=0):
    """ scenario: 'main'; 'S' (the main arm's ODE with single years of age, the rates still on WPP's groups); 'B1'–'B3' (README) """
    s = sc.load(r27.fn.sets)[i]
    if scenario == 'B3': # Exp 27's ODE, no ages, births replacing every death at mu = 1/70, from 1950
        r27.start, r27.stop = start, stop2
        sim, non0 = r27.make_sim(s, dt, 'ode')
        return sim, non0, 0.0
    tb, non0, share_flagged = tb_model(s, dt, 'main' if scenario == 'S' else scenario)
    abm = solver == 'abm'
    kw = dict(people=ss.population('India', year=start, n_agents=n_agents if abm else None), total_pop=n_agents if abm else 1.0)
    if abm:
        kw.update(rand_seed=rand_seed, analyzers=[tbode4.Incidence(share_flagged=share_flagged), tbode4.TxChecks()])
        demog = agent_demographics('matched')
    else:
        demog = ode_demographics(frozen=scenario == 'B2')
    sim = ss.Sim(tb, demographics=demog, start=start, stop=stop_b2 if scenario == 'B2' else stop2, dt=dts[dt], solver=solver, verbose=0, **kw)
    sim.stratify(age=edges1 if scenario == 'S' else edges, sex=sexes)
    return sim, non0, share_flagged


def ode_incidence_demog(sim, non0, frozen=False):
    """ tbode4.ode_incidence, with the background death rate of NON0 each step taken as that of all non-infectious people, from their ages and sexes """
    r = sim.results.tb
    tb = sim.diseases.tb
    get = lambda k: np.asarray(r[k].values, dtype=float) if k in r else np.zeros(len(r['symptomatic'].values))
    rate = lambda name: float(getattr(tb, name).rate.raw.value)
    non_rec, non_asy = rate('recovery'), rate('worsening')
    yv = np.asarray(sim.t.yearvec, dtype=float)
    non = np.stack([by_age_sex(r.non_infectious, s)[0] for s in sexes], axis=1) # [time, sex, age group of the sim]
    lows = by_age_sex(r.non_infectious, sexes[0])[1]
    wpp_idx = np.clip(np.digitize(lows, edges) - 1, 0, len(edges) - 1) # Each of the sim's age groups (e.g. single years) in WPP's group holding it
    mxt = np.array([wpp.mx_at(start if frozen else y)[:, wpp_idx] for y in yv]) # [time, sex, age group of the sim]
    mu = np.divide((non*mxt).sum(axis=(1, 2)), non.sum(axis=(1, 2)), out=np.zeros(len(yv)), where=non.sum(axis=(1, 2)) > 0)
    dt = sim.t.dt_year
    inflow = get('latent_non') + get('remote_non')
    n0 = np.zeros(len(inflow)); first_non = np.zeros(len(inflow))
    n0[0] = non0
    for i in range(len(inflow)): # Index i holds the step from t_i to t_(i+1)
        rr = non_rec + non_asy + mu[i]
        decay = np.exp(-rr*dt)
        a = (1 - decay)/rr
        b = (dt - a)/(rr*dt)
        first_non[i] = non_asy*(n0[i]*a + inflow[i]*b)
        if i + 1 < len(inflow):
            n0[i+1] = n0[i]*decay + inflow[i]*(1 - decay)/(rr*dt)
    return get('latent_asy') + get('remote_asy') + first_non + get('relapse')


def tb_per_step(sim, incidence, by_age=True):
    """ Each output per step, per person (per year for flows); stocks at index i are at t_i, flows the step from t_i """
    out = tbode4.tx_outputs(sim, incidence)
    r = sim.results.tb
    alive = np.asarray(sim.results.alive.values, dtype=float)
    v = lambda k: np.asarray(r[k].values, dtype=float)
    out.prop_asymptomatic = v('asymptomatic')/(v('asymptomatic') + v('symptomatic'))
    out.recent, out.remote, out.cleared = v('infection')/alive, v('infection_remote')/alive, v('cleared')/alive
    out.alive = alive/alive[0] # Population ÷ its starting value (agents count people; the ODE starts at 1)
    if by_age:
        pop = sum(to_groups5(*by_age_sex(r[k])) for k in states_tb if k in r)
        act = to_groups5(*by_age_sex(r.asymptomatic)) + to_groups5(*by_age_sex(r.symptomatic))
        for k, g in enumerate(groups5):
            out[f'prev_{g}'] = act[:, k]/pop[:, k]
            out[f'share_{g}'] = pop[:, k]/pop.sum(axis=1)
        out.pop_check = float(np.max(np.abs(pop.sum(axis=1)/alive - 1))) # The states' sum against everyone alive
    return out


def yearly(out, dt_year):
    """ Yearly means: stocks at each step's end, flows over the step, dated by its middle; per 100,000 except shares and population """
    year = np.floor(out.year[1:] - dt_year/2 + 1e-9).astype(int)
    keys = [k for k in outputs if k in out]
    df = pd.DataFrame({k: (out[k][1:] if k in stocks else out[k][:-1])*(1.0 if k in unscaled else 1e5) for k in keys})
    df['year'] = year
    return df.groupby('year').mean()


def seed_of(i, dt, rep):
    return 28500 + 100*(2*sets_used.index(i) + list(dts).index(dt)) + rep


def run_tb_agents(i, dt, rep):
    T = sc.timer()
    sim, _, _ = make_tb_sim(i, dt, 'abm', rand_seed=seed_of(i, dt, rep))
    sim.run()
    inc, chk = sim.analyzers[0], sim.analyzers[1]
    out = tb_per_step(sim, inc)
    return sc.objdict(set=i, dt=dt, rep=rep, seed=seed_of(i, dt, rep), yearly=yearly(out, sim.t.dt_year), checks=dict(chk.c),
                      pop_check=out.pop_check, secs=T.toc(output=True))


def run_tb_ode(i, dt, scenario='main'):
    T = sc.timer()
    sim, non0, _ = make_tb_sim(i, dt, 'ode', scenario)
    sim.run()
    if scenario == 'B3':
        inc = tbode4.ode_incidence(sim, non0=non0).first
        out = tb_per_step(sim, inc, by_age=False)
    else:
        inc = ode_incidence_demog(sim, non0, frozen=scenario == 'B2')
        out = tb_per_step(sim, inc)
    return sc.objdict(set=i, dt=dt, scenario=scenario, yearly=yearly(out, sim.t.dt_year), pop_check=out.get('pop_check'), secs=T.toc(output=True))


def run_burnin():
    """ Stage 2's ODE-only arms (B1–B3) and the main arm's ODE, at monthly steps (the main arm at weekly steps too) """
    T = sc.timer()
    data = sc.objdict()
    for i in sets_used:
        for scen in ['main', 'B1', 'B2', 'B3']:
            data[f'set{i}_{scen}'] = run_tb_ode(i, 'month', scen)
            print(f'set {i} {scen}: {data[f"set{i}_{scen}"].secs:.0f} s (states sum to alive within {data[f"set{i}_{scen}"].pop_check})', flush=True)
        data[f'set{i}_main_week'] = run_tb_ode(i, 'week', 'main')
    sc.save(fn.burnin, data)
    tbode4.save_versions(outdir, 'v4')
    print(f'Burn-in arms done ({T.toc(output=True):.0f} s)', flush=True)
    summarize_burnin()
    return data


def run_tb():
    T = sc.timer()
    ode = [run_tb_ode(i, dt) for i in sets_used for dt in dts]
    sc.save(fn.ode, ode)
    print(f'ODE runs done ({T.toc(output=True):.0f} s)', flush=True)
    partdir = outdir / 'runs_by_group'
    partdir.mkdir(exist_ok=True)
    groups = [(i, dt) for dt in dts for i in sets_used]
    todo = [(i, dt) for i, dt in groups if not (partdir / f'set{i:02d}_{dt}.obj').exists()] # Saved group by group, so a stopped run resumes
    if todo:
        jobs = [dict(i=i, dt=dt, rep=r) for i, dt in todo for r in range(reps2)]
        res = sc.parallelize(run_tb_agents, iterkwargs=jobs, ncpus=min(ncpus, len(jobs)))
        for i, dt in todo:
            sc.save(partdir / f'set{i:02d}_{dt}.obj', [o for o in res if o.set == i and o.dt == dt])
    runs = sum([sc.load(partdir / f'set{i:02d}_{dt}.obj') for i, dt in groups], [])
    sc.save(fn.runs, runs)
    tbode4.save_versions(outdir, 'v4')
    print(f'{len(runs)} agent runs done ({T.toc(output=True):.0f} s; {np.mean([o.secs for o in runs]):.0f} s each on average; '
          f'states sum to alive within {max(o.pop_check for o in runs):.2e})', flush=True)
    summarize_tb()
    return runs


#%% Stage 2 summaries

def window_mean(df, w):
    lo, hi = windows[w]
    return df[(df.index >= lo) & (df.index < hi)].mean()


def by_window(runs, ode):
    rows = []
    for o in ode:
        rs = [r for r in runs if r.set == o.set and r.dt == o.dt]
        for w in windows:
            ow = window_mean(o.yearly, w)
            ratios = []
            for r in rs:
                aw = window_mean(r.yearly, w)
                d = {k: aw[k]/ow[k] for k in compared}
                d.update({f'{k}_per_case': (aw[k]/aw.prevalence)/(ow[k]/ow.prevalence) for k in per_case})
                ratios.append(d)
            ratios = pd.DataFrame(ratios)
            for k in ratios:
                rows.append(dict(set=o.set, dt=o.dt, window=w, output=k, ratio=ratios[k].mean(), se=ratios[k].sem(), reps=len(ratios)))
    return pd.DataFrame(rows)


def summarize_tb():
    runs, ode = sc.load(fn.runs), sc.load(fn.ode)
    bw = by_window(runs, ode)
    bw.to_csv(outdir / 'by_window.csv', index=False)
    chk = pd.DataFrame([dict(set=r.set, dt=r.dt, **r.checks) for r in runs]).groupby(['set', 'dt']).sum()
    chk.to_csv(outdir / 'pitfall_checks.csv')
    with open(outdir / 'tb.log', 'w') as fh, pd.option_context('display.width', 250, 'display.max_columns', 30, 'display.precision', 4):
        main = ['prevalence', 'prop_asymptomatic', 'incidence_per_case', 'notifications_per_case', 'on_treatment_per_case', 'mortality_per_case',
                'recent', 'remote', 'cleared', 'alive'] + [f'prev_{g}' for g in groups5] + [f'share_{g}' for g in groups5]
        t = bw[bw.output.isin(main)].copy()
        t['value'] = t.apply(lambda x: f'{x.ratio:.3f} ± {x.se:.3f}', axis=1)
        tab = t.pivot_table(index=['output'], columns=['set', 'dt', 'window'], values='value', aggfunc='first').reindex(main)
        txt = f'Agents ÷ ODE (1.00 = match), mean ± SE over {reps2} reps, window means of yearly values; per_case = each ÷ prevalence\n'
        txt += tab.to_string() + '\n\nPitfall checks (summed over reps):\n' + chk.to_string() + '\n'
        fh.write(txt)
        print(txt)
    plot_tb(runs, ode, bw)
    return bw


def plot_tb(runs, ode, bw):
    for i in sets_used:
        os_ = {o.dt: o for o in ode if o.set == i}
        panels = [('prevalence', 'Prevalence per 100,000'), ('incidence', 'Incidence per 100,000 per year'), ('mortality', 'TB mortality per 100,000 per year'),
                  ('notifications', 'Notifications per 100,000 per year')]
        for name, items, ncol in [('burden', panels, 4), ('prev_by_age', [(f'prev_{g}', f'Prevalence per 100,000, aged {g}') for g in groups5], 5),
                                  ('population', [(f'share_{g}', f'Share aged {g}') for g in groups5], 5)]:
            fig, axs = plt.subplots(len(dts), ncol, figsize=(4*ncol, 3.5*len(dts)), squeeze=False)
            for a, dt in enumerate(dts):
                rs = [r for r in runs if r.set == i and r.dt == dt]
                for b, (k, label) in enumerate(items):
                    ax = axs[a, b]
                    for r in rs:
                        ax.plot(r.yearly.index, r.yearly[k], color='tab:blue', alpha=0.15, lw=0.6)
                    m = pd.concat([r.yearly[k] for r in rs], axis=1).mean(axis=1)
                    ax.plot(m.index, m.values, color='tab:blue', lw=1.5, label='Agents, mean')
                    ax.plot(os_[dt].yearly.index, os_[dt].yearly[k], 'k--', lw=1.2, label='ODE')
                    ax.set_title(f'{label} ({dt})', fontsize=9)
                    if a == 0 and b == 0:
                        ax.legend(fontsize=8)
            fig.suptitle(f'Set {i}: agents ({len(rs)} reps) and the ODE'); fig.tight_layout(); fig.savefig(outdir / f'tb_{name}_set{i}.png', dpi=110); plt.close(fig)
    show = ['prevalence', 'prop_asymptomatic', 'incidence_per_case', 'notifications_per_case', 'mortality_per_case', 'recent', 'remote', 'cleared', 'alive'] + \
           [f'prev_{g}' for g in groups5]
    fig, axs = plt.subplots(len(windows), 1, figsize=(16, 3.5*len(windows)), squeeze=False)
    for a, w in enumerate(windows):
        ax = axs[a, 0]
        for c, (i, dt) in enumerate([(i, dt) for i in sets_used for dt in dts]):
            d = bw[(bw.set == i) & (bw.dt == dt) & (bw.window == w)].set_index('output').reindex(show)
            x = np.arange(len(show)) + (c - 1.5)*0.15
            ax.errorbar(x, d.ratio, yerr=2*d.se, fmt='o', ms=3, label=f'set {i}, {dt}')
        ax.axhline(1, color='k', lw=0.6); ax.axhspan(0.99, 1.01, color='0.9')
        ax.set_xticks(np.arange(len(show))); ax.set_xticklabels(show, rotation=30, ha='right', fontsize=8)
        ax.set_title(f'Agents ÷ ODE, {w} ({windows[w][0]}–{windows[w][1] - 1}), mean ± 2 SE'); ax.legend(fontsize=7, ncol=4)
    fig.tight_layout(); fig.savefig(outdir / 'agents_vs_ode.png', dpi=110); plt.close(fig)
    return


def burnin_time(series, tol=0.01):
    """ The first year after which the series stays within tol of its final value for good """
    final = series.iloc[-1]
    off = (np.abs(series/final - 1) > tol).values
    if not off.any():
        return int(series.index[0])
    last = np.flatnonzero(off)[-1]
    return int(series.index[min(last + 1, len(series) - 1)])


def summarize_burnin():
    data = sc.load(fn.burnin)
    rows, eff = [], []
    for i in sets_used:
        b2 = data[f'set{i}_B2'].yearly.iloc[:-1] # The last year is partial
        flat = {k: float(b2[k].iloc[-10:].max()/b2[k].iloc[-10:].min() - 1) for k in burnin_outputs}
        for k in burnin_outputs:
            rows.append(dict(set=i, output=k, burnin_year=burnin_time(b2[k]), steady=float(b2[k].iloc[-1]), last_decade_range=flat[k]))
        main, b3 = data[f'set{i}_main'].yearly, data[f'set{i}_B3'].yearly
        for w in windows:
            mw, bw_ = window_mean(main, w), window_mean(b3, w)
            for k in ['prevalence', 'incidence', 'mortality', 'notifications']:
                eff.append(dict(set=i, window=w, output=k, main_over_b3=mw[k]/bw_[k] if bw_[k] else np.nan))
    df = pd.DataFrame(rows)
    df.to_csv(outdir / 'burnin_times.csv', index=False)
    ef = pd.DataFrame(eff)
    with open(outdir / 'burnin.log', 'w') as fh, pd.option_context('display.width', 250, 'display.max_columns', 30, 'display.precision', 4):
        txt = 'B2 (demography frozen at 1950, untreated): first year after which each output stays within 1% of its value in the last full year\n'
        txt += '(last_decade_range: max ÷ min over the last 10 years, minus 1; should be well under 0.01 if the run has settled)\n'
        txt += df.pivot(index='output', columns='set', values=['burnin_year', 'last_decade_range']).reindex(burnin_outputs).to_string() + '\n\n'
        txt += 'B1 ÷ B2 at selected years:\n'
        for i in sets_used:
            b1, b2 = data[f'set{i}_B1'].yearly, data[f'set{i}_B2'].yearly
            yrs = [y for y in [1960, 1970, 1980, 1994, 2010, 2025, 2040] if y in b1.index and y in b2.index]
            txt += f'set {i}:\n' + (b1.loc[yrs, burnin_outputs]/b2.loc[yrs, burnin_outputs]).T.to_string() + '\n'
        txt += '\nMain arm (with demography) ÷ B3 (exp 27\'s ODE, no ages, mu = 1/70), window means:\n' + ef.pivot_table(index='output', columns=['set', 'window'], values='main_over_b3').to_string() + '\n'
        txt += '\nMain arm, ODE weekly ÷ monthly (step effect in the ODE alone), window means of prevalence:\n'
        for i in sets_used:
            wk, mo = data[f'set{i}_main_week'].yearly, data[f'set{i}_main'].yearly
            txt += f'set {i}: ' + ', '.join(f'{w} {window_mean(wk, w).prevalence/window_mean(mo, w).prevalence:.4f}' for w in windows) + '\n'
        fh.write(txt)
        print(txt)
    for i in sets_used:
        b1, b2 = data[f'set{i}_B1'].yearly, data[f'set{i}_B2'].yearly.iloc[:-1]
        items = ['prevalence'] + [f'prev_{g}' for g in groups5] + ['recent', 'remote']
        fig, axs = plt.subplots(2, 4, figsize=(18, 8))
        for ax, k in zip(axs.flat, items):
            ax.plot(b2.index, b2[k], color='tab:red', label='B2: 1950 demography, frozen')
            ax.plot(b1.index, b1[k], color='tab:blue', label='B1: WPP demography')
            ax.axhline(b2[k].iloc[-1], color='tab:red', ls=':', lw=0.8)
            ax.axvline(df[(df.set == i) & (df.output == k)].burnin_year.iloc[0], color='tab:red', ls='--', lw=0.8)
            ax.set_title(f'{k} (per 100,000)', fontsize=9)
        axs.flat[0].legend(fontsize=8)
        fig.suptitle(f'Set {i}, ODE, untreated: burn-in from the age-independent start (dashed: B2 within 1% for good; dotted: B2 steady state)')
        fig.tight_layout(); fig.savefig(outdir / f'burnin_set{i}.png', dpi=110); plt.close(fig)
        main, b3 = data[f'set{i}_main'].yearly, data[f'set{i}_B3'].yearly
        fig, axs = plt.subplots(1, 4, figsize=(18, 4))
        for ax, k in zip(axs, ['prevalence', 'incidence', 'mortality', 'notifications']):
            ax.plot(main.index, main[k], label='With demography'); ax.plot(b3.index, b3[k], '--', label='B3: no ages, mu = 1/70')
            ax.set_title(f'{k} per 100,000')
        axs[0].legend(fontsize=8)
        fig.suptitle(f'Set {i}, ODE with treatment: the effect of demography'); fig.tight_layout(); fig.savefig(outdir / f'demog_effect_set{i}.png', dpi=110); plt.close(fig)
    return df


#%% Stage 2, arm S: the main arm's ODE with single years of age

def run_single():
    """ Arm S at both sets and steps, compared with the agents and with the main arm's ODE (WPP's 22 groups); needs `tb` first

    Starsim keeps a module's results by group at every step only while they take at most `Engine.group_bytes` (100 MB) in a run, and
    otherwise every few steps; S at weekly steps takes more, so the limit is raised here to keep every step, as for the other arms.
    """
    import resource
    T = sc.timer()
    ss.solvers.Engine.group_bytes = 2e9
    single = []
    for i in sets_used:
        for dt in dts:
            single.append(run_tb_ode(i, dt, 'S'))
            rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1e9 if sys.platform == 'darwin' else 1e6) # Bytes on macOS, kB on Linux
            single[-1].peak_gb = rss
            print(f'set {i} {dt}: {single[-1].secs:.1f} s (states sum to alive within {single[-1].pop_check:.2e}; peak memory so far {rss:.2f} GB)', flush=True)
    sc.save(fn.ode_single, single)
    tbode4.save_versions(outdir, 'v4')
    print(f'Arm S done ({T.toc(output=True):.0f} s)', flush=True)
    summarize_single()
    return single


def summarize_single():
    runs, ode, single = sc.load(fn.runs), sc.load(fn.ode), sc.load(fn.ode_single)
    bw = by_window(runs, single) # Agents ÷ S
    bw.to_csv(outdir / 'by_window_single.csv', index=False)
    rows = [] # S ÷ the main arm's ODE
    for s_ in single:
        o = [o for o in ode if o.set == s_.set and o.dt == s_.dt][0]
        for w in windows:
            sw, ow = window_mean(s_.yearly, w), window_mean(o.yearly, w)
            d = {k: sw[k]/ow[k] for k in compared}
            d.update({f'{k}_per_case': (sw[k]/sw.prevalence)/(ow[k]/ow.prevalence) for k in per_case})
            rows += [dict(set=s_.set, dt=s_.dt, window=w, output=k, ratio=v) for k, v in d.items()]
    sv = pd.DataFrame(rows)
    sv.to_csv(outdir / 'single_vs_groups.csv', index=False)
    show = ['prevalence', 'prop_asymptomatic', 'incidence_per_case', 'notifications_per_case', 'mortality_per_case', 'recent', 'remote', 'cleared', 'alive'] + \
           [f'prev_{g}' for g in groups5] + [f'share_{g}' for g in groups5]
    with open(outdir / 'single.log', 'w') as fh, pd.option_context('display.width', 250, 'display.max_columns', 30, 'display.precision', 4):
        txt = f'Arm S run times: ' + ', '.join(f'set {s_.set} {s_.dt} {s_.secs:.1f} s' for s_ in single)
        txt += ' (main arm, 22 groups: ' + ', '.join(f'set {o.set} {o.dt} {o.secs:.1f} s' for o in ode) + '; the first includes compiling)\n\n'
        t = bw[bw.output.isin(show)].copy()
        t['value'] = t.apply(lambda x: f'{x.ratio:.3f} ± {x.se:.3f}', axis=1)
        txt += f'Agents ÷ S (the ODE with single years of age; 1.00 = match), mean ± SE over {reps2} reps, window means of yearly values\n'
        txt += t.pivot_table(index='output', columns=['set', 'dt', 'window'], values='value', aggfunc='first').reindex(show).to_string() + '\n\n'
        txt += 'S ÷ the main arm\'s ODE (WPP\'s 22 groups), window means\n'
        txt += sv[sv.output.isin(show)].pivot_table(index='output', columns=['set', 'dt', 'window'], values='ratio').reindex(show).to_string() + '\n'
        fh.write(txt)
        print(txt)
    # Prevalence by age and age shares over time: agents' mean, the main arm's ODE and S
    items = [(f'prev_{g}', f'Prevalence per 100,000, aged {g}') for g in groups5] + [(f'share_{g}', f'Share aged {g}') for g in groups5]
    for i in sets_used:
        fig, axs = plt.subplots(2*len(dts), 5, figsize=(20, 3.3*2*len(dts)), squeeze=False)
        for a, dt in enumerate(dts):
            rs = [r for r in runs if r.set == i and r.dt == dt]
            o = [o for o in ode if o.set == i and o.dt == dt][0]
            s_ = [s_ for s_ in single if s_.set == i and s_.dt == dt][0]
            for b, (k, label) in enumerate(items):
                ax = axs[2*a + b//5, b % 5]
                m = pd.concat([r.yearly[k] for r in rs], axis=1).mean(axis=1)
                ax.plot(m.index, m.values, color='tab:blue', lw=1.5, label=f'Agents, mean of {len(rs)}')
                ax.plot(o.yearly.index, o.yearly[k], 'k--', lw=1.2, label="ODE, WPP's 22 groups")
                ax.plot(s_.yearly.index, s_.yearly[k], color='tab:red', ls=':', lw=1.6, label='ODE, single years (S)')
                ax.set_title(f'{label} ({dt})', fontsize=9)
                if a == 0 and b == 0:
                    ax.legend(fontsize=8)
        fig.suptitle(f'Set {i}: agents against the ODE with 22 age groups and with single years'); fig.tight_layout()
        fig.savefig(outdir / f'single_by_age_set{i}.png', dpi=110); plt.close(fig)
    # Agents ÷ ODE by window, for the 22 groups and for S
    bw22 = by_window(runs, ode)
    fig, axs = plt.subplots(len(windows), 1, figsize=(18, 3.6*len(windows)), squeeze=False)
    for a, w in enumerate(windows):
        ax = axs[a, 0]
        for c, (i, dt) in enumerate([(i, dt) for i in sets_used for dt in dts]):
            for e, (d_all, mk, lab) in enumerate([(bw22, 'o', '22 groups'), (bw, 's', 'single years')]):
                d = d_all[(d_all.set == i) & (d_all.dt == dt) & (d_all.window == w)].set_index('output').reindex(show)
                x = np.arange(len(show)) + (2*c + e - 3.5)*0.09
                ax.errorbar(x, d.ratio, yerr=2*d.se, fmt=mk, ms=3, color=f'C{c}', alpha=0.45 if e == 0 else 1.0, label=f'set {i}, {dt}, {lab}')
        ax.axhline(1, color='k', lw=0.6); ax.axhspan(0.99, 1.01, color='0.9')
        ax.set_xticks(np.arange(len(show))); ax.set_xticklabels(show, rotation=30, ha='right', fontsize=8)
        ax.set_title(f'Agents ÷ ODE, {w} ({windows[w][0]}–{windows[w][1] - 1}), mean ± 2 SE: faint circles, 22 groups; squares, single years (S)')
        ax.legend(fontsize=6, ncol=4)
    fig.tight_layout(); fig.savefig(outdir / 'agents_vs_ode_single.png', dpi=110); plt.close(fig)
    return bw, sv


if __name__ == '__main__':
    what = sys.argv[1] if len(sys.argv) > 1 else 'summary'
    T = sc.timer()
    if what == 'inputs':
        run_inputs()
    elif what == 'demog':
        run_demog()
    elif what == 'burnin':
        run_burnin()
    elif what == 'tb':
        run_tb()
    elif what == 'single':
        run_single()
    elif what == 'summary':
        for f, func in [(fn.inputs, None), (fn.demog, summarize_demog), (fn.burnin, summarize_burnin), (fn.runs, summarize_tb), (fn.ode_single, summarize_single)]:
            if func is not None and f.exists():
                func()
    else:
        raise ValueError(f'Unknown step "{what}"')
    print(f'{what} finished in {T.toc(output=True):.0f} s')
