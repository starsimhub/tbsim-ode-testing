"""
Experiment 34, v4 environment: do v4's agents match its ODE with ages and demography, across exp 32's calibrated posterior? See README.md.
Run from the repo root, in order:

    PYTHONPATH=v4 v4/.venv/bin/python experiments/34_calibrated_demography_agents/run.py sets     # draw the 12 sets and their starting tables
    PYTHONPATH=v4 v4/.venv/bin/python experiments/34_calibrated_demography_agents/run.py check    # the mapping check: v4's ODE against fastode_age
    PYTHONPATH=v4 v4/.venv/bin/python experiments/34_calibrated_demography_agents/run.py ode      # v4's ODE at both widths and steps
    PYTHONPATH=v4 v4/.venv/bin/python experiments/34_calibrated_demography_agents/run.py agents   # the agents (saved set by set, so a stopped run resumes)
    PYTHONPATH=v4 v4/.venv/bin/python experiments/34_calibrated_demography_agents/run.py summary  # tables and figures from the saved runs

Environment variables: EXP34_NCPUS (worker processes; default tbode4.ncpus), EXP34_QUICK=1 (a small test of every stage: 1 set, 2 reps of 10,000 agents,
written to quick/).

Both solvers run tbode4.TBTreat with exp 32's parameter mapping (exp 24's calib.to_pars) and schedules, India's WPP demography with Starsim's defaults
(exp 33), from each set's settled start by single year of age and sex (fastode_age.settle(years=100) and start_table(), as exp 32's calibration starts).
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
from tbode4 import fastode_age as fa # noqa: E402

# Settings (README)
quick = os.environ.get('EXP34_QUICK') == '1'
ncpus = int(os.environ.get('EXP34_NCPUS', tbode4.ncpus))
outdir = here / 'quick' if quick else here
outdir.mkdir(exist_ok=True)
exp32 = root / 'experiments' / '32_demography_calibration'
seeds32 = ['real', 'real_seed2'] # Exp 32's two IMIS seeds: 6 sets from each
sets_per_seed, draw_seed = 6, 34
start, stop = 1950, 2041
h = 1/12 # fastode_age's step and age groups (monthly), as exp 32
settle_years = 100
n_agents, reps = 100_000, 20
dts = sc.objdict(month=ss.months(1), week=ss.weeks(1))
widths = sc.objdict(quarter=0.25, eighth=0.125) # The ODE's age groups; exact aging extrapolated as 2 × eighth − quarter
agent_edges = tbode4.age_edges(1.0) # Agents' groups (results only): no coarser than the starting table's single years
table_lows = np.arange(101, dtype=float) # The starting table's single years of age, 100+ open
tol = 1e-3 # Added to agents' ages when counting them by age group, as exps 30 and 33
groups5 = fa.groups5
windows = sc.objdict(pre=(1950, 1995), early=(1995, 2010), late=(2010, 2025), const=(2025, 2041)) # 1950–1994, 1995–2009, 2010–2024, 2025–2040
if quick:
    sets_per_seed, n_agents, reps = 1, 10_000, 2
fn = sc.objdict(sets=outdir / 'sets.obj', ode=outdir / 'ode.obj', runs=outdir / 'runs.obj')
ode_parts, agent_parts = outdir / 'ode_parts', outdir / 'runs_by_group'
wpp = tbode4.WPP()
ss.solvers.Engine.group_bytes = 4e9 # The ODE's fine age groups keep their results at every step (Starsim keeps them only every few steps above 100 MB per run)
tbode4.use_drawn_start() # Agents' starting counts by group drawn, not rounded (proposal 07); every agent sim in this process and its workers

# Outputs
stocks = fa.stock_keys
flows = fa.flow_keys
outputs = ['prevalence', 'prev_15+', 'prop_asymptomatic', 'prop_asymptomatic_15+', 'incidence', 'notifications', 'on_treatment', 'mortality', 'mortality_tx',
           'relapse', 'recent', 'remote', 'cleared', 'alive'] + [f'prev_{g}' for g in groups5] + [f'share_{g}' for g in groups5]
per_case = ['incidence', 'notifications', 'on_treatment', 'mortality', 'mortality_tx'] # Also compared per prevalent case (all ages)
run_targets = ['prevalence', 'prop_asymptomatic', 'prev_to_inc', 'mort_to_inc', 'prev_to_notif', 'inc_ratio', 'mort_ratio', 'prev_2040_ratio']
target_extra = ['prevalence_all_2020', 'prevalence_15plus_1995']


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_cb32 = None


def cb32():
    """ Exp 32's calib.py (its parameter mapping, start and targets; exp 24's underneath), loaded on first use. Only the main process uses it """
    global _cb32
    if _cb32 is None:
        _cb32 = load_module('calib32', exp32 / 'calib.py')
    return _cb32


#%% Sets

def draw_sets():
    """ 6 distinct draws from each of exp 32's IMIS seeds (weighted resampling without replacement), each with its settled start and fastode_age run """
    T = sc.timer()
    cb = cb32()
    assert fa.Demography.birth_female_default == 'young', 'Births at the 1950 share at ages 0–4, as v4\'s ODE at the current pins'
    rng = np.random.default_rng(draw_seed)
    sets, rows = [], []
    for seed32 in seeds32:
        res = sc.load(exp32 / f'{seed32}.obj')
        idx = rng.choice(len(res.weights), size=sets_per_seed, replace=False, p=res.weights)
        for i in idx:
            theta = np.asarray(res.theta[i], dtype=float)
            nat, trend = cb.to_pars(theta)
            init, flagged = cb.untreated_equilibrium(nat)
            p = cb.model_pars(nat, trend)
            st = fa.settle(p, h, init, flagged, years=settle_years, start=start)
            table, share_flagged = fa.start_table(st)
            out = fa.run_cohort(p, h, start, stop, table=table, share_flagged=share_flagged)
            out.trend = trend
            s = sc.objdict(set=len(sets), seed32=seed32, post_index=int(i), weight=float(res.weights[i]), theta=dict(zip(cb.names, map(float, theta))),
                           nat=nat, trend=trend, table=table, share_flagged=share_flagged, drift=st.drift, fast_yearly=fa.yearly(out, h),
                           fast_targets=cb.outputs(out))
            sets.append(s)
            rows.append(dict(set=s.set, seed32=seed32, **s.theta, prev15_1995=1e5*s.fast_targets['prevalence_15plus_1995'], drift=st.drift))
            print(f'Set {s.set} ({seed32}, draw {i}): drift over the frozen start\'s last 20 years {st.drift:.1e}; share flagged {share_flagged:.4f}; '
                  f'15+ prevalence 1995 {rows[-1]["prev15_1995"]:.0f} per 100,000 ({T.toc(output=True):.0f} s)', flush=True)
    sc.save(fn.sets, sets)
    df = pd.DataFrame(rows)
    lines = [f'Exp 34: {len(sets)} sets, {sets_per_seed} from each of exp 32\'s IMIS seeds (weighted resampling without replacement, seed {draw_seed})']
    with pd.option_context('display.width', 300, 'display.max_columns', 40):
        lines += ['\nDrawn sets:', df.round(4).to_string(index=False)]
        for seed32 in seeds32:
            post = pd.read_csv(exp32 / f'{seed32}_posterior.csv').set_index('parameter')
            d = df[df.seed32 == seed32][cb.names]
            inside = ((d >= post.q025.reindex(cb.names)) & (d <= post.q975.reindex(cb.names))).sum()
            lines += [f'\n{seed32}: drawn range against the posterior\'s 95% interval',
                      pd.DataFrame({'min drawn': d.min(), 'max drawn': d.max(), 'post 2.5%': post.q025, 'post 97.5%': post.q975,
                                    'drawn inside': inside}).reindex(cb.names).round(4).to_string()]
        ft = pd.DataFrame([dict(set=s.set, seed32=s.seed32, **{k: s.fast_targets[k] for k in run_targets}) for s in sets])
        lines += ['\nfastode_age (the calibration\'s model) target values by set:', ft.round(4).to_string(index=False)]
    text = '\n'.join(lines)
    sc.savetext(outdir / 'sets.log', text)
    print(text, flush=True)
    tbode4.save_versions(outdir, 'v4')
    return sets


#%% The model in v4

def annual_asfr():
    """ The bundled fertility rates interpolated to every year, as a table for ss.Pregnancy (per 1,000 women), as exps 28 and 33 """
    years = np.arange(1950, 2101)
    return pd.DataFrame([dict(Time=y, AgeGrp=a, ASFR=1000*v) for y in years for a, v in zip(wpp.fert_ages, wpp.asfr_at(y))])


def tb_model(s, dt):
    """ TBTreat for set s (exp 30's tb_model): mu = 0 (WPP's deaths instead), exp 24's schedules, from the set's table by single year of age and sex """
    t = s.trend
    nat = {k: v for k, v in s.nat.items() if k not in ['sym_dead', 'mu']}
    tx = tbode4.lai_tpt_treatment(row=dict(tx_covg_2010=t.tx_covg_2010, tx_covg_2024=t.tx_covg_2024, efficacy_1995=t.efficacy_1995,
                                           efficacy_2024=t.efficacy_2024, p_relapse=t.p_relapse), step=dts[dt].years, end=stop)
    sym_dead = ss.peryear(ss.timeseries({1995.0: t.sym_dead_1995, 2024.0: t.sym_dead_2024}, interpolation='linear'))
    init = tbode4.init_table(s.table, table_lows)
    return tbode4.TBTreat(init=init, course='fixed', tx_rate=tx.tx_rate, tx_success=tx.tx_success, tx_p_relapse=t.p_relapse, dur_relapse_tx=1.5,
                          sym_dead=sym_dead, mu=0, **nat)


class AgeCounts(ss.Analyzer):
    """ Agents alive, with active TB, and with asymptomatic TB in each of exp 28's five age groups, counted with `tol` added to ages (exp 30's,
    with asymptomatic TB added for the 15+ share asymptomatic). Index i is the state at t_i, as Starsim's results of states """
    def init_post(self):
        super().init_post()
        self.lows = np.array([g[0] for g in groups5.values()], dtype=float)
        self.alive = np.zeros((self.sim.t.npts, len(self.lows)))
        self.active = np.zeros_like(self.alive)
        self.asym = np.zeros_like(self.alive)
        self._count(0)
        return

    def _count(self, i):
        ppl, tb = self.sim.people, self.sim.diseases.tb
        if i >= len(self.alive):
            return
        g = lambda u: np.clip(np.digitize(np.asarray(ppl.age.raw[u], dtype=float) + tol, self.lows) - 1, 0, len(self.lows) - 1)
        def living(u): # Leave out those dying in this step: analyzers run before People.step_die() removes them (exps 30 and 33 counted them)
            u = np.asarray(u).astype(int)
            return u[~((ppl.ti_dead.raw[u] <= self.ti) | (ppl.ti_removed.raw[u] <= self.ti))] if i > 0 else u
        n = len(self.lows)
        asy, sym = living(tb.asymptomatic.uids), living(tb.symptomatic.uids)
        self.alive[i] = np.bincount(g(living(ppl.auids)), minlength=n)
        self.asym[i] = np.bincount(g(asy), minlength=n)
        self.active[i] = self.asym[i] + np.bincount(g(sym), minlength=n)
        return

    def step(self): # After this step's transitions: the state at t_(ti+1)
        self._count(self.ti + 1)
        return


def make_sim(s, dt, solver, width='quarter', rand_seed=0):
    """ One sim, 1950–2041. Agents: ss.Pregnancy on WPP's table as given, Starsim's defaults (exp 33). The ODE: births at Starsim's default female
    share, split by sex as its population is at each age """
    abm = solver == 'abm'
    kw = dict(people=ss.population('India', year=start, n_agents=n_agents if abm else None), total_pop=n_agents if abm else 1.0)
    if abm:
        kw.update(rand_seed=rand_seed, analyzers=[tbode4.Incidence(share_flagged=s.share_flagged), tbode4.TxChecks(), AgeCounts()])
        demog = [ss.Pregnancy(fertility_rate=annual_asfr(), rate_per=1000, min_age=10, max_age=55, name='pregnancy'), tbode4.death_rate(wpp)]
    else:
        demog = tbode4.ode_demographics(wpp)
    sim = ss.Sim(tb_model(s, dt), demographics=demog, start=start, stop=stop, dt=dts[dt], solver=solver, verbose=0, **kw)
    sim.stratify(age=agent_edges if abm else tbode4.age_edges(widths[width]), sex=tbode4.sexes)
    return sim


#%% Outputs

def by_age_sex(res, sex=None):
    """ A stratified result by age group (one sex, or both summed): values [time, group] and the groups' lower edges """
    x = res.select(sex=sex) if sex is not None else res
    b = x.by('age')
    return np.asarray(b.values, dtype=float), np.array(ss.data.age_edges(b.columns), dtype=float)


def ode_by_age(sim):
    """ The ODE's people alive, with active TB and asymptomatic TB, in the five groups [time, group], and at 15+ [time, (alive, active, asymptomatic)] """
    r = sim.results.tb
    pop, lows = sum(by_age_sex(r[k])[0] for k in tbode4.tb_states if k in r), by_age_sex(r.susceptible)[1]
    asy = by_age_sex(r.asymptomatic)[0]
    act = asy + by_age_sex(r.symptomatic)[0]
    g = np.clip(np.digitize(lows + 1e-9, [lo for lo, _ in groups5.values()]) - 1, 0, len(groups5) - 1)
    a5 = np.stack([pop[:, g == k].sum(axis=1) for k in range(len(groups5))], axis=1)
    p5 = np.stack([act[:, g == k].sum(axis=1) for k in range(len(groups5))], axis=1)
    adult = lows >= 15 - 1e-9
    return a5, p5, np.stack([pop[:, adult].sum(axis=1), act[:, adult].sum(axis=1), asy[:, adult].sum(axis=1)], axis=1)


def ode_incidence(sim, share_flagged):
    """
    The ODE's first entries to asymptomatic TB per step (plus relapses), as exp 30's: NON0 (non-infectious, not yet asymptomatic this episode) carried
    forward with the exit rate non_rec + non_asy + the background death rate of all non-infectious people (from their ages and sexes, WPP's rates
    that year). NON0 at the start is (1 − share_flagged) × the starting non-infectious TB
    """
    r, tb = sim.results.tb, sim.diseases.tb
    get = lambda k: np.asarray(r[k].values, dtype=float) if k in r else np.zeros(len(r['symptomatic'].values))
    rate = lambda name: float(getattr(tb, name).rate.raw.value)
    non_rec, non_asy = rate('recovery'), rate('worsening')
    yv = np.asarray(sim.t.yearvec, dtype=float)
    non = np.stack([by_age_sex(r.non_infectious, s)[0] for s in tbode4.sexes], axis=1) # [time, sex, age group]
    lows = by_age_sex(r.non_infectious, tbode4.sexes[0])[1]
    idx = np.clip(np.digitize(lows, tbode4.wpp_edges) - 1, 0, len(tbode4.wpp_edges) - 1)
    mxt = np.array([wpp.mx_at(y)[:, idx] for y in yv])
    tot = non.sum(axis=(1, 2))
    mu = np.divide((non*mxt).sum(axis=(1, 2)), tot, out=np.zeros(len(yv)), where=tot > 0)
    dt = sim.t.dt_year
    inflow = get('latent_non') + get('remote_non')
    n0 = np.zeros(len(inflow)); first_non = np.zeros(len(inflow))
    n0[0] = (1 - share_flagged)*get('non_infectious')[0]
    for i in range(len(inflow)): # Index i holds the step from t_i to t_(i+1)
        rr = non_rec + non_asy + mu[i]
        decay = np.exp(-rr*dt)
        a = (1 - decay)/rr
        b = (dt - a)/(rr*dt)
        first_non[i] = non_asy*(n0[i]*a + inflow[i]*b)
        if i + 1 < len(inflow):
            n0[i+1] = n0[i]*decay + inflow[i]*(1 - decay)/(rr*dt)
    return get('latent_asy') + get('remote_asy') + first_non + get('relapse')


def per_step(sim, first, a5, p5, adult):
    """ Each output per step as fastode_age's per_step(): per person (per year for flows); stocks at index i are at t_i, flows the step from t_i """
    out = tbode4.tx_outputs(sim, first)
    r = sim.results.tb
    alive = np.asarray(sim.results.alive.values, dtype=float)
    v = lambda k: np.asarray(r[k].values, dtype=float) if k in r else np.zeros(len(alive))
    out.prop_asymptomatic = v('asymptomatic')/(v('asymptomatic') + v('symptomatic'))
    out.recent, out.remote, out.cleared = v('infection')/alive, v('infection_remote')/alive, v('cleared')/alive
    for k, g in enumerate(groups5):
        out[f'prev_{g}'] = p5[:, k]/a5[:, k]
        out[f'share_{g}'] = a5[:, k]/a5.sum(axis=1)
    out['prev_15+'] = adult[:, 1]/adult[:, 0]
    out['prop_asymptomatic_15+'] = adult[:, 2]/adult[:, 1]
    out.alive = alive/alive[0]
    out.tx_start = v('tx_start')/alive # Per person per step, for the targets
    out.deaths = (v('tb_death') + v('tb_death_tx'))/alive
    out.first = np.asarray(first, dtype=float)/alive
    out.pop_check = float(np.max(np.abs(a5.sum(axis=1)/alive - 1))) # The groups' sum against everyone alive
    return out


def target_values(out, over_detection):
    """ Exp 32's eight run-based targets (calib.outputs(): 15+ prevalence wherever the survey enters; point values at the step nearest 1 January,
    totals over calendar years), and two reported values """
    year = out.year
    at = lambda y: int(np.argmin(np.abs(year - y)))
    total = lambda k, y: out[k][(year >= y - 1e-9) & (year < y + 1 - 1e-9)].sum() # Flows indexed by their step's start
    prev15 = lambda y: out['prev_15+'][at(y)]
    inc = {y: total('first', y) for y in [2015, 2020, 2024]}
    deaths = {y: total('deaths', y) for y in [2015, 2024]}
    notif = (1 + over_detection)*total('tx_start', 2020)
    return dict(prevalence=prev15(2020), prop_asymptomatic=out['prop_asymptomatic_15+'][at(2020)], prev_to_inc=prev15(2020)/inc[2020],
                mort_to_inc=deaths[2024]/inc[2024], prev_to_notif=prev15(2020)/notif if notif > 0 else np.nan, inc_ratio=inc[2024]/inc[2015],
                mort_ratio=deaths[2024]/deaths[2015], prev_2040_ratio=prev15(2040)/prev15(2020), prevalence_all_2020=out.prevalence[at(2020)],
                prevalence_15plus_1995=prev15(1995))


#%% Runs

def run_ode(i, dt, w):
    """ v4's ODE for set i at one step and width, saved to ode_parts/ (reused if there) """
    part = ode_parts / f'set{i:02d}_{dt}_{w}.obj'
    if part.exists():
        return sc.load(part)
    T = sc.timer()
    s = sc.load(fn.sets)[i]
    sim = make_sim(s, dt, 'ode', width=w)
    sim.run()
    out = per_step(sim, ode_incidence(sim, s.share_flagged), *ode_by_age(sim))
    res = sc.objdict(set=i, dt=dt, width=w, yearly=fa.yearly(out, sim.t.dt_year), targets=target_values(out, s.trend.over_detection),
                     pop_check=out.pop_check, secs=T.toc(output=True))
    sc.save(part, res)
    return res


def run_odes(jobs):
    ode_parts.mkdir(exist_ok=True)
    return sc.parallelize(run_ode, iterkwargs=jobs, ncpus=max(1, min(ncpus, len(jobs))))


def exact_ode(parts):
    """ Exact aging extrapolated (2 × eighth − quarter) from each set's and step's ODE runs: yearly values and target values """
    out = []
    for (i, dt) in sorted({(o.set, o.dt) for o in parts}, key=lambda x: (x[0], list(dts).index(x[1]))):
        q = next(o for o in parts if o.set == i and o.dt == dt and o.width == 'quarter')
        e = next(o for o in parts if o.set == i and o.dt == dt and o.width == 'eighth')
        out.append(sc.objdict(set=i, dt=dt, yearly=2*e.yearly - q.yearly, targets={k: 2*e.targets[k] - q.targets[k] for k in q.targets},
                              quarter=q, eighth=e))
    return out


def mapping_check():
    """ For each set: v4's ODE (exact aging extrapolated, monthly) ÷ fastode_age.run_cohort() from the same table, yearly means 1951–2040 """
    T = sc.timer()
    sets = sc.load(fn.sets)
    parts = run_odes([dict(i=s.set, dt='month', w=w) for s in sets for w in widths])
    keys = ['prevalence', 'prev_15+', 'incidence', 'mortality', 'notifications', 'on_treatment', 'prop_asymptomatic', 'alive'] + [f'prev_{g}' for g in groups5]
    rows = []
    for e in exact_ode(parts):
        f = sets[e.set].fast_yearly
        yrs = f.index[(f.index >= start + 1) & (f.index < stop)]
        row = dict(set=e.set)
        for k in keys:
            r = e.yearly.loc[yrs, k]/f.loc[yrs, k]
            if k in ['notifications', 'on_treatment']:
                r = r[r.index >= 1996] # Treatment starts in 1995
            pre, post = r[r.index < 1995], r[r.index >= 1995]
            row[f'{k} pre'] = float(np.nanmax(np.abs(pre - 1))) if len(pre) else np.nan
            row[f'{k} 1995+'] = float(np.nanmax(np.abs(post - 1)))
        tg = {k: e.targets[k]/sets[e.set].fast_targets[k] - 1 for k in run_targets}
        row['targets max'] = float(np.nanmax(np.abs(list(tg.values()))))
        rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(outdir / 'mapping.csv', index=False)
    worst = df.drop(columns='set').max()
    lines = ["Mapping check: v4's ODE (monthly, exact aging = 2 × eighth − quarter) ÷ fastode_age.run_cohort() (monthly, the calibration's model), "
             'from the same settled table. Largest |ratio − 1| over yearly means, 1951–1994 (pre) and 1995–2040 (notifications and on treatment 1996–2040); '
             "'targets max': the largest |ratio − 1| of the eight targets",
             df.set_index('set').map(lambda v: f'{v:.1e}').to_string(),
             '\nLargest over sets:', worst.map(lambda v: f'{v:.1e}').to_string(),
             f'\nODE run times: ' + ', '.join(f'set {o.set} {o.width} {o.secs:.0f} s' for o in parts) + f'; total {T.toc(output=True):.0f} s']
    text = '\n'.join(lines)
    sc.savetext(outdir / 'mapping.log', text)
    print(text, flush=True)
    return df


def run_all_odes():
    T = sc.timer()
    sets = sc.load(fn.sets)
    parts = run_odes([dict(i=s.set, dt=dt, w=w) for s in sets for dt in dts for w in widths])
    sc.save(fn.ode, parts)
    tbode4.save_versions(outdir, 'v4')
    print(f'{len(parts)} ODE runs ({T.toc(output=True):.0f} s); run times ' + ', '.join(f'{dt} {w} {np.mean([o.secs for o in parts if o.dt == dt and o.width == w]):.0f} s'
                                                                             for dt in dts for w in widths), flush=True)
    return parts


def seed_of(i, dt, rep):
    """ New seeds for this experiment: 34000 + 100·(2·set + step) + rep """
    return 34000 + 100*(2*i + list(dts).index(dt)) + rep


def run_agents(i, dt, rep):
    T = sc.timer()
    s = sc.load(fn.sets)[i]
    sim = make_sim(s, dt, 'abm', rand_seed=seed_of(i, dt, rep))
    sim.run()
    inc, chk, ac = sim.analyzers[0], sim.analyzers[1], sim.analyzers[2]
    adult = np.stack([ac.alive[:, 1:].sum(axis=1), ac.active[:, 1:].sum(axis=1), ac.asym[:, 1:].sum(axis=1)], axis=1) # Groups 15–29 onwards
    out = per_step(sim, inc.first, ac.alive, ac.active, adult)
    alive = np.asarray(sim.results.alive.values, dtype=float)
    return sc.objdict(set=i, dt=dt, rep=rep, seed=seed_of(i, dt, rep), yearly=fa.yearly(out, sim.t.dt_year), targets=target_values(out, s.trend.over_detection),
                      checks=dict(chk.c), course_steps=dict(chk.course_steps), pop_check=out.pop_check, dt_year=sim.t.dt_year,
                      alive_range=(float(alive.min()), float(alive.max())), secs=T.toc(output=True))


def run_all_agents():
    """ Every agent run, in one pool; each set and step saved to runs_by_group/ as it is collected, so a stopped run resumes """
    T = sc.timer()
    sets = sc.load(fn.sets)
    agent_parts.mkdir(exist_ok=True)
    groups = [(s.set, dt) for dt in dts for s in sets]
    part = lambda i, dt: agent_parts / f'set{i:02d}_{dt}.obj'
    todo = [(i, dt) for i, dt in groups if not part(i, dt).exists()]
    if todo:
        jobs = [dict(i=i, dt=dt, rep=r) for i, dt in todo for r in range(reps)]
        print(f'{len(jobs)} agent runs on {min(ncpus, len(jobs))} workers', flush=True)
        res = sc.parallelize(run_agents, iterkwargs=jobs, ncpus=min(ncpus, len(jobs)))
        for i, dt in todo:
            sc.save(part(i, dt), [o for o in res if o.set == i and o.dt == dt])
    runs = sum([sc.load(part(i, dt)) for i, dt in groups], [])
    sc.save(fn.runs, runs)
    tbode4.save_versions(outdir, 'v4')
    print(f'{len(runs)} agent runs done ({T.toc(output=True):.0f} s); ' +
          ', '.join(f'{dt} {np.mean([o.secs for o in runs if o.dt == dt]):.0f} s each' for dt in dts) +
          f'; groups sum to alive within {max(o.pop_check for o in runs):.1e}', flush=True)
    return runs


#%% Summaries

def target_interval(t):
    """ Exp 24's: the target's mean and 95% interval (Jeffreys for a Binomial, ±1.96σ in the log for a log-Normal) """
    if t['kind'] == 'binomial':
        from scipy.stats import beta
        return t['x']/t['n'], *beta.ppf([0.025, 0.975], t['x'] + 0.5, t['n'] - t['x'] + 0.5)
    return t['mean'], t['mean']*np.exp(-1.96*t['sigma']), t['mean']*np.exp(1.96*t['sigma'])


def window_mean(df, w):
    lo, hi = windows[w]
    return df[(df.index >= lo) & (df.index < hi)].mean()


def by_set(runs, ode):
    """ Agents ÷ ODE (exact aging) for each set, step and window: each output, and the per-case outputs ÷ all-ages prevalence; mean ± SE across reps """
    rows = []
    for o in ode:
        rs = [r for r in runs if r.set == o.set and r.dt == o.dt]
        for w in windows:
            ow = window_mean(o.yearly, w)
            ratios = []
            for r in rs:
                aw = window_mean(r.yearly, w)
                d = {k: aw[k]/ow[k] for k in outputs}
                d.update({f'{k}_per_case': (aw[k]/aw.prevalence)/(ow[k]/ow.prevalence) for k in per_case})
                ratios.append(d)
            ratios = pd.DataFrame(ratios)
            for k in ratios:
                rows.append(dict(set=o.set, dt=o.dt, window=w, output=k, ratio=ratios[k].mean(), se=ratios[k].sem(), reps=len(ratios)))
    return pd.DataFrame(rows)


def pooled(bs):
    """ The mean over sets of each set's mean ratio, the SE between sets, min–max, the largest |z| of a set against 1, and sets beyond 2 SE """
    bs = bs.assign(z=(bs.ratio - 1)/bs.se)
    g = bs.groupby(['dt', 'window', 'output'], sort=False)
    return g.agg(mean=('ratio', 'mean'), se_between=('ratio', 'sem'), lo=('ratio', 'min'), hi=('ratio', 'max'), sets=('ratio', 'size'),
                 max_abs_z=('z', lambda x: x.abs().max()), sets_beyond_2se=('z', lambda x: int((x.abs() > 2).sum()))).reset_index()


def target_table(runs, ode):
    """ Each set's target values: the ODE (exact aging), agents (mean ± SE over reps), agents ÷ ODE, and whether each is inside the target's interval """
    tg = cb32().targets
    rows = []
    for o in ode:
        rs = pd.DataFrame([r.targets for r in runs if r.set == o.set and r.dt == o.dt])
        for k in run_targets + target_extra:
            mean, lo, hi = target_interval(tg[k]) if k in tg else (np.nan, np.nan, np.nan)
            a, ov = rs[k].mean(), o.targets[k]
            rows.append(dict(set=o.set, dt=o.dt, target=k, ode=ov, agents=a, agents_se=rs[k].sem(), ratio=a/ov, ratio_se=rs[k].sem()/ov,
                             target_mean=mean, target_lo=lo, target_hi=hi, ode_in=lo <= ov <= hi, agents_in=lo <= a <= hi))
    return pd.DataFrame(rows)


def pitfall_checks(runs):
    """ Exp 23's pitfall checks, pooled over sets and reps for each step """
    rows = []
    for dt in dts:
        rs = [r for r in runs if r.dt == dt]
        if not rs:
            continue
        c = sc.objdict({k: sum(r.checks[k] for r in rs) for k in rs[0].checks})
        steps = sc.ddict(int)
        for r in rs:
            for k, v in r.course_steps.items():
                steps[int(k)] += v
        dty = rs[0].dt_year
        mode_steps = max(steps, key=steps.get)
        rows.append(dict(dt=dt, starts=c.starts, start_while_on_tx=c.start_while_on_tx, start_not_symptomatic=c.start_not_symptomatic,
                         relapse_not_from_success=c.relapse_not_from_success, tb_death_from_success=c.tb_death_from_success,
                         course_years_mean=c.course_steps_sum/c.courses_timed*dty, share_courses_at_mode=steps[mode_steps]/c.courses_timed,
                         mode_years=mode_steps*dty))
    return pd.DataFrame(rows)


main_outputs = ['prevalence', 'prev_15+', 'prop_asymptomatic', 'prop_asymptomatic_15+', 'incidence_per_case', 'notifications_per_case',
                'on_treatment_per_case', 'mortality_per_case']
age_outputs = [f'prev_{g}' for g in groups5] + ['alive'] + [f'share_{g}' for g in groups5]
latent_outputs = ['recent', 'remote', 'cleared', 'mortality_tx_per_case']
window_colors = dict(pre='#7f7f7f', early='#1f77b4', late='#ff7f0e', const='#2ca02c')


def order_by_ridge(sets):
    return list(np.argsort([s.theta['k_latent'] for s in sets]))


def set_label(s):
    return f"{s.theta['k_latent']:.1f}\n{s.theta['inf_asy']:.2f}\n{1e5*s.fast_targets['prevalence_15plus_1995']:.0f}"


def plot_by_set(bs, sets, keys, name, title):
    """ Agents ÷ ODE by set and window for each output, at each step; sets ordered along the ridge, labelled by k_latent, inf_asy and 15+ prevalence in 1995 """
    order = order_by_ridge(sets)
    ncol = len(keys)
    fig, axes = plt.subplots(len(dts), ncol, figsize=(2.9*ncol + 1, 7.5), sharey='row', squeeze=False)
    for row, dt in enumerate(dts):
        for col, k in enumerate(keys):
            ax = axes[row, col]
            for j, w in enumerate(windows):
                d = bs[(bs.dt == dt) & (bs.output == k) & (bs.window == w)].set_index('set').reindex([sets[i].set for i in order])
                x = np.arange(len(order)) + (j - 1.5)*0.18
                ax.errorbar(x, d.ratio, yerr=2*d.se, fmt='o', ms=2.5, capsize=1, lw=0.8, color=window_colors[w], label=f'{windows[w][0]}–{windows[w][1] - 1}')
            ax.axhline(1, color='k', lw=0.6)
            ax.axhspan(0.99, 1.01, color='0.92', zorder=0)
            ax.set_xticks(range(len(order)), [set_label(sets[i]) for i in order], fontsize=5)
            if row == 0:
                ax.set_title(k, fontsize=9)
            if col == 0:
                ax.set_ylabel(f'{dt.capitalize()}ly steps\nagents ÷ ODE (mean ± 2 SE, {reps} reps)', fontsize=8)
    axes[0, 0].legend(fontsize=6, frameon=False)
    fig.supxlabel('Sets ordered by k_latent (labels: k_latent; inf_asy; 15+ prevalence per 100,000 in 1995); grey band ±1%', fontsize=8)
    fig.suptitle(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(outdir / f'{name}.png', dpi=120)
    plt.close(fig)
    return


def plot_over_time(runs, ode, sets, dt):
    """ Prevalence (all ages, 15+), incidence, notifications and mortality, 1950–2040, at three sets along the ridge: reps, their mean and the ODE """
    order = order_by_ridge(sets)
    pick = sorted({order[0], order[len(order)//2], order[-1]}, key=order.index)
    keys = sc.objdict({'prevalence': 'Prevalence, all ages', 'prev_15+': 'Prevalence, 15+', 'incidence': 'Incidence per year',
                       'notifications': 'Treatment starts per year', 'mortality': 'TB mortality per year'})
    fig, axes = plt.subplots(len(pick), len(keys), figsize=(19, 3.2*len(pick)), sharex=True, squeeze=False)
    for r_i, i in enumerate(pick):
        o = next(x for x in ode if x.set == sets[i].set and x.dt == dt).yearly
        for c_i, (k, label) in enumerate(keys.items()):
            ax = axes[r_i, c_i]
            d = pd.concat([x.yearly[k] for x in runs if x.set == sets[i].set and x.dt == dt], axis=1)
            ax.plot(d.index, d.values, color='#7fa7d6', lw=0.5, alpha=0.4)
            ax.plot(d.index, d.mean(axis=1), color='#1f4e8c', lw=1.6, label=f'Agents (mean of {d.shape[1]})')
            ax.plot(o.index, o[k], color='k', lw=1.1, ls='--', label='ODE, exact aging')
            for y in [1995, 2024]:
                ax.axvline(y, color='k', lw=0.5, ls=':')
            ax.set_ylim(bottom=0)
            if r_i == 0:
                ax.set_title(f'{label} (per 100,000)', fontsize=9)
            if c_i == 0:
                ax.set_ylabel(f"Set {sets[i].set}\nk_latent {sets[i].theta['k_latent']:.1f}, inf_asy {sets[i].theta['inf_asy']:.2f}", fontsize=8)
    axes[0, 0].legend(fontsize=7, frameon=False)
    fig.suptitle(f'Exp 34: three sets along the ridge, {dt}ly steps, {n_agents:,} agents at the start', fontsize=10)
    fig.tight_layout()
    fig.savefig(outdir / f'over_time_{dt}.png', dpi=120)
    plt.close(fig)
    return


def plot_targets(tt):
    """ Each set's target values, agents (mean ± 2 SE) and the ODE, against the target's 95% interval, at each step """
    labels = dict(prevalence='15+ prevalence 2020', prop_asymptomatic='15+ share asymptomatic 2020', prev_to_inc='15+ prevalence ÷ incidence 2020 (y)',
                  mort_to_inc='Mortality ÷ incidence 2024', prev_to_notif='15+ prevalence ÷ notifications 2020 (y)', inc_ratio='Incidence 2024 ÷ 2015',
                  mort_ratio='Mortality 2024 ÷ 2015', prev_2040_ratio='15+ prevalence 2040 ÷ 2020')
    fig, axes = plt.subplots(2, 4, figsize=(18, 7.5))
    for ax, k in zip(axes.flat, run_targets):
        t = tt[tt.target == k]
        ax.axhspan(t.target_lo.iloc[0], t.target_hi.iloc[0], color='grey', alpha=0.2, lw=0, label='Target 95% interval')
        ax.axhline(t.target_mean.iloc[0], color='grey', lw=0.8)
        for j, (dt, c) in enumerate(zip(dts, ['#1f4e8c', '#d62728'])):
            d = t[t.dt == dt].sort_values('set')
            x = d.set.values + (j - 0.5)*0.3
            ax.plot(x, d.ode, '_', color='k', ms=7, mew=1.5, label='ODE' if j == 0 else None)
            ax.errorbar(x, d.agents, yerr=2*d.agents_se, fmt='o', ms=2.5, color=c, capsize=1, label=f'Agents, {dt}ly')
        ax.set_title(labels[k], fontsize=9)
        ax.set_xticks(sorted(t.set.unique()))
        ax.set_xlabel('Set', fontsize=8)
    axes.flat[0].legend(fontsize=7, frameon=False)
    fig.suptitle(f'Exp 34: target values by set ({reps} reps per set and step; the ODE with exact aging at the same step)', fontsize=10)
    fig.tight_layout()
    fig.savefig(outdir / 'targets.png', dpi=120)
    plt.close(fig)
    return


def summary():
    sets, runs, parts = sc.load(fn.sets), sc.load(fn.runs), sc.load(fn.ode)
    ode = exact_ode(parts)
    bs = by_set(runs, ode)
    pl = pooled(bs)
    tt = target_table(runs, ode)
    pc = pitfall_checks(runs)
    bs.to_csv(outdir / 'by_set.csv', index=False)
    pl.to_csv(outdir / 'pooled.csv', index=False)
    tt.to_csv(outdir / 'targets.csv', index=False)
    pc.to_csv(outdir / 'pitfall_checks.csv', index=False)
    tp = tt.groupby(['dt', 'target'], sort=False).agg(ratio=('ratio', 'mean'), se_between=('ratio', 'sem'), lo=('ratio', 'min'), hi=('ratio', 'max'),
                                                      ode_in=('ode_in', 'sum'), agents_in=('agents_in', 'sum')).reset_index()
    tp.to_csv(outdir / 'targets_pooled.csv', index=False)
    plot_by_set(bs, sets, main_outputs, 'agents_vs_ode', f'Exp 34: agents ÷ ODE (exact aging) by set and window; {n_agents:,} agents at the start')
    plot_by_set(bs, sets, age_outputs, 'agents_vs_ode_ages', 'Exp 34: agents ÷ ODE by set and window: prevalence by age, population and its age shares')
    for dt in dts:
        plot_over_time(runs, ode, sets, dt)
    plot_targets(tt)
    widths_rows = [] # Quarter- and eighth-year groups ÷ exact aging (deterministic), by window
    for e in ode:
        for w in windows:
            ew = window_mean(e.yearly, w)
            for wd in widths:
                ow = window_mean(e[wd].yearly, w)
                widths_rows.append(dict(set=e.set, dt=e.dt, window=w, width=wd, **{k: ow[k]/ew[k] for k in ['prevalence', 'prev_15+', 'alive', 'prev_65+']}))
    wdf = pd.DataFrame(widths_rows)
    wdf.to_csv(outdir / 'ode_widths.csv', index=False)
    fmt = lambda df: df.to_string(index=False, float_format=lambda v: f'{v:.4f}')
    nrep = sorted({sum(1 for r in runs if r.set == o.set and r.dt == o.dt) for o in ode})
    lines = [f'Exp 34: agents ÷ the ODE with exact aging (2 × eighth − quarter) at the same step. {len(sets)} sets, {nrep} reps each, {n_agents:,} agents at the start. '
             'Window means of yearly means; per_case = each ÷ all-ages prevalence. Windows: ' + ', '.join(f'{k} {a}–{b - 1}' for k, (a, b) in windows.items())]
    with pd.option_context('display.width', 300, 'display.max_columns', 40, 'display.max_rows', 1000):
        for name, keys in [('Main outputs', main_outputs), ('By age and population', age_outputs), ('Latent states and deaths on treatment', latent_outputs)]:
            lines += [f'\n{name}: mean over sets ± SE between sets (min–max of set means; largest |z| of a set against 1; sets beyond 2 SE):']
            for dt in dts:
                d = pl[(pl.dt == dt) & pl.output.isin(keys)].copy()
                d['value'] = d.apply(lambda x: f'{x["mean"]:.4f} ± {x.se_between:.4f} ({x.lo:.3f}–{x.hi:.3f}; z {x.max_abs_z:.1f}; {x.sets_beyond_2se})', axis=1)
                lines += [f'{dt}:', d.pivot(index='output', columns='window', values='value').reindex(index=keys, columns=list(windows)).to_string()]
        lines += ['\nTargets: agents ÷ ODE, mean over sets ± SE between sets (min–max), and sets inside the target\'s 95% interval:', fmt(tp)]
        lines += ['\nThe ODE\'s widths ÷ exact aging (window means; the range over sets):',
                  wdf.groupby(['dt', 'width', 'window'], sort=False)[['prevalence', 'prev_15+', 'alive', 'prev_65+']].agg(['min', 'max']).round(4).to_string()]
        lines += ['\nPitfall checks (pooled over sets and reps):', fmt(pc)]
    lines += [f'\nAgent run times: ' + ', '.join(f'{dt} {np.mean([r.secs for r in runs if r.dt == dt]):.0f} s' for dt in dts) +
              f'; ODE run times: ' + ', '.join(f'{dt} {w} {np.mean([o.secs for o in parts if o.dt == dt and o.width == w]):.0f} s' for dt in dts for w in widths),
              f'Agents alive: {min(r.alive_range[0] for r in runs):,.0f} to {max(r.alive_range[1] for r in runs):,.0f}; '
              f'groups sum to alive within {max(r.pop_check for r in runs):.1e} (agents), {max(o.pop_check for o in parts):.1e} (ODE)']
    text = '\n'.join(lines)
    sc.savetext(outdir / 'summary.log', text)
    print(text, flush=True)
    return


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'summary'
    dict(sets=draw_sets, check=mapping_check, ode=run_all_odes, agents=run_all_agents, summary=summary)[cmd]()
