"""
Experiment 35, v4 environment: death and loss to follow-up on treatment. See README.md. Run from the repo root (after data.py):

    PYTHONPATH=v4 v4/.venv/bin/python experiments/35_treatment_outcomes/run.py phases     # the phase width against continuous time (closed form)
    PYTHONPATH=v4 v4/.venv/bin/python experiments/35_treatment_outcomes/run.py schedules  # stage 2's sets and their carve-out schedules
    PYTHONPATH=v4 v4/.venv/bin/python experiments/35_treatment_outcomes/run.py check      # fastode_age against v4's ODE with the new treatment
    PYTHONPATH=v4 v4/.venv/bin/python experiments/35_treatment_outcomes/run.py ode        # stage 2: v4's ODE, both widths and steps
    PYTHONPATH=v4 v4/.venv/bin/python experiments/35_treatment_outcomes/run.py agents     # stage 2: agents (saved by set, arm and step)
    PYTHONPATH=v4 v4/.venv/bin/python experiments/35_treatment_outcomes/run.py summary    # stage 2: tables and figures
    PYTHONPATH=v4 v4/.venv/bin/python experiments/35_treatment_outcomes/run.py targets    # stage 3: exp 32's posterior draws, five ways (fastode_age)

Environment variables: EXP35_NCPUS (worker processes; default tbode4.ncpus), EXP35_QUICK=1 (a small test of every stage, written to quick/).

Stage 2 follows exp 34's setup (India's WPP demography from 1950, each set's settled start by single year of age and sex, 100,000 agents), with
TBTreat's death and loss to follow-up on treatment (tbode4.treat; README "The design") in place of exp 23's treatment.
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
quick = os.environ.get('EXP35_QUICK') == '1'
ncpus = int(os.environ.get('EXP35_NCPUS', tbode4.ncpus))
outdir = here / 'quick' if quick else here
outdir.mkdir(exist_ok=True)
exp32, exp34 = root / 'experiments' / '32_demography_calibration', root / 'experiments' / '34_calibrated_demography_agents'
start, stop = 1950, 2041
h = 1/12 # fastode_age's step (monthly), as exp 32
tx_dur, phase_dur = 0.5, 1/6 # A 6-month course in 2-month phases
knot_years = [1995, 1999, 2003, 2014, 2023, 2024] # The carve-out is solved at efficacy's knots (1995, 2024) and the loss schedule's (data.py), linear between
n_agents, reps = 100_000, 20
dts = sc.objdict(month=ss.months(1), week=ss.weeks(1))
widths = sc.objdict(quarter=0.25, eighth=0.125)
agent_edges = tbode4.age_edges(1.0)
table_lows = np.arange(101, dtype=float)
tol = 1e-3
groups5 = fa.groups5
windows = sc.objdict(pre=(1950, 1995), early=(1995, 2010), late=(2010, 2025), const=(2025, 2041))
arms2 = sc.objdict(A=dict(omega=0.5, q='f', tx_from=False), B=dict(omega=0.5, q='f', tx_from=True)) # Stage 2
ways3 = sc.objdict(base=dict(omega=0.0, q='f', phase=1/6), q0=dict(omega=0.0, q=0.0, phase=1/6), front=dict(omega=0.5, q='f', phase=1/6),
                   fine=dict(omega=0.0, q='f', phase=1/12)) # Stage 3, besides the current model
n_draws3 = 200 # Per IMIS seed
asym_share, susc_rate, latent_rate, clear_latent = 0.1, 0.001, 0.01, 0.9 # Arm B's placeholder starts: from asymptomatic TB at a tenth of symptomatic TB's rate,
# from susceptible at 0.001 per year, from recent latent infection at 0.01 per year with a 0.9 chance of clearing it (a trial-like efficacy)
if quick:
    n_agents, reps, n_draws3, stop = 10_000, 2, 3, 2041
fn = sc.objdict(sets=outdir / 'sets.obj', ode=outdir / 'ode.obj', runs=outdir / 'runs.obj', targets=outdir / 'targets.obj')
ode_parts, agent_parts = outdir / 'ode_parts', outdir / 'runs_by_group'
wpp = tbode4.WPP()
ss.solvers.Engine.group_bytes = 3e10 # Results by fine age group at every step; exp 35's treatment states and flows take about 2.5× exp 34's
ncpus_ode = int(os.environ.get('EXP35_ODE_NCPUS', min(ncpus, 12))) # The ODEs at weekly steps and eighth-year ages take up to ~10 GB each
tbode4.use_drawn_start()

outputs = ['prevalence', 'prev_15+', 'prop_asymptomatic', 'prop_asymptomatic_15+', 'incidence', 'notifications', 'on_treatment', 'mortality', 'mortality_tx',
           'mortality_tx_cure', 'ltfu', 'ltfu_cured', 'relapse', 'recent', 'remote', 'cleared', 'alive'] + [f'prev_{g}' for g in groups5] + [f'share_{g}' for g in groups5]
per_case = ['incidence', 'notifications', 'on_treatment', 'mortality', 'mortality_tx', 'mortality_tx_cure', 'ltfu', 'ltfu_cured']
target_keys = ['prevalence', 'prop_asymptomatic', 'prev_to_inc', 'mort_to_inc', 'prev_to_notif', 'inc_ratio', 'mort_ratio', 'prev_2040_ratio']


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_cb32 = None


def cb32():
    """ Exp 32's calib.py (parameter mapping, start, targets, likelihood), loaded on first use """
    global _cb32
    if _cb32 is None:
        _cb32 = load_module('calib32', exp32 / 'calib.py')
    return _cb32


#%% WHO's outcomes and the carve-out

def who_knots():
    """ d (died) and ℓ (lost) at the knot years, from data.py's fitted schedule (outcome_schedule.csv: died constant, lost linear between fitted knots) """
    sch = pd.read_csv(here / 'outcome_schedule.csv').set_index('year')
    return pd.DataFrame({k: [float(np.interp(y, sch.index.values, sch[k].values)) for y in knot_years] for k in ['died', 'lost']}, index=knot_years)


def q_rule(q):
    return (lambda f: f) if q == 'f' else q


def efficacy_at(trend, y):
    """ The calibrated overall cure among starts, s(y): exp 24's schedule (linear 1995–2024, constant after) """
    return float(np.interp(y, [1995.0, 2024.0], [trend.efficacy_1995, trend.efficacy_2024]))


def background_at(out, years):
    """ The age-weighted background death rate of people with symptomatic TB, the mean over each year, from a fastode_age run (its bg_sym) """
    yr = out.year
    return {y: float(np.nanmean(out.bg_sym[(yr >= y - 1e-9) & (yr < y + 1 - 1e-9)])) for y in years}


def carve_schedules(trend, bg, omega, q, phase, dt=h):
    """ The carve-out at each knot (tbode4.carve_out), for a course in phases of whole steps of dt: c, mu_f, mu_c (by phase), h, and the solution """
    kn = who_knots()
    widths_ph = tbode4.phase_lengths(tx_dur, phase, dt)
    M = len(widths_ph)
    ends = np.cumsum(widths_ph)
    mid = [float((e - w/2)/sum(widths_ph)) for e, w in zip(ends, widths_ph)]
    qs = [q_rule(q)(f) for f in mid] if q == 'f' else [q]*M
    sol = {y: tbode4.carve_out(efficacy_at(trend, y), kn.loc[y, 'died'], kn.loc[y, 'lost'], widths_ph, q=qs, omega=omega, bg=bg[y]) for y in knot_years}
    return sc.objdict(widths=widths_ph, q=qs, sol=sol, c={y: sol[y].c for y in knot_years}, mu_f={y: sol[y].mu_f for y in knot_years},
                      mu_c=[{y: sol[y].mu_c[k] for y in knot_years} for k in range(M)], h={y: sol[y].h for y in knot_years})


def interp(knots):
    """ A function of the year from {year: value} (linear between, constant outside) """
    x = np.array(sorted(knots), dtype=float)
    v = np.array([knots[k] for k in sorted(knots)], dtype=float)
    return lambda t: float(np.interp(t, x, v))


#%% The phase width against continuous time (closed form)

def run_phases():
    """
    The expected share of a cure-track course lost to follow-up and cured, with q(f) = f, for phases of 6, 3, 2 and 1 months (q at each phase's
    middle) against continuous time (∫ h e^(−λt) t/D dt over the course, λ = h + background), at each knot's carve-out h (a typical set)
    """
    cb = cb32()
    draws = pd.read_csv(exp32 / 'real_draws.csv')
    _, trend = cb.to_pars(draws[cb.names].values[0])
    rows = []
    for y in knot_years:
        sch = carve_schedules(trend, {k: 0.01 for k in knot_years}, 0.0, 'f', phase_dur)
        hh, bg, D = sch.h[y], 0.01, tx_dur
        lam = hh + bg
        tt = np.linspace(0, D, 20001)
        exact = float(np.trapezoid(hh*np.exp(-lam*tt)*tt/D, tt))
        for w in [6, 3, 2, 1]:
            wl = tbode4.phase_lengths(D, w/12, 1/12)
            o = tbode4.outcome_shares(1.0, 0.0, 0.0, hh, [((np.cumsum(wl)[k] - wl[k]/2)/D) for k in range(len(wl))], wl, bg)
            rows.append(dict(knot=y, h=hh, phase_months=w, phases=len(wl), lost_cured=o.cured_ltfu, exact=exact, ratio=o.cured_ltfu/exact,
                             cured=o.cured, cured_exact=float(np.exp(-lam*D) + exact)))
    df = pd.DataFrame(rows)
    df.to_csv(outdir / 'phases.csv', index=False)
    text = ('Partial benefit with q(f) = f on the cure track: the share of courses lost and cured, by phase width, against continuous time (exact), '
            'at the carve-out h of exp 32\'s first real draw (q(f) = f, ω = 0) and a background death rate of 0.01 per year\n' + df.round(5).to_string(index=False))
    sc.savetext(outdir / 'phases.log', text)
    print(text)
    return df


#%% Stage 2: sets, schedules and the model in v4

def pick_sets():
    """ 3 of exp 34's sets along the ridge (lowest, middle, highest k_latent), leaving out its set 0 (a badly fitting prior draw; exp 34) """
    s34 = sc.load(exp34 / 'sets.obj')
    ok = [s for s in s34 if s.set != 0]
    order = sorted(ok, key=lambda s: s.theta['k_latent'])
    return [order[0], order[len(order)//2], order[-1]] if not quick else [order[len(order)//2]]


def make_schedules():
    """ Stage 2's sets, each with its current-model fastode_age run (for the background death rate) and the carve-out schedules for arm A/B,
    at monthly and weekly steps (phases of whole steps) """
    T = sc.timer()
    cb = cb32()
    sets = []
    for s in pick_sets():
        p = cb.model_pars(s.nat, s.trend)
        base = fa.run_cohort(p, h, start, stop, table=s.table, share_flagged=s.share_flagged)
        bg = background_at(base, knot_years)
        s = sc.dcp(s)
        s.bg = bg
        s.sched = sc.objdict({dt: carve_schedules(s.trend, bg, arms2.A['omega'], arms2.A['q'], phase_dur, dts[dt].years) for dt in dts})
        sets.append(s)
        for y in knot_years:
            sol = s.sched.month.sol[y]
            print(f'set {s.set} {y}: s = {efficacy_at(s.trend, y):.3f}, bg {bg[y]:.4f}; c = {sol.c:.3f}{" (capped; raw " + f"{sol.c_raw:.3f})" if sol.capped else ""}, '
                  f'mu_f = {sol.mu_f:.3f}, mu_c1 = {sol.mu_c[0]:.3f}, h = {sol.h:.3f} (sym_dead {np.interp(y, [1995, 2024], [s.trend.sym_dead_1995, s.trend.sym_dead_2024]):.3f})', flush=True)
    sc.save(fn.sets, sets)
    tbode4.save_versions(outdir, 'v4')
    print(f'Schedules done ({T.toc(output=True):.0f} s)', flush=True)
    return sets


def annual_asfr():
    years = np.arange(1950, 2101)
    return pd.DataFrame([dict(Time=y, AgeGrp=a, ASFR=1000*v) for y in years for a, v in zip(wpp.fert_ages, wpp.asfr_at(y))])


def tb_model(s, dt, arm):
    """ TBTreat for set s with exp 35's treatment (arm A, or B with starts from asymptomatic TB and susceptible too), mu = 0, from the set's table """
    t, sch = s.trend, s.sched[dt]
    nat = {k: v for k, v in s.nat.items() if k not in ['sym_dead', 'mu']}
    tx = tbode4.lai_tpt_treatment(row=dict(tx_covg_2010=t.tx_covg_2010, tx_covg_2024=t.tx_covg_2024, efficacy_1995=t.efficacy_1995,
                                           efficacy_2024=t.efficacy_2024, p_relapse=t.p_relapse), step=dts[dt].years, end=stop)
    sym_dead = ss.peryear(ss.timeseries({1995.0: t.sym_dead_1995, 2024.0: t.sym_dead_2024}, interpolation='linear'))
    tx_from = None
    if arms2[arm]['tx_from']:
        times = np.round(np.arange(1995, stop + dts[dt].years/2, dts[dt].years), 6)
        asym = ss.timeseries({float(y): ss.peryear(asym_share*float(tbode4.theta_from_coverage(tx.coverage(y)))) for y in times}, interpolation='linear')
        tx_from = {'symptomatic': tx.tx_rate, 'asymptomatic': asym, 'infection': {1994.99: 0.0, 1995.0: latent_rate}, 'susceptible': {1994.99: 0.0, 1995.0: susc_rate}}
    return tbode4.TBTreat(init=tbode4.init_table(s.table, table_lows), course='fixed', tx_rate=tx.tx_rate, tx_from=tx_from, tx_success=dict(sch.c),
                          tx_p_relapse=t.p_relapse, dur_relapse_tx=1.5, tx_phase_dur=sch.widths, tx_death_fail=dict(sch.mu_f),
                          tx_death_cure=[dict(m) for m in sch.mu_c], tx_ltfu=dict(sch.h), tx_ltfu_cure=sch.q, tx_clear_latent=clear_latent, sym_dead=sym_dead, mu=0, **nat)


class AgeCounts(ss.Analyzer):
    """ Exp 34's: agents alive, with active TB and asymptomatic TB in exp 28's five age groups (ages + tol), leaving out those dying in the step """
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
        def living(u):
            u = np.asarray(u).astype(int)
            return u[~((ppl.ti_dead.raw[u] <= self.ti) | (ppl.ti_removed.raw[u] <= self.ti))] if i > 0 else u
        n = len(self.lows)
        asy, sym = living(tb.asymptomatic.uids), living(tb.symptomatic.uids)
        self.alive[i] = np.bincount(g(living(ppl.auids)), minlength=n)
        self.asym[i] = np.bincount(g(asy), minlength=n)
        self.active[i] = self.asym[i] + np.bincount(g(sym), minlength=n)
        return

    def step(self):
        self._count(self.ti + 1)
        return


def make_sim(s, dt, solver, arm, width='quarter', rand_seed=0):
    abm = solver == 'abm'
    kw = dict(people=ss.population('India', year=start, n_agents=n_agents if abm else None), total_pop=n_agents if abm else 1.0)
    if abm:
        kw.update(rand_seed=rand_seed, analyzers=[tbode4.Incidence(share_flagged=s.share_flagged), tbode4.TxChecks(), AgeCounts()])
        demog = [ss.Pregnancy(fertility_rate=annual_asfr(), rate_per=1000, min_age=10, max_age=55, name='pregnancy'), tbode4.death_rate(wpp)]
    else:
        demog = tbode4.ode_demographics(wpp)
    sim = ss.Sim(tb_model(s, dt, arm), demographics=demog, start=start, stop=stop, dt=dts[dt], solver=solver, verbose=0, **kw)
    sim.stratify(age=agent_edges if abm else tbode4.age_edges(widths[width]), sex=tbode4.sexes)
    return sim


#%% Outputs (exp 34's, with every treatment state counted)

def by_age_sex(res, sex=None):
    x = res.select(sex=sex) if sex is not None else res
    b = x.by('age')
    return np.asarray(b.values, dtype=float), np.array(ss.data.age_edges(b.columns), dtype=float)


def all_states(tb):
    return list(tbode4.tb_states) + [k for k in tb.tx_info.states if k not in tbode4.tb_states]


def ode_by_age(sim):
    r, tb = sim.results.tb, sim.diseases.tb
    pop, lows = sum(by_age_sex(r[k])[0] for k in all_states(tb) if k in r), by_age_sex(r.susceptible)[1]
    asy = by_age_sex(r.asymptomatic)[0]
    act = asy + by_age_sex(r.symptomatic)[0]
    g = np.clip(np.digitize(lows + 1e-9, [lo for lo, _ in groups5.values()]) - 1, 0, len(groups5) - 1)
    a5 = np.stack([pop[:, g == k].sum(axis=1) for k in range(len(groups5))], axis=1)
    p5 = np.stack([act[:, g == k].sum(axis=1) for k in range(len(groups5))], axis=1)
    adult = lows >= 15 - 1e-9
    return a5, p5, np.stack([pop[:, adult].sum(axis=1), act[:, adult].sum(axis=1), asy[:, adult].sum(axis=1)], axis=1)


def ode_incidence(sim, share_flagged):
    """ Exp 34's: the ODE's first entries to asymptomatic TB per step (plus relapses), NON0 carried forward with its background deaths """
    r, tb = sim.results.tb, sim.diseases.tb
    get = lambda k: np.asarray(r[k].values, dtype=float) if k in r else np.zeros(len(r['symptomatic'].values))
    rate = lambda name: float(getattr(tb, name).rate.raw.value)
    non_rec, non_asy = rate('recovery'), rate('worsening')
    yv = np.asarray(sim.t.yearvec, dtype=float)
    non = np.stack([by_age_sex(r.non_infectious, s)[0] for s in tbode4.sexes], axis=1)
    lows = by_age_sex(r.non_infectious, tbode4.sexes[0])[1]
    idx = np.clip(np.digitize(lows, tbode4.wpp_edges) - 1, 0, len(tbode4.wpp_edges) - 1)
    mxt = np.array([wpp.mx_at(y)[:, idx] for y in yv])
    tot = non.sum(axis=(1, 2))
    mu = np.divide((non*mxt).sum(axis=(1, 2)), tot, out=np.zeros(len(yv)), where=tot > 0)
    dt = sim.t.dt_year
    inflow = get('latent_non') + get('remote_non')
    n0 = np.zeros(len(inflow)); first_non = np.zeros(len(inflow))
    n0[0] = (1 - share_flagged)*get('non_infectious')[0]
    for i in range(len(inflow)):
        rr = non_rec + non_asy + mu[i]
        decay = np.exp(-rr*dt)
        a = (1 - decay)/rr
        b = (dt - a)/(rr*dt)
        first_non[i] = non_asy*(n0[i]*a + inflow[i]*b)
        if i + 1 < len(inflow):
            n0[i+1] = n0[i]*decay + inflow[i]*(1 - decay)/(rr*dt)
    return get('latent_asy') + get('remote_asy') + first_non + get('relapse')


def per_step(sim, first, a5, p5, adult):
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
    out.tx_start = out.notifications*sim.t.dt_year
    out.deaths = out.mortality*sim.t.dt_year
    out.first = np.asarray(first, dtype=float)/alive
    out.pop_check = float(np.max(np.abs(a5.sum(axis=1)/alive - 1)))
    return out


def target_values(out, over_detection):
    """ Exp 32's eight run-based targets (15+ prevalence where the survey enters) """
    year = out.year
    at = lambda y: int(np.argmin(np.abs(year - y)))
    total = lambda k, y: out[k][(year >= y - 1e-9) & (year < y + 1 - 1e-9)].sum()
    prev15 = lambda y: out['prev_15+'][at(y)]
    inc = {y: total('first', y) for y in [2015, 2020, 2024]}
    deaths = {y: total('deaths', y) for y in [2015, 2024]}
    notif = (1 + over_detection)*total('tx_start', 2020)
    return dict(prevalence=prev15(2020), prop_asymptomatic=out['prop_asymptomatic_15+'][at(2020)], prev_to_inc=prev15(2020)/inc[2020],
                mort_to_inc=deaths[2024]/inc[2024], prev_to_notif=prev15(2020)/notif if notif > 0 else np.nan, inc_ratio=inc[2024]/inc[2015],
                mort_ratio=deaths[2024]/deaths[2015], prev_2040_ratio=prev15(2040)/prev15(2020))


def yearly(out, dt):
    """ fastode_age's yearly() over the outputs both solvers have """
    year = np.floor(out.year[1:] - dt/2 + 1e-9).astype(int)
    df = pd.DataFrame({k: (out[k][1:] if k in fa.stock_keys else out[k][:-1])*(1.0 if k in fa.unscaled else 1e5) for k in outputs})
    df['year'] = year
    return df.groupby('year').mean()


#%% Stage 2 runs

def run_ode(i, dt, w, arm):
    part = ode_parts / f'set{i:02d}_{arm}_{dt}_{w}.obj'
    if part.exists():
        return sc.load(part)
    T = sc.timer()
    s = next(x for x in sc.load(fn.sets) if x.set == i)
    sim = make_sim(s, dt, 'ode', arm, width=w)
    sim.run()
    out = per_step(sim, ode_incidence(sim, s.share_flagged), *ode_by_age(sim))
    res = sc.objdict(set=i, dt=dt, width=w, arm=arm, yearly=yearly(out, sim.t.dt_year), targets=target_values(out, s.trend.over_detection),
                     pop_check=out.pop_check, secs=T.toc(output=True))
    sc.save(part, res)
    return res


def run_odes(jobs):
    ode_parts.mkdir(exist_ok=True)
    return sc.parallelize(run_ode, iterkwargs=jobs, ncpus=max(1, min(ncpus_ode, len(jobs))))


def exact_ode(parts):
    out = []
    for key in sorted({(o.set, o.arm, o.dt) for o in parts}):
        i, arm, dt = key
        q = next(o for o in parts if (o.set, o.arm, o.dt) == key and o.width == 'quarter')
        e = next(o for o in parts if (o.set, o.arm, o.dt) == key and o.width == 'eighth')
        out.append(sc.objdict(set=i, arm=arm, dt=dt, yearly=2*e.yearly - q.yearly, targets={k: 2*e.targets[k] - q.targets[k] for k in q.targets}))
    return out


def mapping_check():
    """ fastode_age (run_cohort, monthly, with tx_options() from the same schedules) against v4's ODE with exact aging extrapolated, arm A """
    T = sc.timer()
    sets = sc.load(fn.sets)
    cb = cb32()
    parts = run_odes([dict(i=s.set, dt='month', w=w, arm='A') for s in sets for w in widths])
    keys = ['prevalence', 'prev_15+', 'incidence', 'mortality', 'mortality_tx', 'mortality_tx_cure', 'notifications', 'on_treatment', 'ltfu', 'ltfu_cured',
            'alive'] + [f'prev_{g}' for g in groups5]
    rows, fast = [], {}
    for e in exact_ode(parts):
        s = next(x for x in sets if x.set == e.set)
        sch = s.sched.month
        p = cb.model_pars(s.nat, s.trend)
        p['tx_success'] = interp(sch.c)
        tx = fa.tx_options(h, tx_dur=tx_dur, phase_dur=phase_dur, death_fail=interp(sch.mu_f), death_cure=[interp(m) for m in sch.mu_c], ltfu=interp(sch.h),
                           ltfu_cure=sch.q)
        f = fa.run_cohort(p, h, start, stop, table=s.table, share_flagged=s.share_flagged, tx=tx)
        fy = fa.yearly(f, h)
        fast[e.set] = fy
        yrs = fy.index[(fy.index >= start + 1) & (fy.index < stop)]
        row = dict(set=e.set)
        for k in keys:
            r = e.yearly.loc[yrs, k]/fy.loc[yrs, k]
            if k in ['notifications', 'on_treatment', 'mortality_tx', 'mortality_tx_cure', 'ltfu', 'ltfu_cured']:
                r = r[r.index >= 1996]
            pre, post = r[r.index < 1995], r[r.index >= 1995]
            row[f'{k} pre'] = float(np.nanmax(np.abs(pre - 1))) if len(pre) else np.nan
            row[f'{k} 1995+'] = float(np.nanmax(np.abs(post - 1)))
        ft = cb.outputs(sc.objdict(f, trend=s.trend))
        row['targets max'] = float(np.nanmax([abs(e.targets[k]/ft[k] - 1) for k in target_keys]))
        rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(outdir / 'check.csv', index=False)
    sc.save(outdir / 'fastode.obj', fast) # fastode_age's yearly values at each set (arm A, monthly), for agents ÷ fastode_age in the summary
    text = ("fastode_age with exp 35's treatment (arm A: ω = 0.5, q(f) = f) against v4's ODE (monthly, exact aging = 2 × eighth − quarter): largest |v4 ÷ "
            "fastode_age − 1| over yearly means, 1951–1994 (pre) and 1995–2040 (treatment outputs from 1996)\n" + df.set_index('set').map(lambda v: f'{v:.1e}').T.to_string() +
            f'\n({T.toc(output=True):.0f} s)')
    sc.savetext(outdir / 'check.log', text)
    print(text, flush=True)
    return df


def run_all_odes():
    T = sc.timer()
    sets = sc.load(fn.sets)
    parts = run_odes([dict(i=s.set, dt=dt, w=w, arm=a) for s in sets for a in arms2 for dt in dts for w in widths])
    sc.save(fn.ode, parts)
    tbode4.save_versions(outdir, 'v4')
    print(f'{len(parts)} ODE runs ({T.toc(output=True):.0f} s)', flush=True)
    return parts


def seed_of(i, arm, dt, rep):
    return 35000 + 1000*list(arms2).index(arm) + 100*(2*i + list(dts).index(dt)) + rep


def run_agents(i, arm, dt, rep):
    T = sc.timer()
    s = next(x for x in sc.load(fn.sets) if x.set == i)
    sim = make_sim(s, dt, 'abm', arm, rand_seed=seed_of(i, arm, dt, rep))
    sim.run()
    inc, chk, ac = sim.analyzers[0], sim.analyzers[1], sim.analyzers[2]
    adult = np.stack([ac.alive[:, 1:].sum(axis=1), ac.active[:, 1:].sum(axis=1), ac.asym[:, 1:].sum(axis=1)], axis=1)
    out = per_step(sim, inc.first, ac.alive, ac.active, adult)
    return sc.objdict(set=i, arm=arm, dt=dt, rep=rep, seed=seed_of(i, arm, dt, rep), yearly=yearly(out, sim.t.dt_year),
                      targets=target_values(out, s.trend.over_detection), checks=dict(chk.c), course_steps=dict(chk.course_steps), dt_year=sim.t.dt_year,
                      secs=T.toc(output=True))


def run_all_agents():
    T = sc.timer()
    sets = sc.load(fn.sets)
    agent_parts.mkdir(exist_ok=True)
    groups = [(s.set, a, dt) for dt in dts for a in arms2 for s in sets]
    part = lambda i, a, dt: agent_parts / f'set{i:02d}_{a}_{dt}.obj'
    todo = [g for g in groups if not part(*g).exists()]
    if todo:
        jobs = [dict(i=i, arm=a, dt=dt, rep=r) for i, a, dt in todo for r in range(reps)]
        print(f'{len(jobs)} agent runs on {min(ncpus, len(jobs))} workers', flush=True)
        res = sc.parallelize(run_agents, iterkwargs=jobs, ncpus=min(ncpus, len(jobs)))
        for g in todo:
            sc.save(part(*g), [o for o in res if (o.set, o.arm, o.dt) == g])
    runs = sum([sc.load(part(*g)) for g in groups], [])
    sc.save(fn.runs, runs)
    tbode4.save_versions(outdir, 'v4')
    print(f'{len(runs)} agent runs done ({T.toc(output=True):.0f} s); ' + ', '.join(f'{dt} {np.mean([o.secs for o in runs if o.dt == dt]):.0f} s each' for dt in dts), flush=True)
    return runs


#%% Stage 2 summary

def window_mean(df, w):
    lo, hi = windows[w]
    return df[(df.index >= lo) & (df.index < hi)].mean()


def by_set(runs, ode):
    rows = []
    for o in ode:
        rs = [r for r in runs if (r.set, r.arm, r.dt) == (o.set, o.arm, o.dt)]
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
                if w == 'pre' and any(k.startswith(t) for t in ['notifications', 'on_treatment', 'mortality_tx', 'ltfu']):
                    continue # No treatment before 1995 (only the edge of the scale-up, a step apart between solvers)
                rows.append(dict(set=o.set, arm=o.arm, dt=o.dt, window=w, output=k, ratio=ratios[k].mean(), se=ratios[k].sem(), reps=len(ratios)))
    return pd.DataFrame(rows)


def pooled(bs):
    bs = bs.assign(z=(bs.ratio - 1)/bs.se)
    g = bs.groupby(['arm', 'dt', 'window', 'output'], sort=False)
    return g.agg(mean=('ratio', 'mean'), se_between=('ratio', 'sem'), lo=('ratio', 'min'), hi=('ratio', 'max'), sets=('ratio', 'size'),
                 max_abs_z=('z', lambda x: x.abs().max()), sets_beyond_2se=('z', lambda x: int((x.abs() > 2).sum()))).reset_index()


def pitfalls(runs):
    rows = []
    for arm in arms2:
        for dt in dts:
            rs = [r for r in runs if r.arm == arm and r.dt == dt]
            if not rs:
                continue
            c = sc.objdict({k: sum(r.checks[k] for r in rs) for k in rs[0].checks})
            dty = rs[0].dt_year
            rows.append(dict(arm=arm, dt=dt, starts=c.starts, start_while_on_tx=c.start_while_on_tx, start_not_symptomatic=c.start_not_symptomatic,
                             relapse_not_from_success=c.relapse_not_from_success, tb_death_from_success=c.tb_death_from_success, notb_cured=c.notb_cured,
                             tb_death_notb=c.tb_death_notb, latent_cleared=c.latent_cleared, latent_to_relapsing=c.latent_to_relapsing, tb_death_latent=c.tb_death_latent,
                             left_to_origin=c.left_to_origin, left_elsewhere=c.left_elsewhere,
                             died_share=(c.tx_deaths_fail + c.tx_deaths_cure)/c.starts, cure_track_death_share=c.tx_deaths_cure/max(c.tx_deaths_fail + c.tx_deaths_cure, 1),
                             lost_share=c.ltfu/c.starts, lost_cured_share=c.ltfu_cured/max(c.ltfu, 1),
                             completed_course_years=c.course_steps_sum/max(c.courses_timed, 1)*dty))
    return pd.DataFrame(rows)


main_outputs = ['prevalence', 'prev_15+', 'prop_asymptomatic', 'incidence_per_case', 'notifications_per_case', 'on_treatment_per_case', 'mortality_per_case',
                'mortality_tx_per_case', 'mortality_tx_cure_per_case', 'ltfu_per_case', 'ltfu_cured_per_case']


def plot_by_set(bs, sets):
    fig, axes = plt.subplots(len(arms2)*len(dts), len(main_outputs), figsize=(2.6*len(main_outputs) + 1, 2.6*len(arms2)*len(dts)), sharey=True, squeeze=False)
    colors = dict(pre='#7f7f7f', early='#1f77b4', late='#ff7f0e', const='#2ca02c')
    for row, (arm, dt) in enumerate([(a, d) for a in arms2 for d in dts]):
        for col, k in enumerate(main_outputs):
            ax = axes[row, col]
            for j, w in enumerate(windows):
                d = bs[(bs.arm == arm) & (bs.dt == dt) & (bs.output == k) & (bs.window == w)].set_index('set').reindex([s.set for s in sets])
                ax.errorbar(np.arange(len(sets)) + (j - 1.5)*0.15, d.ratio, yerr=2*d.se, fmt='o', ms=2.5, capsize=1, lw=0.8, color=colors[w], label=f'{windows[w][0]}–{windows[w][1] - 1}')
            ax.axhline(1, color='k', lw=0.6); ax.axhspan(0.99, 1.01, color='0.92', zorder=0)
            ax.set_xticks(range(len(sets)), [f"set {s.set}\n{s.theta['k_latent']:.1f}" for s in sets], fontsize=6)
            ax.set_ylim(0.9, 1.1)
            if row == 0:
                ax.set_title(k.replace('_per_case', '\nper case'), fontsize=8)
            if col == 0:
                ax.set_ylabel(f'Arm {arm}, {dt}ly\nagents ÷ ODE (± 2 SE)', fontsize=8)
    axes[0, 0].legend(fontsize=6, frameon=False)
    fig.suptitle(f'Exp 35, stage 2: agents ÷ ODE (exact aging) by set and window; {reps} reps of {n_agents:,} agents; sets labelled by k_latent', fontsize=10)
    fig.tight_layout(); fig.savefig(outdir / 'agents_vs_ode.png', dpi=120); plt.close(fig)
    return


def summary():
    sets, runs, parts = sc.load(fn.sets), sc.load(fn.runs), sc.load(fn.ode)
    ode = exact_ode(parts)
    bs = by_set(runs, ode)
    pl = pooled(bs)
    pc = pitfalls(runs)
    bs.to_csv(outdir / 'by_set.csv', index=False); pl.to_csv(outdir / 'pooled.csv', index=False); pc.to_csv(outdir / 'pitfall_checks.csv', index=False)
    plot_by_set(bs, sets)
    lines = [f'Exp 35, stage 2: agents ÷ the ODE with exact aging, {len(sets)} sets ({", ".join(str(s.set) for s in sets)}), {reps} reps of {n_agents:,} agents. '
             'Mean over sets ± SE between sets (min–max; largest |z|; sets beyond 2 SE). Windows: ' + ', '.join(f'{k} {a}–{b - 1}' for k, (a, b) in windows.items())]
    with pd.option_context('display.width', 300, 'display.max_columns', 40, 'display.max_rows', 1000):
        for arm in arms2:
            for dt in dts:
                d = pl[(pl.arm == arm) & (pl.dt == dt) & pl.output.isin(main_outputs + [f'prev_{g}' for g in groups5] + ['alive'])].copy()
                d['value'] = d.apply(lambda x: f'{x["mean"]:.4f} ± {x.se_between:.4f} ({x.lo:.3f}–{x.hi:.3f}; z {x.max_abs_z:.1f}; {x.sets_beyond_2se})', axis=1)
                lines += [f'\nArm {arm}, {dt}ly:', d.pivot(index='output', columns='window', values='value').reindex(
                    index=main_outputs + [f'prev_{g}' for g in groups5] + ['alive'], columns=list(windows)).to_string()]
        lines += ['\nPitfall checks and outcome shares (agents, pooled over sets and reps):', pc.round(5).to_string(index=False)]
    fast = sc.load(outdir / 'fastode.obj')
    fref = [sc.objdict(set=i, arm='A', dt=dt, yearly=y) for i, y in fast.items() for dt in dts] # The calibration's model (monthly) as the reference for both steps
    bf = by_set([r for r in runs if r.arm == 'A'], fref)
    bf.to_csv(outdir / 'by_set_fastode.csv', index=False)
    pf = pooled(bf)
    show = ['prevalence', 'prev_15+', 'prop_asymptomatic', 'incidence_per_case', 'notifications_per_case', 'mortality_per_case', 'mortality_tx_per_case', 'ltfu_per_case']
    with pd.option_context('display.width', 300, 'display.max_columns', 40):
        for dt in dts:
            d = pf[(pf.dt == dt) & pf.output.isin(show)].copy()
            d['value'] = d.apply(lambda x: f'{x["mean"]:.4f} ± {x.se_between:.4f} ({x.lo:.3f}–{x.hi:.3f}; z {x.max_abs_z:.1f}; {x.sets_beyond_2se})', axis=1)
            lines += [f'\nAgents ({dt}ly, arm A) ÷ fastode_age (monthly, the calibration\'s model) directly:', d.pivot(index='output', columns='window', values='value').reindex(index=show, columns=list(windows)).to_string()]
    lines += [f'\nAgent run times: ' + ', '.join(f'{dt} {np.mean([r.secs for r in runs if r.dt == dt]):.0f} s' for dt in dts)]
    text = '\n'.join(lines)
    sc.savetext(outdir / 'summary.log', text)
    print(text, flush=True)
    return


#%% Stage 3: exp 32's posterior draws, five ways

def run_draw3(seed32, index, theta):
    """ One draw: the current model (fastode_age, exp 32's), then the carve-out four ways (ways3) from the same start, each with its targets """
    cb = cb32()
    nat, trend = cb.to_pars(np.asarray(theta, dtype=float))
    init, flagged = cb.untreated_equilibrium(nat)
    out = sc.objdict(seed32=seed32, draw=int(index), theta=dict(zip(cb.names, map(float, theta))), ways=sc.objdict())
    if init is None:
        return out
    p = cb.model_pars(nat, trend)
    st = fa.settle(p, h, init, flagged, years=100, start=start)
    table, share_flagged = fa.start_table(st)
    prog = cb.progression(nat)
    base = fa.run_cohort(p, h, start, stop, table=table, share_flagged=share_flagged)
    base.trend = trend
    o = cb.outputs(base); o.update(prog)
    out.ways.current = sc.objdict(outputs=o, ll=float(cb.loglik(o)))
    bg = background_at(base, knot_years)
    out.bg = bg
    for name, w in ways3.items():
        sch = carve_schedules(trend, bg, w['omega'], w['q'], w['phase'])
        q = dict(p); q['tx_success'] = interp(sch.c)
        tx = fa.tx_options(h, tx_dur=tx_dur, phase_dur=w['phase'], death_fail=interp(sch.mu_f), death_cure=[interp(m) for m in sch.mu_c], ltfu=interp(sch.h),
                           ltfu_cure=sch.q)
        r = fa.run_cohort(q, h, start, stop, table=table, share_flagged=share_flagged, tx=tx)
        r.trend = trend
        o = cb.outputs(r); o.update(prog)
        sols = {y: dict(c=sch.sol[y].c, c_raw=sch.sol[y].c_raw, capped=sch.sol[y].capped, mu_f=sch.sol[y].mu_f, mu_c1=sch.sol[y].mu_c[0], h=sch.sol[y].h,
                        residual=sch.sol[y].residual, **{f'share_{k}': v for k, v in sch.sol[y].shares.items()}) for y in knot_years}
        out.ways[name] = sc.objdict(outputs=o, ll=float(cb.loglik(o)), sol=sols)
    return out


def run_targets():
    T = sc.timer()
    cb = cb32()
    jobs = []
    for seed32 in ['real', 'real_seed2']:
        d = pd.read_csv(exp32 / f'{seed32}_draws.csv').iloc[:n_draws3]
        jobs += [dict(seed32=seed32, index=int(r['index']), theta=r[cb.names].values.astype(float)) for _, r in d.iterrows()]
    print(f'{len(jobs)} draws on {min(ncpus, len(jobs))} workers', flush=True)
    res = sc.parallelize(run_draw3, iterkwargs=jobs, ncpus=max(1, min(ncpus, len(jobs))))
    sc.save(fn.targets, res)
    tbode4.save_versions(outdir, 'v4')
    print(f'Stage 3 runs done ({T.toc(output=True):.0f} s)', flush=True)
    summarize_targets()
    return res


def summarize_targets():
    res = [r for r in sc.load(fn.targets) if 'current' in r.ways]
    cb = cb32()
    keys = target_keys + ['ProgressionY1', 'ProgressionRatio']
    rows, sols = [], []
    for r in res:
        cur = r.ways.current
        for name in ways3:
            w = r.ways[name]
            rows.append(dict(seed32=r.seed32, draw=r.draw, way=name, d_ll=w.ll - cur.ll, ll=w.ll, ll_current=cur.ll,
                             **{f'{k}_ratio': w.outputs[k]/cur.outputs[k] for k in keys}, **{f'{k}': w.outputs[k] for k in keys}))
            for y, s in w.sol.items():
                sols.append(dict(seed32=r.seed32, draw=r.draw, way=name, knot=y, bg=r.bg[y], efficacy=efficacy_at(cb.to_pars(np.array([r.theta[k] for k in cb.names]))[1], y), **s))
    df, sd = pd.DataFrame(rows), pd.DataFrame(sols)
    df.to_csv(outdir / 'targets.csv', index=False); sd.to_csv(outdir / 'carve_out.csv', index=False)
    q = lambda x: f'{np.median(x):.3f} ({np.percentile(x, 5):.3f}–{np.percentile(x, 95):.3f})'
    lines = [f'Exp 35, stage 3: {len(res)} of exp 32\'s posterior draws ({n_draws3} per seed), fastode_age monthly. Each way ÷ the current model: median (5–95%) over draws']
    t = pd.DataFrame({name: {k: q(df[df.way == name][f'{k}_ratio']) for k in target_keys} for name in ways3})
    t.loc['Δ log-likelihood'] = [q(df[df.way == name].d_ll) for name in ways3]
    lines += [t.to_string()]
    tg = cb.targets
    def interval(k):
        tt = tg[k]
        if tt['kind'] == 'binomial':
            from scipy.stats import beta
            return beta.ppf([0.025, 0.975], tt['x'] + 0.5, tt['n'] - tt['x'] + 0.5)
        return tt['mean']*np.exp(-1.96*tt['sigma']), tt['mean']*np.exp(1.96*tt['sigma'])
    inside = {name: {k: f"{np.mean([(interval(k)[0] <= v <= interval(k)[1]) for v in df[df.way == name][k]]):.2f}" for k in target_keys} for name in ways3}
    cur = {k: f"{np.mean([(interval(k)[0] <= r.ways.current.outputs[k] <= interval(k)[1]) for r in res]):.2f}" for k in target_keys}
    lines += ['\nShare of draws inside each target\'s 95% interval:', pd.DataFrame(dict(current=cur, **inside)).to_string()]
    lines += ['\nThe carve-out by knot (median (5–95%) over draws): c, its uncapped value, the share capped, mu_f, mu_c1, h; the fail track\'s death rate against sym_dead']
    for name in ways3:
        for y in knot_years:
            x = sd[(sd.way == name) & (sd.knot == y)]
            bad = x.residual > 1e-6 # No rates reproduce the shares: with deaths only on the fail track, too few failing courses to hold WHO's deaths
            lines.append(f'{name} {y}: efficacy {q(x.efficacy)}, c {q(x.c)}, c_raw {q(x.c_raw)}, capped {x.capped.mean():.3f}, unsolved {bad.mean():.3f} '
                         f'(their efficacy {q(x.efficacy[bad]) if bad.any() else "-"}), mu_f {q(x.mu_f[~bad])}, mu_c1 {q(x.mu_c1[~bad])}, h {q(x.h)}, bg {q(x.bg)}, '
                         f'cured on being lost {q(x.share_cured_ltfu)}, failed {q(x.share_failed)}')
    fig, axes = plt.subplots(2, 5, figsize=(18, 6.5))
    for ax, k in zip(axes.flat, target_keys + ['d_ll']):
        for name, c in zip(ways3, ['#1f4e8c', '#7f7f7f', '#d62728', '#2ca02c']):
            v = df[df.way == name][k if k == 'd_ll' else f'{k}_ratio']
            ax.hist(v, bins=40, histtype='step', color=c, lw=1.2, label=name)
        ax.axvline(0 if k == 'd_ll' else 1, color='k', lw=0.6)
        ax.set_title('Δ log-likelihood' if k == 'd_ll' else f'{k}: new ÷ current', fontsize=9)
    axes.flat[-1].axis('off')
    axes.flat[0].legend(fontsize=7, frameon=False)
    fig.suptitle(f"Exp 35, stage 3: exp 32's targets with death and loss to follow-up on treatment (carved out of calibrated efficacy) ÷ the current model; {len(res)} posterior draws", fontsize=10)
    fig.tight_layout(); fig.savefig(outdir / 'targets.png', dpi=120); plt.close(fig)
    text = '\n'.join(lines)
    sc.savetext(outdir / 'targets.log', text)
    print(text, flush=True)
    return df


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'summary'
    dict(phases=run_phases, schedules=make_schedules, check=mapping_check, ode=run_all_odes, agents=run_all_agents, summary=summary,
         targets=run_targets, targets_summary=summarize_targets)[cmd]()
