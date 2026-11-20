"""
A standalone copy of v4's TB ODE with ages and demography, for calibration (exp 31): `tbode4.fastode`'s TB model and treatment, on a
population by age and sex with WPP's fertility and mortality, as exp 30 runs v4's ODE (CONVENTIONS.md 10).

The TB model is fastode's (the two-state latent model, every rate the same at every age, homogeneous mixing with exposure the share of
the living population in each state, treatment as CONVENTIONS.md 9 with a fixed course), with `mu` = 0 and no replacement of TB deaths.
Demography: deaths by age and sex at WPP's rates, the same in every TB state; births by women's age group into susceptible in the first
age group, as female as the starting population (as v4's ODE with `joint_sex()`); WPP's rates interpolated linearly between its decades.

Two ways of aging (exp 31's README):

- `run_rate()`: a copy of v4's ODE. Age groups of `width` years, aged out of at rate 1/width (as Starsim's compartmental solver does);
  solved as one system of ODEs (`solve_ivp`); the fixed course as Erlang stages; time-varying inputs in continuous time.
- `run_cohort()`: exact aging, for calibration. Age groups one step h wide, each moved up one group at the end of every step. Within a
  step, the TB transitions are expm(Q·h) applied to every age and sex group at once (Q, the transition rates, is the same at every age),
  with the force of infection at mid-step by a predictor-corrector (`foi='mid'`) or at the step's start (`foi='start'`), and inputs at
  mid-step. Deaths are survival exp(−m·h) by group. The fixed course is a conveyor of whole steps, as in Starsim's ODE: K = the course in
  whole steps (nearest), slots 0..K; people starting treatment enter slot 0 during a step, every slot moves along one at its end, and the
  last is released over the next step (half at its start, half at its end), so a course lasts K steps on average.

First-entry incidence (CONVENTIONS.md 5) is counted exactly: non-infectious TB is NON0 (not yet asymptomatic this episode) and NON1, by age.

Inputs: `fastode.make_pars()` parameters, where tx_rate, tx_success and sym_dead may be functions of calendar year (`schedules()` makes
exp 24's); a `Demography` (WPP rates, frozen or not); a starting state from `initial()` (a table of shares by single year of age and sex,
as exp 30's `start.obj`, on WPP's population). Outputs per step as exp 28's `tb_per_step()` (stocks at t_i, flows over the step from t_i),
`yearly()` as exp 28's yearly means, and `targets()` as exp 24's calibration targets.
"""
import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss
from scipy.linalg import expm
from scipy.integrate import solve_ivp
from . import fastode as fo
from .demog import WPP, wpp_edges, sexes, female_share

__all__ = ['Demography', 'Layout', 'schedules', 'initial', 'run_cohort', 'run_rate', 'settle', 'start_table', 'yearly', 'calendar', 'targets', 'groups5']

base_states = ['susceptible', 'infection', 'infection_remote', 'cleared', 'recovered', 'non0', 'non1', 'asymptomatic', 'symptomatic', 'treated', 'relapsing']
S, L, M, C, R, N0, N1, A, Y, TR, RL = range(len(base_states))
counter_names = ['first', 'tx_start', 'tb_death', 'tb_death_tx', 'relapse'] # Per-step flows: incidence (with relapses), notifications, TB deaths, on failing courses, relapses
groups5 = sc.objdict({'0-14': (0, 15), '15-29': (15, 30), '30-49': (30, 50), '50-64': (50, 65), '65+': (65, 999)}) # Exp 28's age groups
t_tx = 1994 + 364/365 # tb_LAI_TPT's first treatment knot, 1994-12-31


#%% Inputs

class Demography:
    """
    WPP's rates on an age grid (lower edges `lows`, the last open): death rates by sex and fertility by group, at any year.

    Args:
        lows (array): the groups' lower edges (years)
        country (str): WPP country
        frozen (float): if given, every rate is held at its value in this year (e.g. 1950, for the starting state)
        start_year (int): the year whose population sets the female share at birth (as v4's ODE with joint_sex())
        cohort (bool): rates for exact aging by cohort, groups one step wide: over a step, a group's people age from [g, g+1) to [g+1, g+2)
            (in steps), spending exactly half their person-time in each, so each group takes the mean of its own rate and the next group's.
            With WPP's rate boundaries on the grid, this is exact; the lower edge's rate alone applies changes with age half a step late
    """
    def __init__(self, lows, country='India', frozen=None, start_year=1950, cohort=False):
        self.lows = np.asarray(lows, dtype=float)
        self.wpp = WPP(country)
        self.frozen = frozen
        self.jd = np.clip(np.digitize(self.lows + 1e-9, wpp_edges) - 1, 0, len(wpp_edges) - 1) # Each group's WPP mortality group
        fa = np.array(self.wpp.fert_ages, dtype=float)
        jf = np.clip(np.digitize(self.lows + 1e-9, fa) - 1, 0, len(fa) - 1)
        self.fertile = (self.lows >= fa[0] - 1e-9) & (self.lows < fa[-1] + 5 - 1e-9) # Ages 10–54 (by the group's lower edge, as v4's ODE)
        self.jf = np.where(self.fertile, jf, 0)
        self.p_female = female_share(country, start_year)[2]
        self.country, self.start_year, self.cohort = country, start_year, cohort
        return

    def year(self, t):
        return self.frozen if self.frozen is not None else t

    def _over_step(self, v):
        """ With cohort=True, each group's rate over a step: the mean of its own and the next group's (the last, open, keeps its own) """
        if not self.cohort:
            return v
        nxt = np.concatenate([v[..., 1:], v[..., -1:]], axis=-1)
        return 0.5*(v + nxt)

    def mx(self, t):
        """ Death rates [sex, group] in year t """
        return self._over_step(self.wpp.mx_at(self.year(t))[:, self.jd])

    def fert(self, t):
        """ Fertility per woman per year, by group, in year t """
        return self._over_step(np.where(self.fertile, self.wpp.asfr_at(self.year(t))[self.jf], 0.0))

    def mx_first(self, t):
        """ The first group's own death rates [sex] (for births during a step, who stay in it) """
        return self.wpp.mx_at(self.year(t))[:, self.jd[0]]


class Layout:
    """
    The compartments of one age and sex group: the 11 base states, each course's compartments (Erlang stages, or conveyor slots 0..K),
    then the counters (per-step flows). `course` is ('stages', k) or ('slots', K)
    """
    def __init__(self, course):
        self.kind, self.k = course
        n = self.k if self.kind == 'stages' else self.k + 1
        self.succ = np.arange(len(base_states), len(base_states) + n)
        self.fail = np.arange(self.succ[-1] + 1, self.succ[-1] + 1 + n)
        self.ns = int(self.fail[-1] + 1) # People's compartments
        self.cnt = {name: self.ns + i for i, name in enumerate(counter_names)}
        self.n = self.ns + len(counter_names)
        self.names = base_states + [f'succ_{i}' for i in range(n)] + [f'fail_{i}' for i in range(n)] + counter_names
        return


def generator_parts(p, lay):
    """
    The transition rates within one group as matrices, dx/dt = Q x: Q = Qc + λ·Ql + θ·s·Qs + θ·(1 − s)·Qf + d·Qd, with λ the force of
    infection on a susceptible (beta × (trans_asymp × A + Y) ÷ the living), θ the rate of starting treatment, s the chance a course cures, and
    d = sym_dead. Rows for counters accumulate flows. In 'slots' courses, people move along the conveyor outside Q
    """
    n = lay.n
    parts = sc.objdict({k: np.zeros((n, n)) for k in ['c', 'l', 's', 'f', 'd']})
    def flow(Q, src, dst, rate, counters=()):
        Q[src, src] -= rate
        if dst is not None:
            Q[dst, src] += rate
        for c in counters:
            Q[lay.cnt[c], src] += rate
    Qc, Ql = parts.c, parts.l
    for src, rr in [(S, 1.0), (C, p['rr_reinfection_cleared']), (R, p['rr_reinfection_rec']), (TR, p['rr_reinfection_treat']), (M, p['rr_reinfection_inf'])]:
        flow(Ql, src, L, rr)
    flow(Qc, L, C, p['inf_cle']); flow(Qc, L, N0, p['inf_non']); flow(Qc, L, A, p['inf_asy'], ['first']); flow(Qc, L, M, p['k_latent'])
    flow(Qc, M, C, p['inf_cle']); flow(Qc, M, N0, p['inf_non_remote']); flow(Qc, M, A, p['inf_asy_remote'], ['first'])
    for src in [N0, N1]:
        flow(Qc, src, R, p['non_rec']); flow(Qc, src, A, p['non_asy'], ['first'] if src == N0 else [])
    flow(Qc, A, N1, p['asy_non']); flow(Qc, A, Y, p['asy_sym'])
    flow(Qc, Y, A, p['sym_asy'])
    flow(Qc, RL, Y, 1/p['dur_relapse'], ['first', 'relapse'])
    flow(parts.s, Y, lay.succ[0], 1.0, ['tx_start'])
    flow(parts.f, Y, lay.fail[0], 1.0, ['tx_start'])
    flow(parts.d, Y, None, 1.0, ['tb_death'])
    for j in lay.fail:
        flow(parts.d, j, None, 1.0, ['tb_death_tx'])
    if lay.kind == 'stages': # Erlang stages: k/D from each to the next; out of the last, cures to treated or relapsing, failures back to symptomatic TB
        rate = lay.k/p['tx_dur']
        for course in [lay.succ, lay.fail]:
            for a, b in zip(course[:-1], course[1:]):
                flow(Qc, a, b, rate)
        pr = p['tx_p_relapse']
        Qc[lay.succ[-1], lay.succ[-1]] -= rate
        Qc[TR, lay.succ[-1]] += (1 - pr)*rate
        Qc[RL, lay.succ[-1]] += pr*rate
        flow(Qc, lay.fail[-1], Y, rate)
    return parts


def value(p, key, t):
    return fo.value(p, key, t)


def schedules(trend):
    """ Exp 24's treatment and sym_dead over calendar time (tb_LAI_TPT's knots, linear between them): tx_rate, tx_success, sym_dead functions of the year """
    cov = lambda y: np.interp(y, [t_tx, 2010.0, 2024.0], [0.0, trend.tx_covg_2010, trend.tx_covg_2024], left=0.0)
    eff = lambda y: np.interp(y, [t_tx, 1995.0, 2024.0], [0.0, trend.efficacy_1995, trend.efficacy_2024], left=0.0)
    sd = lambda y: np.interp(y, [1995.0, 2024.0], [trend.sym_dead_1995, trend.sym_dead_2024])
    return dict(tx_rate=lambda y: float(-np.log1p(-cov(y))), tx_success=lambda y: float(eff(y)), sym_dead=lambda y: float(sd(y)), tx_p_relapse=trend.p_relapse)


#%% The starting state

def population(lows, country='India', year=1950):
    """ WPP's population [sex, group] on the grid `lows` (the last open), each of its groups spread evenly over its years, summing to 1 """
    pyr = ss.data.population(country, year=year)
    plows = np.array(ss.data.age_edges(pyr.index), dtype=float)
    pw = np.append(np.diff(plows), np.inf) # The last population group is open
    widths = np.append(np.diff(lows), np.inf)
    out = np.zeros((2, len(lows)))
    for sx, s in enumerate(sexes):
        vals = pyr[s].values.astype(float)
        for j, lo in enumerate(lows):
            hi = lo + widths[j]
            for q, (pl, w) in enumerate(zip(plows, pw)): # The overlap of the group with each population group, as a share of that group
                ov = min(hi, pl + w) - max(lo, pl)
                if ov > 0:
                    out[sx, j] += vals[q]*(ov/w if np.isfinite(w) else (1.0 if np.isinf(hi) else 0.0))
    return out/out.sum()


def initial(lay, lows, table, share_flagged, pop=None, country='India', year=1950):
    """
    The starting state [compartment, sex, group]: WPP's population on `lows`, split by `table` ({state: shares [sex, single year 0..100]},
    as exp 30's start.obj; every state but susceptible, which takes the rest; or {state: share} the same at every age). Non-infectious TB
    is split into NON1 (a share `share_flagged`, already asymptomatic this episode) and NON0, at every age
    """
    pop = population(lows, country, year) if pop is None else pop
    yr = np.clip(np.floor(lows + 1e-9).astype(int), 0, 100)
    X = np.zeros((lay.n, 2, len(lows)))
    named = {k: i for i, k in enumerate(base_states)}
    rest = np.ones((2, len(lows)))
    for k, v in table.items():
        if k == 'susceptible':
            continue
        sh = np.asarray(v, dtype=float)
        sh = sh[:, yr] if sh.ndim == 2 else np.full((2, len(lows)), float(sh))
        if not np.any(sh > 0):
            continue
        if k == 'non_infectious':
            X[N0] += (1 - share_flagged)*sh*pop
            X[N1] += share_flagged*sh*pop
        elif k in named:
            X[named[k]] += sh*pop
        else:
            raise KeyError(f'A starting share for {k}, which this start does not take')
        rest -= sh
    X[S] = rest*pop
    return X


#%% Outputs

def _stocks(X, lay, lows):
    """ What is recorded of a state: each compartment's total, and the living and active TB by exp 28's five groups and at 15+ """
    tot = X[:lay.ns].sum(axis=(1, 2))
    alive_g = X[:lay.ns].sum(axis=(0, 1))
    active_g = (X[A] + X[Y]).sum(axis=0)
    g5 = np.clip(np.digitize(lows + 1e-9, [g[0] for g in groups5.values()]) - 1, 0, len(groups5) - 1)
    a5 = np.bincount(g5, weights=alive_g, minlength=len(groups5))
    p5 = np.bincount(g5, weights=active_g, minlength=len(groups5))
    adult = lows >= 15 - 1e-9
    return tot, a5, p5, np.array([alive_g[adult].sum(), active_g[adult].sum()])


def per_step(year, tot, a5, p5, adult, flows, lay, dt):
    """ Each output per step, as exp 28's tb_per_step(): per person (per year for flows); stocks at index i are at t_i, flows the step from t_i """
    alive = tot.sum(axis=1)
    get = lambda i: tot[:, i]
    fl = lambda k: flows[:, counter_names.index(k)]
    out = sc.objdict(year=year)
    out.prevalence = (get(A) + get(Y))/alive
    out.prop_asymptomatic = get(A)/(get(A) + get(Y))
    out.incidence = fl('first')/alive/dt
    out.notifications = fl('tx_start')/alive/dt
    out.on_treatment = tot[:, np.concatenate([lay.succ, lay.fail])].sum(axis=1)/alive
    out.mortality = (fl('tb_death') + fl('tb_death_tx'))/alive/dt
    out.mortality_tx = fl('tb_death_tx')/alive/dt
    out.relapse = fl('relapse')/alive/dt
    out.recent, out.remote, out.cleared = get(L)/alive, get(M)/alive, get(C)/alive
    for k, g in enumerate(groups5):
        out[f'prev_{g}'] = p5[:, k]/a5[:, k]
        out[f'share_{g}'] = a5[:, k]/a5.sum(axis=1)
    out['prev_15+'] = adult[:, 1]/adult[:, 0]
    out.alive = alive/alive[0]
    out.tx_start = fl('tx_start')/alive # Per person per step, for targets()
    out.deaths = (fl('tb_death') + fl('tb_death_tx'))/alive
    out.first = fl('first')/alive
    out.population = alive # Shares of the starting population
    return out


stock_keys = ['prevalence', 'prop_asymptomatic', 'on_treatment', 'recent', 'remote', 'cleared', 'alive', 'prev_15+'] + [f'prev_{g}' for g in groups5] + [f'share_{g}' for g in groups5]
flow_keys = ['incidence', 'notifications', 'mortality', 'mortality_tx', 'relapse']
unscaled = ['prop_asymptomatic', 'alive'] + [f'share_{g}' for g in groups5]


def yearly(out, dt):
    """ Yearly means, as exp 28's yearly(): stocks at each step's end, flows over the step, dated by its middle; per 100,000 except shares and population """
    year = np.floor(out.year[1:] - dt/2 + 1e-9).astype(int)
    df = pd.DataFrame({k: (out[k][1:] if k in stock_keys else out[k][:-1])*(1.0 if k in unscaled else 1e5) for k in stock_keys + flow_keys})
    df['year'] = year
    return df.groupby('year').mean()


def calendar(out):
    """
    By calendar year: stocks on 1 January, and flows as totals over the year per 100,000 of the population on its 1 January (shares unscaled).
    Unlike yearly(), these don't depend on the step, on grids that put a step's start on each 1 January (quarter-year, monthly, 1/52 year)
    """
    year = out.year
    first = [int(np.argmin(np.abs(year - y))) for y in range(int(np.ceil(year[0] - 1e-9)), int(np.floor(year[-1] + 1e-9)) + 1)]
    yrs = np.round(year[first]).astype(int)
    df = pd.DataFrame({k: out[k][first]*(1.0 if k in unscaled else 1e5) for k in stock_keys}, index=pd.Index(yrs, name='year'))
    dt = year[1] - year[0]
    pop = out.population
    for k in flow_keys:
        counts = out[k]*pop*dt # People per step (shares of the starting population)
        df[k] = [counts[a:b].sum()/pop[a]*1e5 if b > a else np.nan for a, b in zip(first[:-1], first[1:])] + [np.nan]
    return df


def targets(out, over_detection):
    """ Exp 24's targets from a run's per-step outputs (point values on 1 January, totals over calendar years; as exp 27's target_values()),
    with prevalence also at 15+ (the survey's ages) """
    year = out.year
    at = lambda y: int(np.argmin(np.abs(year - y)))
    total = lambda k, y: out[k][(year >= y - 1e-9) & (year < y + 1 - 1e-9)].sum()
    prev = lambda y: out.prevalence[at(y)]
    inc = {y: total('first', y) for y in [2015, 2020, 2024]}
    deaths = {y: total('deaths', y) for y in [2015, 2024]}
    notif = (1 + over_detection)*total('tx_start', 2020)
    return dict(prevalence=prev(2020), prevalence_15plus=out['prev_15+'][at(2020)], prop_asymptomatic=out.prop_asymptomatic[at(2020)],
                prev_to_inc=prev(2020)/inc[2020], mort_to_inc=deaths[2024]/inc[2024], prev_to_notif=prev(2020)/notif, inc_ratio=inc[2024]/inc[2015],
                mort_ratio=deaths[2024]/deaths[2015], prev_2040_ratio=prev(2040)/prev(2020))


#%% Exact aging by cohort

def cohort_grid(h, top=100):
    """ Lower edges of groups one step h wide up to `top`, and the open group from there """
    n = int(np.ceil(top/h - 1e-9))
    return np.arange(n + 1)*h


def run_cohort(p, h, start, stop, table=None, share_flagged=0.0, X0=None, demog=None, foi='mid', keep_final=False):
    """
    Exact aging (exp 31's option 2): groups one step h wide, from calendar year `start` to `stop`.

    Args:
        p (dict): fastode.make_pars() parameters (mu = 0); tx_rate, tx_success, sym_dead may be functions of calendar year
        h (float): the step and the groups' width (years)
        table, share_flagged: the starting state (see initial()); or X0, a state on this grid
        demog (Demography): on this grid (default: WPP, not frozen)
        foi (str): 'mid' (predictor-corrector) or 'start' (the force of infection at the step's start)
        keep_final (bool): also return the state at the end

    Returns:
        objdict: per-step outputs (per_step()), with lows, h, and optionally X (the final state)
    """
    lows = cohort_grid(h)
    ng = len(lows)
    demog = Demography(lows, cohort=True) if demog is None else demog
    assert demog.cohort and len(demog.lows) == ng, 'run_cohort() needs a Demography(cohort=True) on its grid'
    K = max(1, int(round(p['tx_dur']/h)))
    lay = Layout(('slots', K))
    X = initial(lay, lows, table, share_flagged) if X0 is None else X0.copy()
    parts = generator_parts(p, lay)
    nsteps = int(round((stop - start)/h))
    year = start + np.arange(nsteps + 1)*h
    tot = np.zeros((nsteps + 1, lay.ns)); a5 = np.zeros((nsteps + 1, 5)); p5 = np.zeros_like(a5); adult = np.zeros((nsteps + 1, 2))
    flows = np.zeros((nsteps + 1, len(counter_names)))
    people = np.arange(lay.ns)
    cnt = np.array([lay.cnt[c] for c in counter_names])
    beta, kappa, pr = p['beta'], p['trans_asymp'], p['tx_p_relapse']
    sK, fK = lay.succ[-1], lay.fail[-1]
    Xf = X.reshape(lay.n, 2*ng) # A view: [compartment, sex × group]

    def release(frac):
        """ A share frac of the conveyors' last slots leaves: cures to treated or relapsing, failures back to symptomatic TB """
        s, f = frac*Xf[sK], frac*Xf[fK]
        Xf[TR] += (1 - pr)*s; Xf[RL] += pr*s; Xf[Y] += f
        Xf[sK] -= s; Xf[fK] -= f
        return

    def force(v, extra=0.0):
        return beta*(kappa*v[A] + v[Y])/(v[:lay.ns].sum() + extra)

    for i in range(nsteps):
        tot[i], a5[i], p5[i], adult[i] = _stocks(X, lay, lows)
        tm = year[i] + h/2
        th, sc_, sd = value(p, 'tx_rate', tm), value(p, 'tx_success', tm), value(p, 'sym_dead', tm)
        Qb = parts.c + th*sc_*parts.s + th*(1 - sc_)*parts.f + sd*parts.d
        mx = demog.mx(tm) # [sex, group]
        surv = np.exp(-mx*h).ravel()
        mh = (mx*h).ravel()
        cfac = np.where(mh > 1e-12, -np.expm1(-mh)/np.where(mh > 1e-12, mh, 1.0), 1.0) # Counters: flows over the step thinned by deaths during it
        women0 = Xf[people][:, :ng].sum(axis=0)

        release(0.5)
        lam0 = force(Xf.sum(axis=1))
        if foi == 'mid':
            P0 = expm((Qb + lam0*parts.l)*h)
            born_pred = h*np.sum(demog.fert(tm)*women0) # The step's births, all susceptible, in the living at its end
            lam1 = force(P0 @ (Xf @ surv), born_pred) # The end of the step, predicted: survival commutes with Q, which acts within a group
            lam = 0.5*(lam0 + lam1)
        else:
            lam = lam0
        P = expm((Qb + lam*parts.l)*h)
        Xf[:] = P @ Xf
        Xf[people] *= surv[None, :]
        Xf[cnt] *= cfac[None, :]
        release(1.0)
        flows[i] = Xf[cnt].sum(axis=1)
        Xf[cnt] = 0.0

        # The conveyors move along one slot
        for course in [lay.succ, lay.fail]:
            Xf[course[1:]] = Xf[course[:-1]]
            Xf[course[0]] = 0.0

        # Births over the step (women at its start and end), surviving half a step in the first group; then everyone ages one group
        women1 = Xf[people][:, :ng].sum(axis=0)
        births = h*np.sum(demog.fert(tm)*0.5*(women0 + women1))
        X[:, :, -1] += X[:, :, -2]
        X[:, :, 1:-1] = X[:, :, :-2].copy()
        X[:, :, 0] = 0.0
        m0 = demog.mx_first(tm)
        born = expm((Qb + lam*parts.l)*h/2)[:lay.ns, S] # Born on average half a step before its end: half a step of its TB transitions from susceptible
        X[:lay.ns, 0, 0] = demog.p_female*births*np.exp(-m0[0]*h/2)*born
        X[:lay.ns, 1, 0] = (1 - demog.p_female)*births*np.exp(-m0[1]*h/2)*born
    tot[-1], a5[-1], p5[-1], adult[-1] = _stocks(X, lay, lows)
    out = per_step(year, tot, a5, p5, adult, flows, lay, h)
    out.lows, out.h, out.layout = lows, h, lay
    if keep_final:
        out.X = X.copy()
    return out


#%% Aging at a rate: a copy of v4's ODE

def run_rate(p, width, start, stop, t_eval_step, table=None, share_flagged=0.0, stages=12, rtol=1e-8, atol=1e-13):
    """
    Constant-rate aging (exp 31's option 1), as v4's ODE: groups of `width` years up to 100 (100+ open), aged out of at rate 1/width,
    the course as `stages` Erlang stages, solved by solve_ivp (RK45) in continuous time. Outputs at every `t_eval_step` years from `start`
    (e.g. the step of the Starsim run compared with), with flows integrated exactly between them
    """
    lows = np.round(np.arange(0, 100 + 1e-9, width), 6)
    ng = len(lows)
    demog = Demography(lows)
    lay = Layout(('stages', stages))
    X0 = initial(lay, lows, table, share_flagged)[:lay.ns]
    parts = generator_parts(p, lay)
    ns = lay.ns
    cnt = [lay.cnt[c] for c in counter_names]
    beta, kappa = p['beta'], p['trans_asymp']
    age_rate = np.full(ng, 1/width); age_rate[-1] = 0.0

    def rhs(t, z):
        X = z[:ns*2*ng].reshape(ns, 2, ng)
        tot = X.sum(axis=(1, 2))
        lam = beta*(kappa*tot[A] + tot[Y])/tot.sum()
        th, s, sd = value(p, 'tx_rate', t), value(p, 'tx_success', t), value(p, 'sym_dead', t)
        Q = parts.c + lam*parts.l + th*s*parts.s + th*(1 - s)*parts.f + sd*parts.d
        Xf = X.reshape(ns, 2*ng)
        QX = Q[:, :ns] @ Xf
        dX = QX[:ns].reshape(ns, 2, ng) - demog.mx(t)[None]*X
        out_age = X*age_rate
        dX -= out_age
        dX[:, :, 1:] += out_age[:, :, :-1]
        births = np.sum(demog.fert(t)*X[:, 0].sum(axis=0))
        dX[S, 0, 0] += demog.p_female*births
        dX[S, 1, 0] += (1 - demog.p_female)*births
        return np.concatenate([dX.ravel(), QX[cnt].sum(axis=1)])

    n = int(round((stop - start)/t_eval_step))
    year = start + np.arange(n + 1)*t_eval_step
    z0 = np.concatenate([X0.ravel(), np.zeros(len(cnt))])
    sol = solve_ivp(rhs, (start, year[-1]), z0, t_eval=year, method='RK45', rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(sol.message)
    tot = np.zeros((len(year), ns)); a5 = np.zeros((len(year), 5)); p5 = np.zeros_like(a5); adult = np.zeros((len(year), 2))
    for i in range(len(year)):
        X = np.zeros((lay.n, 2, ng)); X[:ns] = sol.y[:ns*2*ng, i].reshape(ns, 2, ng)
        tot[i], a5[i], p5[i], adult[i] = _stocks(X, lay, lows)
    cum = sol.y[ns*2*ng:].T
    flows = np.vstack([np.diff(cum, axis=0), np.zeros((1, len(cnt)))]) # Index i: the step from t_i
    out = per_step(year, tot, a5, p5, adult, flows, lay, t_eval_step)
    out.lows, out.layout, out.nfev = lows, lay, sol.nfev
    return out


#%% The settled start with ages (exp 30's, in this model)

def settle(p, h, init, share_flagged, years=200, start=1950, foi='mid', drift_years=20):
    """
    Exp 30's starting state in this model: from shares the same at every age (`init`, e.g. exp 27's untreated equilibrium), run `years`
    with `start`'s demography frozen, no treatment and sym_dead at p's value at `start` (pass its 1995 value), by exact aging at step h.
    Returns the final state, and the drift (the largest change in any group's share of any state over the last `drift_years`)
    """
    q = dict(p); q.update(tx_rate=0.0, sym_dead=value(p, 'sym_dead', start))
    lows = cohort_grid(h)
    demog = Demography(lows, frozen=start, cohort=True)
    first = run_cohort(q, h, start, start + years - drift_years, table=init, share_flagged=share_flagged, demog=demog, foi=foi, keep_final=True)
    last = run_cohort(q, h, start + years - drift_years, start + years, X0=first.X, demog=demog, foi=foi, keep_final=True)
    sh = lambda X: X[:last.layout.ns]/np.maximum(X[:last.layout.ns].sum(axis=0), 1e-300)
    drift = float(np.max(np.abs(sh(last.X) - sh(first.X))))
    return sc.objdict(X=last.X, lows=lows, layout=last.layout, drift=drift, years=years)


def start_table(st):
    """ A settled run's shares by single year of age and sex (weighted by its people), as exp 30's start.obj table, and the share of
    non-infectious TB already asymptomatic this episode """
    X, lows, lay = st.X, st.lows, st.layout
    yr = np.clip(np.floor(lows + 1e-9).astype(int), 0, 100)
    pop = X[:lay.ns].sum(axis=0)
    den = np.zeros((2, 101))
    for sx in range(2):
        np.add.at(den[sx], yr, pop[sx])
    def avg(rows):
        num = np.zeros((2, 101))
        v = X[rows].sum(axis=0) if np.ndim(rows) else X[rows]
        for sx in range(2):
            np.add.at(num[sx], yr, v[sx])
        return np.divide(num, den, out=np.zeros_like(num), where=den > 0)
    table = {k: avg(i) for i, k in enumerate(base_states) if k not in ['non0', 'non1']}
    table['non_infectious'] = avg([N0, N1])
    nontot = X[[N0, N1]].sum()
    return table, float(X[N1].sum()/nontot) if nontot > 0 else 0.0
