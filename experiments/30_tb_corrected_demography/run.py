"""
Experiment 30, v4 environment: do v4's agents match its ODE with TB, treatment and the corrected demography (CONVENTIONS.md 10)? See README.md.
Run from the repo root:

    PYTHONPATH=v4 v4/.venv/bin/python experiments/30_tb_corrected_demography/run.py check    # quick test of the starting state in both solvers
    PYTHONPATH=v4 v4/.venv/bin/python experiments/30_tb_corrected_demography/run.py timing   # one ODE run with TB at each age-group width and step
    PYTHONPATH=v4 v4/.venv/bin/python experiments/30_tb_corrected_demography/run.py start    # the starting states with ages (frozen 1950 demography)
    PYTHONPATH=v4 v4/.venv/bin/python experiments/30_tb_corrected_demography/run.py tb       # the ODE at each width, then the agents (after start)
    PYTHONPATH=v4 v4/.venv/bin/python experiments/30_tb_corrected_demography/run.py agents   # only the agents (after start and tb)
    PYTHONPATH=v4 v4/.venv/bin/python experiments/30_tb_corrected_demography/run.py more     # 40 more agent reps for set 6 (new seeds), merged into the saved runs
    PYTHONPATH=v4 v4/.venv/bin/python experiments/30_tb_corrected_demography/run.py summary  # tables and figures from the saved runs

Environment variables: EXP30_NCPUS (cores for the agent runs; default tbode4.ncpus), EXP30_WIDTHS (the ODE's widths for tb, e.g. 'single,quarter'),
EXP30_QUICK=1 (a small, short test of every step, written to quick/).
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


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


r27 = load_module('run27', root / 'experiments' / '27_calibrated_treatment' / 'run.py') # Exp 27's sets and untreated equilibrium without ages
r28 = load_module('run28', root / 'experiments' / '28_age_demography' / 'run.py') # Exp 28's TB outputs (per step, yearly) and helpers by age

# Settings (README)
quick = os.environ.get('EXP30_QUICK') == '1'
ncpus = int(os.environ.get('EXP30_NCPUS', tbode4.ncpus))
outdir = here / 'quick' if quick else here
outdir.mkdir(exist_ok=True)
start, stop = 1950, 2041
burn_years = 200 # The frozen run for the starting state
n_agents, reps = 100_000, 20
sets_used = [0, 6]
dts = sc.objdict(month=ss.months(1), week=ss.weeks(1))
widths = sc.objdict(single=1.0, quarter=0.25, eighth=0.125) # The ODE's age groups (years), each up to 100+
start_width = 'eighth' # The frozen run for the starting state uses the finest width; its shares are averaged to single years of age, the table both solvers start from
agent_edges = tbode4.age_edges(1.0) # Agents' groups (results only): their starting table needs groups no coarser than its bins (a finer table gives wrong shares)
ode_widths = os.environ.get('EXP30_WIDTHS', ','.join(widths)).split(',')
if quick:
    stop, burn_years, n_agents, reps = 1998, 30, 20_000, 2
    sets_used = [0]
tol = 1e-3 # Added to agents' ages when counting them by age group (their 32-bit ages; CONVENTIONS.md 10)
groups5 = r28.groups5
states_tb = r28.states_tb
fn = sc.objdict(timing=outdir / 'timing.csv', start=outdir / 'start.obj', ode=outdir / 'ode.obj', runs=outdir / 'runs.obj')
wpp = tbode4.WPP()
sets = sc.load(r27.fn.sets)
ss.solvers.Engine.group_bytes = 4e9 # Fine age groups keep their results at every step (Starsim keeps them only every few steps above 100 MB per run)


def _drawn_start(init_post):
    """ Agents' starting states, drawn at random within each group (multinomial, seeded by the run's seed) instead of Starsim's exact counts by group,
    round(share × agents in the group): with single years of age by sex, a rare state (e.g. 0.3% with active TB) rounds to nobody in small old-age
    groups, so 65+ started with half the ODE's prevalence. Swaps the rounding only while agents' states are set at the start """
    def wrapped(self):
        Comp = ss.solvers.compartmental.Compartmental
        orig = Comp.__dict__['round_preserving']
        rng = np.random.default_rng([int(self.sim.pars.rand_seed), 30])
        def drawn(y):
            y = np.asarray(y, dtype=float)
            n = int(round(y.sum()))
            return rng.multinomial(n, y/y.sum()).astype(float) if n > 0 else np.zeros_like(y)
        Comp.round_preserving = staticmethod(drawn)
        try:
            return init_post(self)
        finally:
            Comp.round_preserving = orig
    return wrapped


ss.solvers.agents.Agents.init_post = _drawn_start(ss.solvers.agents.Agents.init_post)


#%% The model

def tb_model(s, dt, init, treat=True):
    """ Exp 27's TBTreat for set s with mu = 0 (demography from outside). Without treatment: none, and sym_dead held at its 1995 value """
    t = s.trend
    nat = {k: v for k, v in s.nat.items() if k not in ['sym_dead', 'mu']}
    if not treat:
        return tbode4.TBTreat(init=init, course='fixed', sym_dead=float(t.sym_dead_1995), mu=0, **nat)
    tx = tbode4.lai_tpt_treatment(row=dict(tx_covg_2010=t.tx_covg_2010, tx_covg_2024=t.tx_covg_2024, efficacy_1995=t.efficacy_1995,
                                           efficacy_2024=t.efficacy_2024, p_relapse=t.p_relapse), step=dts[dt].years, end=stop)
    sym_dead = ss.peryear(ss.timeseries({1995.0: t.sym_dead_1995, 2024.0: t.sym_dead_2024}, interpolation='linear'))
    return tbode4.TBTreat(init=init, course='fixed', tx_rate=tx.tx_rate, tx_success=tx.tx_success, tx_p_relapse=t.p_relapse, dur_relapse_tx=1.5,
                          sym_dead=sym_dead, mu=0, **nat)


class AgeCounts(ss.Analyzer):
    """ Agents alive and with active TB (asymptomatic + symptomatic) in each of exp 28's five age groups, counted with `tol` added to ages.
    Index i is the state at t_i, as Starsim's results of states """
    def init_post(self):
        super().init_post()
        self.lows = np.array([g[0] for g in groups5.values()], dtype=float)
        self.alive = np.zeros((self.sim.t.npts, len(self.lows)))
        self.active = np.zeros_like(self.alive)
        self._count(0)
        return

    def _count(self, i):
        ppl, tb = self.sim.people, self.sim.diseases.tb
        if i >= len(self.alive):
            return
        g = lambda u: np.clip(np.digitize(np.asarray(ppl.age.raw[u], dtype=float) + tol, self.lows) - 1, 0, len(self.lows) - 1)
        self.alive[i] = np.bincount(g(np.asarray(ppl.auids)), minlength=len(self.lows))
        act = np.concatenate([np.asarray(tb.asymptomatic.uids), np.asarray(tb.symptomatic.uids)]).astype(int)
        self.active[i] = np.bincount(g(act), minlength=len(self.lows))
        return

    def step(self): # After this step's transitions: the state at t_(ti+1)
        self._count(self.ti + 1)
        return


def make_sim(i, dt, solver, width='quarter', start_state=None, treat=True, frozen=False, years=None, rand_seed=0):
    """ One sim with the corrected demography. start_state: from starting_state() (None: exp 27's equilibrium without ages, as exp 28) """
    s = sets[i]
    abm = solver == 'abm'
    init0, non0, flag0 = r27.starting_state(s)
    if start_state is None:
        init, share_flagged = init0, flag0
    else:
        init, share_flagged = start_state.init, start_state.share_flagged
    tb = tb_model(s, dt, init, treat=treat)
    kw = dict(people=ss.population('India', year=start, n_agents=n_agents if abm else None), total_pop=n_agents if abm else 1.0)
    if abm:
        kw.update(rand_seed=rand_seed, analyzers=[tbode4.Incidence(share_flagged=share_flagged), tbode4.TxChecks(), AgeCounts()])
        demog = tbode4.agent_demographics(wpp, year=start)
    else:
        demog = tbode4.ode_demographics(wpp, frozen=frozen)
    t0, t1 = (start, stop) if years is None else years
    sim = ss.Sim(tb, demographics=demog, start=t0, stop=t1, dt=dts[dt], solver=solver, verbose=0, **kw)
    if abm:
        sim.stratify(age=agent_edges, sex=tbode4.sexes) # Agents age exactly; their groups only report results, and read the starting table
    else:
        e = tbode4.age_edges(widths[width])
        sim.stratify(age=e, sex=tbode4.joint_sex(e, year=start))
    non0_start = non0 if start_state is None else start_state.non0.get(width) # The ODE's; agents use share_flagged
    return sim, non0_start, share_flagged


#%% The ODE's first-entry incidence with demography (exp 28's, also returning NON0 over time)

def ode_incidence(sim, non0, frozen=False):
    """ First entries to asymptomatic TB per step (plus relapses), and NON0 (non-infectious, not yet asymptomatic this episode) over time """
    r, tb = sim.results.tb, sim.diseases.tb
    get = lambda k: np.asarray(r[k].values, dtype=float) if k in r else np.zeros(len(r['symptomatic'].values))
    rate = lambda name: float(getattr(tb, name).rate.raw.value)
    non_rec, non_asy = rate('recovery'), rate('worsening')
    yv = np.asarray(sim.t.yearvec, dtype=float)
    non = np.stack([r28.by_age_sex(r.non_infectious, s)[0] for s in tbode4.sexes], axis=1) # [time, sex, age group]
    lows = r28.by_age_sex(r.non_infectious, tbode4.sexes[0])[1]
    idx = np.clip(np.digitize(lows, tbode4.wpp_edges) - 1, 0, len(tbode4.wpp_edges) - 1)
    mxt = np.array([wpp.mx_at(start if frozen else y)[:, idx] for y in yv])
    tot = non.sum(axis=(1, 2))
    mu = np.divide((non*mxt).sum(axis=(1, 2)), tot, out=np.zeros(len(yv)), where=tot > 0)
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
    return get('latent_asy') + get('remote_asy') + first_non + get('relapse'), n0, tot


#%% The starting state with ages

def shares_by_group(sim, idx=-1):
    """ Each TB state's share of each age and sex group at time index idx: {state: [sex, age]}, and the groups' lower edges and people """
    r = sim.results.tb
    vals = {k: np.stack([r28.by_age_sex(r[k], s)[0][idx] for s in tbode4.sexes]) for k in states_tb if k in r}
    pop = sum(vals.values())
    lows = r28.by_age_sex(r['susceptible'], tbode4.sexes[0])[1]
    return {k: v/pop for k, v in vals.items()}, lows, pop


def init_table(shares, lows):
    """ init= for TBTreat: each state but susceptible (which takes the rest; agents refuse a value by group for every state) as ss.by(age, sex) """
    out = {}
    for k, v in shares.items():
        if k == 'susceptible' or not np.any(v > 0):
            continue
        out[k] = ss.by(pd.DataFrame(dict(age=np.tile(lows, 2), sex=np.repeat(tbode4.sexes, len(lows)), value=np.concatenate([v[0], v[1]]))))
    return out


def coarsen(shares, lows, pop, edges):
    """ Shares averaged over the fine groups in each coarser group (weighted by the frozen run's people) """
    j = np.clip(np.digitize(lows + 1e-9, edges) - 1, 0, len(edges) - 1)
    out = {}
    for k, v in shares.items():
        num = np.zeros((2, len(edges))); den = np.zeros((2, len(edges)))
        for sx in range(2):
            np.add.at(num[sx], j, v[sx]*pop[sx])
            np.add.at(den[sx], j, pop[sx])
        out[k] = np.divide(num, den, out=np.zeros_like(num), where=den > 0)
    return out


def starting_state(i):
    """ The ODE's untreated steady state at 1950's demography (frozen), by age and sex, from exp 27's state without ages; burn_years at monthly steps """
    T = sc.timer()
    sim, non0, _ = make_sim(i, 'month', 'ode', width=start_width, treat=False, frozen=True, years=(start, start + burn_years))
    sim.run()
    _, n0, nontot = ode_incidence(sim, non0, frozen=True)
    sh, lows, pop = shares_by_group(sim)
    yv = np.asarray(sim.t.yearvec, dtype=float)
    j20 = int(np.argmin(np.abs(yv - (yv[-1] - 20))))
    sh20, _, _ = shares_by_group(sim, j20)
    drift = {k: float(np.max(np.abs(sh[k] - sh20[k]))) for k in sh} # Largest change in a group's share over the last 20 years
    unflagged = n0[-1]/nontot[-1] # NON0 ÷ NON at the end: the share of non-infectious people not yet asymptomatic this episode
    e = np.array(agent_edges)
    out = sc.objdict(set=i, lows=lows, shares=sh, people=pop, drift=drift, share_flagged=float(1 - unflagged), edges=e, table=coarsen(sh, lows, pop, e),
                     non0=sc.objdict(), secs=T.toc(output=True)) # non0 is set from each width's starting population (set_non0)
    out.init = init_table(out.table, e) # The starting table of both solvers, by single year of age and sex
    return out


def set_non0(st):
    """ NON0 at the start of each solver: the unflagged share × the starting share in non-infectious TB (the ODE: of a population of 1) """
    for w in widths:
        sim, _, _ = make_sim(st.set, 'month', 'ode', width=w, start_state=st, years=(start, start + 1))
        sim.run()
        non = float(np.asarray(sim.results.tb.non_infectious.values, dtype=float)[0])
        st.non0[w] = non*(1 - st.share_flagged)
    return st


def run_start():
    T = sc.timer()
    data = sc.objdict()
    for i in sets_used:
        st = starting_state(i)
        data[f'set{i}'] = set_non0(st)
        top = sorted(st.drift.items(), key=lambda kv: -kv[1])[:3]
        print(f'set {i}: frozen run {st.secs:.0f} s; largest change in a group\'s share over its last 20 years: {top}; share flagged {st.share_flagged:.4f}', flush=True)
    sc.save(fn.start, data)
    tbode4.save_versions(outdir, 'v4')
    print(f'Starting states done ({T.toc(output=True):.0f} s)', flush=True)
    return data


#%% Checks and timing

def group5_shares(sim, abm):
    """ At the start: each TB state's share of each of the five age groups (agents: counted from their ages with tol; the ODE: from results) """
    lows5 = np.array([g[0] for g in groups5.values()], dtype=float)
    if abm:
        ppl, tb = sim.people, sim.diseases.tb
        g = lambda u: np.clip(np.digitize(np.asarray(ppl.age.raw[u], dtype=float) + tol, lows5) - 1, 0, len(lows5) - 1)
        alive = np.bincount(g(np.asarray(ppl.auids)), minlength=5)
        return {k: np.bincount(g(np.asarray(getattr(tb, k).uids)), minlength=5)/alive for k in states_tb if hasattr(tb, k)}
    sh, lows, pop = shares_by_group(sim, 0)
    return {k: r28.regroup((v*pop).sum(0)[None, :], lows, lows5)[0]/r28.regroup(pop.sum(0)[None, :], lows, lows5)[0] for k, v in sh.items()}


def run_check():
    """ Quick test (use with EXP30_QUICK=1): the starting state from a short frozen run; then each solver's shares at the start against the table """
    T = sc.timer()
    i = sets_used[0]
    st = set_non0(starting_state(i))
    print(f'Frozen run ({burn_years} years, {start_width} ages): {st.secs:.0f} s; share flagged {st.share_flagged:.4f}; NON0 at the start {dict(st.non0)}')
    rows, ref = [], None
    for solver, w in [('ode', 'eighth'), ('ode', 'quarter'), ('ode', 'single'), ('abm', None)]:
        sim, _, _ = make_sim(i, 'month', solver, width=w or 'quarter', start_state=st, years=(start, start + 1), rand_seed=30001)
        if solver == 'abm':
            sim.init()
            got = group5_shares(sim, True) # Agents at the start, before the first step
        sim.run()
        if solver == 'ode':
            got = group5_shares(sim, False) # Index 0: the start
        ref = got if ref is None else ref # The eighth-year ODE's start, the reference for the others
        if solver == 'ode': # Each fine group's shares at t0 against the single-year table holding it
            sh, lows, _ = shares_by_group(sim, 0)
            j = np.clip(np.digitize(lows + 1e-9, st.edges) - 1, 0, len(st.edges) - 1)
            err = max(float(np.max(np.abs(sh[k] - st.table[k][:, j]))) for k in st.table if k in sh)
            print(f'ODE, {w} ages: largest difference from the table in any group and state at the start: {err:.2e}')
        for k in ['susceptible', 'infection', 'infection_remote', 'cleared', 'non_infectious', 'asymptomatic', 'symptomatic']:
            rows.append(dict(solver=f'{solver} {w or ""}'.strip(), state=k, **{g: got[k][j]/ref[k][j] if ref[k][j] > 0 else np.nan for j, g in enumerate(groups5)}))
        if solver == 'abm':
            a = sim.analyzers[2]
            sx = [r28.by_age_sex(sim.results.tb[k])[0] for k in states_tb if k in sim.results.tb]
            alive_res = r28.to_groups5(sum(x for x in sx), r28.by_age_sex(sim.results.tb['susceptible'])[1])
            print(f'Agents: AgeCounts ÷ Starsim results by age, alive, first and last index: {np.round(a.alive[0]/alive_res[0], 4)}, {np.round(a.alive[-1]/alive_res[-1], 4)}')
            fem = np.asarray(sim.results.pregnancy.birth.select(sex='female').values).sum()/max(np.asarray(sim.results.pregnancy.birth.values).sum(), 1)
            print(f'Agents: female share at the start {float(np.asarray(sim.people.female.raw[sim.people.auids]).mean()):.4f}; births female in the first year {fem:.4f}')
    df = pd.DataFrame(rows)
    with pd.option_context('display.width', 200, 'display.max_columns', 20):
        print(f'Share of each age group in each state at the start, ÷ the eighth-year ODE\'s; agents: one run of {n_agents}')
        print(df.round(3).to_string(index=False))
    df.to_csv(outdir / 'check.csv', index=False)
    print(f'check finished in {T.toc(output=True):.0f} s')
    return df


def run_timing():
    """ One ODE run with TB and treatment at each width and step, over the full period, from exp 27's state without ages """
    rows = []
    for dt in dts:
        for w in widths:
            T = sc.timer()
            sim, _, _ = make_sim(sets_used[0], dt, 'ode', width=w)
            sim.run()
            rows.append(dict(dt=dt, width=w, secs=T.toc(output=True)))
            print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(fn.timing, index=False)
    return rows


#%% The main runs

def seed_of(i, dt, rep):
    return 30500 + 100*(2*sets_used.index(i) + list(dts).index(dt)) + rep


def run_agents(i, dt, rep, st):
    T = sc.timer()
    sim, _, _ = make_sim(i, dt, 'abm', start_state=st, rand_seed=seed_of(i, dt, rep))
    sim.run()
    inc, chk, ac = sim.analyzers[0], sim.analyzers[1], sim.analyzers[2]
    out = r28.tb_per_step(sim, inc)
    for k, g in enumerate(groups5): # By age from our own counts (ages + tol), not Starsim's results by age group
        out[f'prev_{g}'] = ac.active[:, k]/ac.alive[:, k]
        out[f'share_{g}'] = ac.alive[:, k]/ac.alive.sum(axis=1)
    return sc.objdict(set=i, dt=dt, rep=rep, seed=seed_of(i, dt, rep), yearly=r28.yearly(out, sim.t.dt_year), checks=dict(chk.c), pop_check=out.pop_check,
                      secs=T.toc(output=True))


def run_ode(i, dt, w, st):
    T = sc.timer()
    sim, non0, _ = make_sim(i, dt, 'ode', width=w, start_state=st)
    sim.run()
    inc, _, _ = ode_incidence(sim, non0)
    out = r28.tb_per_step(sim, inc)
    return sc.objdict(set=i, dt=dt, width=w, yearly=r28.yearly(out, sim.t.dt_year), pop_check=out.pop_check, secs=T.toc(output=True))


def run_tb(agents_only=False):
    """ The ODE at each width (unless agents_only), then the agents """
    T = sc.timer()
    starts = sc.load(fn.start)
    if not agents_only:
        ode = sc.objdict()
        for i in sets_used:
            for dt in dts:
                for w in ode_widths:
                    ode[f'set{i}_{dt}_{w}'] = run_ode(i, dt, w, starts[f'set{i}'])
                    print(f'ODE set {i} {dt} {w}: {ode[f"set{i}_{dt}_{w}"].secs:.0f} s', flush=True)
        sc.save(fn.ode, ode)
    jobs = [dict(i=i, dt=dt, rep=r, st=starts[f'set{i}']) for i in sets_used for dt in dts for r in range(reps)]
    runs = sc.parallelize(run_agents, iterkwargs=jobs, ncpus=ncpus)
    sc.save(fn.runs, runs)
    tbode4.save_versions(outdir, 'v4')
    print(f'tb finished in {T.toc(output=True):.0f} s; agent runs {np.mean([o.secs for o in runs]):.0f} s each on average', flush=True)
    return runs


def run_more(sets_more=(6,), extra=40):
    """ More agent reps with new seeds (reps 20 onwards) for sets_more at both steps, merged into the saved runs """
    T = sc.timer()
    starts = sc.load(fn.start)
    runs = sc.load(fn.runs)
    jobs = []
    for i in sets_more:
        for dt in dts:
            done = max(r.rep for r in runs if r.set == i and r.dt == dt) + 1
            jobs += [dict(i=i, dt=dt, rep=r, st=starts[f'set{i}']) for r in range(done, done + extra)]
    new = sc.parallelize(run_agents, iterkwargs=jobs, ncpus=ncpus)
    sc.save(fn.runs, runs + new)
    tbode4.save_versions(outdir, 'v4')
    print(f'more finished in {T.toc(output=True):.0f} s: {len(new)} runs added; agent runs {np.mean([o.secs for o in new]):.0f} s each on average', flush=True)
    return new


#%% Summaries and figures

colors = {(0, 'month'): 'tab:blue', (0, 'week'): 'tab:cyan', (6, 'month'): 'tab:red', (6, 'week'): 'tab:orange'}


def ode_refs(ode):
    """ The ODE runs as exp 28's by_window() takes them, for each reference: each width, and exact aging extrapolated (2 × eighth − quarter) """
    refs = sc.objdict()
    for ref in ['exact'] + list(widths):
        out = []
        for i in sets_used:
            for dt in dts:
                if ref == 'exact':
                    y = 2*ode[f'set{i}_{dt}_eighth'].yearly - ode[f'set{i}_{dt}_quarter'].yearly
                else:
                    y = ode[f'set{i}_{dt}_{ref}'].yearly
                out.append(sc.objdict(set=i, dt=dt, yearly=y))
        refs[ref] = out
    return refs


def table(bw, outputs, cols=None):
    t = bw[bw.output.isin(outputs)].copy()
    t['value'] = t.apply(lambda x: f'{x.ratio:.3f} ± {x.se:.3f}' if x.se == x.se else f'{x.ratio:.3f}', axis=1)
    return t.pivot_table(index='output', columns=cols or ['set', 'dt', 'window'], values='value', aggfunc='first').reindex(outputs)


def summarize():
    runs, old, ode, starts = sc.load(fn.runs), sc.load(outdir / 'runs_rounded.obj'), sc.load(fn.ode), sc.load(fn.start)
    refs = ode_refs(ode)
    bw = r28.by_window(runs, refs.exact)
    bw.to_csv(outdir / 'by_window.csv', index=False)
    bw1 = r28.by_window(runs, refs.single)
    bw1.to_csv(outdir / 'by_window_single.csv', index=False)
    bwr = r28.by_window(old, refs.exact)
    bwr.to_csv(outdir / 'by_window_rounded.csv', index=False)
    rows = [] # The ODE's widths against exact aging (deterministic)
    for w in widths:
        for o, e in zip(refs[w], refs.exact):
            for win in r28.windows:
                ow, ew = r28.window_mean(o.yearly, win), r28.window_mean(e.yearly, win)
                rows.append(dict(set=o.set, dt=o.dt, window=win, width=w, **{k: ow[k]/ew[k] for k in ['alive', 'prevalence'] + [f'prev_{g}' for g in groups5]}))
    widths_df = pd.DataFrame(rows)
    widths_df.to_csv(outdir / 'ode_widths.csv', index=False)
    chk = pd.DataFrame([dict(set=r.set, dt=r.dt, **r.checks) for r in runs]).groupby(['set', 'dt']).sum()
    chk.to_csv(outdir / 'pitfall_checks.csv')
    main = ['prevalence', 'prop_asymptomatic', 'incidence_per_case', 'notifications_per_case', 'on_treatment_per_case', 'mortality_per_case',
            'recent', 'remote', 'cleared', 'alive'] + [f'prev_{g}' for g in groups5] + [f'share_{g}' for g in groups5]
    age = ['prevalence', 'alive'] + [f'prev_{g}' for g in groups5]
    b28 = {k: pd.read_csv(root / 'experiments' / '28_age_demography' / f) for k, f in [('5-year', 'by_window.csv'), ('single', 'by_window_single.csv')]}
    first = [] # 65+ prevalence in 1950, agents ÷ the exactly aged ODE, drawn start and Starsim's rounded one
    for name, rr in [('drawn', runs), ('rounded', old)]:
        for e in refs.exact:
            v = np.array([r.yearly.loc[start, 'prev_65+'] for r in rr if r.set == e.set and r.dt == e.dt])/e.yearly.loc[start, 'prev_65+']
            first.append(dict(start=name, set=e.set, dt=e.dt, ratio=v.mean(), se=v.std(ddof=1)/np.sqrt(len(v))))
    first = pd.DataFrame(first)
    with open(outdir / 'summary.log', 'w') as fh, pd.option_context('display.width', 300, 'display.max_columns', 40):
        nrep = ', '.join(f'set {i} {dt}: {sum(1 for r in runs if r.set == i and r.dt == dt)}' for i in sets_used for dt in dts)
        txt = f'Agents ÷ ODE (1.00 = match), mean ± SE over reps ({nrep}) of window means of yearly values; per_case = each ÷ prevalence. Windows: ' + \
              ', '.join(f'{k} {a}–{b - 1}' for k, (a, b) in r28.windows.items()) + '\n\n'
        txt += 'Agents ÷ the ODE with exact aging (extrapolated: 2 × eighth-year − quarter-year)\n' + table(bw, main).to_string() + '\n\n'
        txt += 'Agents ÷ the ODE with single years of age\n' + table(bw1, age).to_string() + '\n\n'
        txt += "The ODE at each width ÷ the ODE with exact aging (window means)\n" + widths_df.round(4).to_string(index=False) + '\n\n'
        txt += "Starsim's rounded start (each state's count in each age and sex group rounded to whole agents) against the drawn start:\n"
        txt += '65+ prevalence in 1950, agents ÷ exact ODE:\n' + first.round(3).to_string(index=False) + '\n'
        txt += 'Agents with the rounded start ÷ exact ODE, by window\n' + table(bwr, age).to_string() + '\n\n'
        for k, b in b28.items():
            txt += f"Exp 28 (agents with exp 28's demography and start ÷ its ODE with {k} age groups), for comparison\n" + table(b, age).to_string() + '\n\n'
        txt += 'Starting states: largest change in an age and sex group\'s share of any state over the frozen run\'s last 20 years: ' + \
               ', '.join(f'set {st.set} {max(st.drift.values()):.1e}' for st in starts.values()) + '\n'
        txt += 'ODE run times (s): ' + ', '.join(f'{k} {o.secs:.0f}' for k, o in ode.items()) + '\n'
        txt += f'Agent run times: {np.mean([r.secs for r in runs if r.dt == "month"]):.0f} s monthly, {np.mean([r.secs for r in runs if r.dt == "week"]):.0f} s weekly (zebra, 80 in parallel)\n\n'
        txt += 'Pitfall checks (summed over reps):\n' + chk.to_string() + '\n'
        fh.write(txt)
        print(txt)
    plot(runs, old, refs, bw, b28)
    return bw


def plot(runs, old, refs, bw, b28):
    years = np.arange(start, stop)
    # Burden and prevalence by age over time, for each set: agents' reps and mean, the exactly aged ODE and the single-year one
    for i in sets_used:
        for name, items, ncol in [('burden', [('prevalence', 'Prevalence'), ('incidence', 'Incidence per year'), ('mortality', 'TB mortality per year'),
                                              ('notifications', 'Notifications per year')], 4),
                                  ('prev_by_age', [(f'prev_{g}', f'Prevalence, aged {g}') for g in groups5], 5)]:
            fig, axs = plt.subplots(len(dts), ncol, figsize=(4*ncol, 3.4*len(dts)), squeeze=False)
            for a, dt in enumerate(dts):
                rs = [r for r in runs if r.set == i and r.dt == dt]
                ex = [e for e in refs.exact if e.set == i and e.dt == dt][0].yearly
                o1 = [e for e in refs.single if e.set == i and e.dt == dt][0].yearly
                for b, (k, label) in enumerate(items):
                    ax = axs[a, b]
                    for r in rs:
                        ax.plot(r.yearly.index, r.yearly[k], color='tab:blue', alpha=0.12, lw=0.6)
                    m = pd.concat([r.yearly[k] for r in rs], axis=1).mean(axis=1)
                    ax.plot(m.index, m.values, color='tab:blue', lw=1.6, label=f'Agents, mean of {len(rs)}')
                    ax.plot(ex.index, ex[k], 'k--', lw=1.2, label='ODE, exact aging')
                    ax.plot(o1.index, o1[k], ':', color='0.5', lw=1.2, label='ODE, single years')
                    ax.set_title(f'{label} per 100,000 ({dt})', fontsize=9)
                    if a == 0 and b == 0:
                        ax.legend(fontsize=7)
            fig.suptitle(f'Set {i}: agents and the ODE, with the corrected demography'); fig.tight_layout()
            fig.savefig(outdir / f'tb_{name}_set{i}.png', dpi=110); plt.close(fig)
    # Agents ÷ exact ODE over time: prevalence, 65+ prevalence and population, with the drawn and rounded starts
    items = [('prevalence', 'Prevalence, all ages'), ('prev_0-14', 'Prevalence, 0–14'), ('prev_65+', 'Prevalence, 65+'), ('alive', 'Population')]
    fig, axs = plt.subplots(2, len(items), figsize=(4.5*len(items), 7), squeeze=False)
    for a, (label_r, rr) in enumerate([('drawn start', runs), ("Starsim's rounded start", old)]):
        for b, (k, label) in enumerate(items):
            ax = axs[a, b]
            for e in refs.exact:
                v = np.array([r.yearly[k].loc[years].values for r in rr if r.set == e.set and r.dt == e.dt])/e.yearly[k].loc[years].values
                m, se = v.mean(0), v.std(0, ddof=1)/np.sqrt(len(v))
                c = colors[(e.set, e.dt)]
                ax.plot(years, m, color=c, lw=1.2, label=f'set {e.set}, {e.dt}')
                ax.fill_between(years, m - 2*se, m + 2*se, color=c, alpha=0.15)
            ax.axhline(1, color='k', lw=0.6); ax.axhspan(0.99, 1.01, color='0.9', zorder=0)
            ax.set_ylim(0.85, 1.15) if k != 'alive' else ax.set_ylim(0.98, 1.02)
            ax.set_title(f'{label}: agents ÷ ODE ({label_r})', fontsize=9)
            if a == 0 and b == 0:
                ax.legend(fontsize=7)
    fig.suptitle('Agents (mean ± 2 SE over reps: set 0, 20; set 6, 60 with the drawn start, 20 with the rounded one) ÷ the ODE with exact aging, by year; grey band ±1%'); fig.tight_layout()
    fig.savefig(outdir / 'ratio_over_time.png', dpi=110); plt.close(fig)
    # By window: this experiment against exp 28's
    show = ['prevalence', 'prop_asymptomatic', 'incidence_per_case', 'notifications_per_case', 'mortality_per_case', 'alive'] + [f'prev_{g}' for g in groups5]
    fig, axs = plt.subplots(len(r28.windows), 1, figsize=(16, 3.4*len(r28.windows)), squeeze=False)
    for a, win in enumerate(r28.windows):
        ax = axs[a, 0]
        for c, (i, dt) in enumerate([(i, dt) for i in sets_used for dt in dts]):
            x = np.arange(len(show)) + (c - 1.5)*0.12
            d = bw[(bw.set == i) & (bw.dt == dt) & (bw.window == win)].set_index('output').reindex(show)
            ax.errorbar(x, d.ratio, yerr=2*d.se, fmt='o', ms=3.5, color=colors[(i, dt)], label=f'exp 30: set {i}, {dt}')
            d28 = b28['single'][(b28['single'].set == i) & (b28['single'].dt == dt) & (b28['single'].window == win)].set_index('output').reindex(show)
            ax.plot(x + 0.05, d28.ratio, 'x', ms=4, color=colors[(i, dt)], alpha=0.6, label='exp 28 (single-year ODE)' if c == 0 else None)
        ax.axhline(1, color='k', lw=0.6); ax.axhspan(0.99, 1.01, color='0.9', zorder=0)
        ax.set_xticks(np.arange(len(show))); ax.set_xticklabels(show, fontsize=8); ax.set_ylim(0.88, 1.08)
        ax.set_title(f'{win}: {r28.windows[win][0]}–{r28.windows[win][1] - 1}', fontsize=9)
        if a == 0:
            ax.legend(fontsize=7, ncol=5)
    fig.suptitle('Agents ÷ ODE by window (exp 30: exact aging, mean ± 2 SE over reps, set 0 20, set 6 60; exp 28: its own demography and start, single-year ODE)')
    fig.tight_layout(); fig.savefig(outdir / 'by_window.png', dpi=110); plt.close(fig)
    return


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'check'
    dict(check=run_check, timing=run_timing, start=run_start, tb=run_tb, agents=lambda: run_tb(agents_only=True), summary=summarize, more=run_more)[cmd]()
