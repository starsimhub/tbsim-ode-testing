"""
Experiment 30, v4 environment: do v4's agents match its ODE with TB, treatment and the corrected demography (CONVENTIONS.md 10)? See README.md.
Run from the repo root:

    PYTHONPATH=v4 v4/.venv/bin/python experiments/30_tb_corrected_demography/run.py check    # quick test of the starting state in both solvers
    PYTHONPATH=v4 v4/.venv/bin/python experiments/30_tb_corrected_demography/run.py timing   # one ODE run with TB at each age-group width and step
    PYTHONPATH=v4 v4/.venv/bin/python experiments/30_tb_corrected_demography/run.py start    # the starting states with ages (frozen 1950 demography)
    PYTHONPATH=v4 v4/.venv/bin/python experiments/30_tb_corrected_demography/run.py tb       # the ODE at each width, then the agents (after start)

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


def run_tb():
    T = sc.timer()
    starts = sc.load(fn.start)
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


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'check'
    dict(check=run_check, timing=run_timing, start=run_start, tb=run_tb)[cmd]()
