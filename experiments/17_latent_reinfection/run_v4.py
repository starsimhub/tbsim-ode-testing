"""
Experiment 17, v4 environment. See README.md.

Question 1: v4 agents vs v4's ODE for the two-state latent model (remote progression = phi x recent) with latent reinfection on,
for set B at phi = 0.003, 0.01 and 0.03. Each arm's ODE beta is calibrated to 500/100k active prevalence at equilibrium; the agents
start at the ODE's equilibrium shares and run 1900-2000 at monthly and weekly dt (40 reps each, new seeds), reduced to yearly means.

Question 1, population size (q1n): phi = 0.01, monthly dt, at 6,250-200,000 agents, to test whether the small prevalence
shortfall at 50,000 agents is a finite-population effect (which would shrink like 1/N).

Question 2 (q2): the two-state model as v4 agents in a cohort at a constant force of infection, for tb_LAI_TPT's 12 face-validity
sets at each phi, with the per-person outputs of cohort.py; and the same cohorts as exact chains (tbode4.latent), two-state and
tb_LAI_TPT's model solved exactly (continuous time).

Builds on Cliff's v4 helpers (v4/common.py, v4/run_equilibrium.py), used unchanged. Run from the repo root:

    v4/.venv/bin/python experiments/17_latent_reinfection/run_v4.py q1    # or q1n (needs q1 first), q2, summary, summary_n

Rerunning resumes: calibrations, ODE runs and agent arms already saved are kept.
"""
import sys
import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss

outdir = sc.thispath(__file__)
v4dir = outdir.parent.parent / 'v4'
sys.path.insert(0, str(v4dir))
import common as c # noqa: E402 (Cliff's v4 helpers)
import run_equilibrium as req # noqa: E402
import tbode4 # noqa: E402
import tbsim # noqa: E402
sys.path.insert(0, str(outdir))
import cohort as co # noqa: E402

phis = [0.003, 0.01, 0.03]
dts = sc.objdict(month=ss.months(1), week=ss.weeks(1))
reps = 40
seed0 = 1000 # New seeds: v4's draws are keyed by seed, so seeds used in exp 13 would repeat its runs
start, stop = 1900, 2000
n_agents = 50_000
pi = c.parsets.B['rr_reinfection_rec'] # tb_LAI_TPT's reinfection risk ratio, the same for everyone who can be reinfected
fits_file = outdir.parent / '16_two_state_vs_lai_tpt' / 'fits.csv'
fn_q1 = outdir / 'q1.obj'
fn_q1n = outdir / 'q1_n.obj'
pop_sizes = {6250: 160, 12_500: 160, 25_000: 80, 200_000: 40} # Agents: reps, for q1n; 50,000 is q1's arm
seed0_n = 5000 # q1n's seeds: 5000 + 1000*i + rep, for the i-th population size

# Question 2: cohorts
samples_csv = outdir.parents[2] / 'tb_LAI_TPT' / 'results/hm_real_20260728/trajectory_selection/collapsed_samples_eps1.2.csv'
n_sets = 12              # tb_LAI_TPT's face-validity parameter sets: the first 12 rows
n_cohort = 250_000
seeds_q2 = [0, 1, 2, 3]
years_q2 = 30
lams = [0.0, 0.10, 0.60] # Force of infection, per year
fn_q2 = outdir / 'q2_v4.obj'
cascade_keys = ['inf_cle', 'non_rec', 'non_asy', 'asy_non', 'asy_sym', 'sym_asy']
events_v4 = dict( # Each cohort.py event, as v4 transitions (those that exist in the model)
    new      = ['foi_cle', 'foi_rec'],
    reset    = ['foi_remote'],
    progress = ['latent_non', 'latent_asy', 'remote_non', 'remote_asy'],
    clear    = ['clearance', 'clearance_remote'],
    asy      = ['latent_asy', 'remote_asy', 'worsening', 'remission'],
    sym      = ['symptoms'],
)

# Each state's share of everyone alive, as recorded per step
shares_of = sc.objdict(
    recent         = ['infection'],
    remote         = ['infection_remote'],
    latent         = ['infection', 'infection_remote'],
    non_infectious = ['non_infectious'],
    asymptomatic   = ['asymptomatic'],
    symptomatic    = ['symptomatic'],
    prevalence     = ['asymptomatic', 'symptomatic'],
    cleared        = ['cleared'],
    recovered      = ['recovered'],
)


def two_state_pars(phi, sample=0):
    """ The two-state progression parameters for one parameter set and phi: exp 16's arm-1 fit, plus latent reinfection at pi """
    fits = pd.read_csv(fits_file)
    row = fits[(fits['sample'] == sample) & (fits.arm == 1) & np.isclose(fits.phi, phi)]
    assert len(row) == 1, f'No unique exp 16 fit for sample {sample}, phi {phi}'
    row = row.iloc[0]
    pars = {k: float(row[k]) for k in ['inf_non', 'inf_asy', 'inf_non_remote', 'inf_asy_remote', 'k_latent']}
    assert np.isclose(row.inf_cle, c.parsets.B['inf_cle']) # Arm 1 keeps the set's own cascade
    pars['rr_reinfection_inf'] = pi
    pars['fit'] = dict(z_Y1=float(row.z_Y1), z_Ratio=float(row.z_Ratio))
    return pars


def outputs(sim):
    """ Shares of everyone alive by state, and TB deaths per person per year """
    r = sim.results.tb
    alive = sim.results.alive.values
    out = sc.objdict(year=np.array(sim.t.yearvec))
    for key, states in shares_of.items():
        out[key] = sum(r[s].values for s in states)/alive
    out.tb_deaths = r.tb_death.values/alive/sim.t.dt_year
    return out


def share(o, key):
    """ One output's yearly series; 'susceptible' is derived as everyone not in another recorded state (no one is in treatment) """
    if key == 'susceptible':
        return 1 - sum(o[k] for k in ['latent', 'non_infectious', 'prevalence', 'cleared', 'recovered'])
    return o[key]


def run_one(phi, dt, seed, beta, init, solver='abm', n=n_agents):
    """ One run of set B's two-state model at phi, from the ODE's equilibrium, reduced to yearly means """
    pars = {k: v for k, v in two_state_pars(phi).items() if k != 'fit'}
    tb = c.make_tb('B', beta, init=init, model='tb', **pars)
    sim = c.make_sim(tb, solver=solver, dt=dts[dt], start=start, stop=stop, n_agents=n, seed=seed)
    T = sc.timer()
    sim.run()
    out = req.yearly(outputs(sim))
    out.time = T.toc(output=True)
    return out


def ratio(runs, ode, key, years=(1950, 2000)):
    """ Agents ÷ ODE for one output: each run's mean over the years, divided by the ODE's """
    sel = lambda o: (o.year >= years[0]) & (o.year < years[1])
    ref = share(ode, key)[sel(ode)].mean()
    return np.array([share(o, key)[sel(o)].mean() for o in runs])/ref


def q1():
    data = sc.load(fn_q1) if fn_q1.exists() else sc.objdict(calib=sc.objdict(), ode=sc.objdict(), abm=sc.objdict())
    for phi in phis:
        key = f'phi{phi}'
        if key not in data.calib:
            pars = two_state_pars(phi)
            T = sc.timer()
            cal = c.calibrate('B', model='tb', **{k: v for k, v in pars.items() if k != 'fit'})
            data.calib[key] = sc.objdict(phi=phi, pars=pars, beta=cal.beta, equilibrium=cal.equilibrium, prevalence=cal.prevalence)
            sc.save(fn_q1, data)
            print(f'{key}: beta = {cal.beta:.5f}, prevalence = {1e5*cal.prevalence:.2f}/100k, fit z = {pars["fit"]} ({T.toc(output=True):.0f} s)', flush=True)
        cal = data.calib[key]
        for dt in dts:
            akey = f'{key}_{dt}'
            if akey not in data.ode:
                data.ode[akey] = run_one(phi, dt, 0, cal.beta, cal.equilibrium, solver='ode')
            if akey in data.abm and len(data.abm[akey]) >= reps:
                continue
            T = sc.timer()
            data.abm[akey] = sc.parallelize(run_one, iterkwargs=[dict(phi=phi, dt=dt, seed=seed0 + s, beta=cal.beta, init=cal.equilibrium) for s in range(reps)], ncpus=tbode4.ncpus)
            sc.save(fn_q1, data)
            r = ratio(data.abm[akey], data.ode[akey], 'prevalence')
            print(f'{akey:14s} {reps} reps in {T.toc(output=True):.0f} s; prevalence 1950–2000, agents ÷ ODE: {r.mean():.3f} ± {r.std(ddof=1)/np.sqrt(reps):.3f}', flush=True)
    tbode4.save_versions(outdir, 'v4')
    (outdir / 'versions.json').replace(outdir / 'versions_v4.json') # Kept apart from question 2's v0.10 versions
    summarise(data)
    return data


def summarise(data):
    """ Agents ÷ ODE (mean ± SE, 10–90% range across reps) for each arm and output, over 1950–2000 """
    rows = []
    for akey, runs in data.abm.items():
        for key in ['prevalence', 'latent', 'recent', 'remote', 'susceptible', 'cleared', 'recovered', 'non_infectious', 'tb_deaths']:
            rows.append(dict(arm=akey, output=key, **stats(ratio(runs, data.ode[akey], key))))
    df = pd.DataFrame(rows)
    df.to_csv(outdir / 'q1_summary.csv', index=False)
    with pd.option_context('display.width', 200, 'display.max_rows', 200):
        print(df.round(4).to_string(index=False))
    pooled = np.mean([ratio(runs, data.ode[akey], 'prevalence') for akey, runs in data.abm.items()], axis=0) # Each seed averaged across the arms
    print(f'Prevalence, pooled over arms by seed: {pooled.mean():.4f} ± {pooled.std(ddof=1)/np.sqrt(len(pooled)):.4f}')
    return df


def stats(r):
    """ Mean, SE and 10–90% range of agents ÷ ODE across reps """
    return dict(reps=len(r), mean=r.mean(), se=r.std(ddof=1)/np.sqrt(len(r)), p10=np.percentile(r, 10), p90=np.percentile(r, 90))


def q1n():
    """ Population size: phi = 0.01, monthly dt, at each size in pop_sizes, against the same ODE run as q1's 50,000-agent arm """
    q1data = sc.load(fn_q1)
    cal = q1data.calib['phi0.01']
    data = sc.load(fn_q1n) if fn_q1n.exists() else sc.objdict()
    for i, (n, nreps) in enumerate(pop_sizes.items()):
        key = f'N{n}'
        if key in data and len(data[key]) >= nreps:
            continue
        T = sc.timer()
        seeds = [seed0_n + 1000*i + s for s in range(nreps)]
        data[key] = sc.parallelize(run_one, iterkwargs=[dict(phi=0.01, dt='month', seed=seed, beta=cal.beta, init=cal.equilibrium, n=n) for seed in seeds], ncpus=tbode4.ncpus)
        sc.save(fn_q1n, data)
        print(f'{key}: {nreps} reps in {T.toc(output=True):.0f} s', flush=True)
    tbode4.save_versions(outdir, 'v4')
    (outdir / 'versions.json').replace(outdir / 'versions_v4.json')
    summarise_n(data, q1data)
    return data


def summarise_n(data, q1data):
    """ Agents ÷ ODE by population size (q1's 50,000-agent arm included), and the new runs pooled by inverse variance """
    ode = q1data.ode['phi0.01_month']
    arms = dict(data)
    arms['N50000'] = q1data.abm['phi0.01_month']
    arms = dict(sorted(arms.items(), key=lambda kv: int(kv[0][1:])))
    rows = []
    for key, runs in arms.items():
        for out in ['prevalence', 'susceptible', 'latent', 'cleared', 'recovered', 'tb_deaths']:
            rows.append(dict(n_agents=int(key[1:]), output=out, **stats(ratio(runs, ode, out))))
    df = pd.DataFrame(rows)
    df.to_csv(outdir / 'q1_n_summary.csv', index=False)
    with pd.option_context('display.width', 200, 'display.max_rows', 200):
        print(df.round(4).to_string(index=False))
    new = df[(df.output == 'prevalence') & (df.n_agents != n_agents)]
    w = 1/new.se**2
    print(f'Prevalence, the new population sizes pooled (inverse variance): {(w*new["mean"]).sum()/w.sum():.4f} ± {1/np.sqrt(w.sum()):.4f}')
    return df


class TBCohort(tbsim.TB):
    """
    The two-state model in a cohort: everyone in recent latent infection at the start, no transmission, births or deaths (sym_dead = 0),
    and infection at a constant force of infection `lam` (per year) times the risk ratio `pi`: of cleared and recovered people (back to
    recent latent infection), and, with `latent_reinf`, of remote latent infections (back to recent). Declared as transitions, so it
    runs under either solver.
    """
    def __init__(self, lam=0.0, latent_reinf=True, pi=1.0, **kwargs):
        super().__init__(beta=0.0, mu=0.0, sym_dead=0.0, init=dict(infection=1.0), **kwargs)
        if lam > 0:
            rate = ss.peryear(pi*lam)
            self.define_transitions(foi_cle=self.cleared >> self.infection @ rate, foi_rec=self.recovered >> self.infection @ rate)
            if latent_reinf:
                self.define_transitions(foi_remote=self.infection_remote >> self.infection @ rate)
        return


class CohortTracker(ss.Analyzer):
    """ Feeds each step's transitions to cohort.Tracker, and records the recent and remote shares (index i is time i·dt) """

    def init_post(self):
        super().init_post()
        self.tr = co.Tracker(len(self.sim.people), self.sim.t.npts, self.sim.t.dt_year, offset=0)
        self.recent, self.remote = np.zeros(self.sim.t.npts), np.zeros(self.sim.t.npts)
        return

    def step(self):
        tb = self.sim.diseases.tb
        ev = {}
        for key, names in events_v4.items():
            us = [np.asarray(getattr(tb, name).uids) for name in names if hasattr(tb, name)]
            ev[key] = np.concatenate(us) if us else []
        self.tr.step(self.sim.ti, ev)
        n = len(self.sim.people)
        self.recent[self.sim.ti] = len(tb.infection.uids)/n
        self.remote[self.sim.ti] = len(tb.infection_remote.uids)/n
        return


def cohort_pars(sample, phi):
    """ One set's two-state parameters at phi (exp 16's arm-1 fit, with the set's own cascade), and its pi; rates per year """
    fits = pd.read_csv(fits_file)
    row = fits[(fits['sample'] == sample) & (fits.arm == 1) & np.isclose(fits.phi, phi)]
    assert len(row) == 1, f'No unique exp 16 fit for sample {sample}, phi {phi}'
    row = row.iloc[0]
    post = pd.read_csv(samples_csv).iloc[sample]
    p = {k: float(post[k]) for k in cascade_keys}
    p.update({k: float(row[k]) for k in ['inf_non', 'inf_asy', 'inf_non_remote', 'inf_asy_remote', 'k_latent']})
    assert np.isclose(row.inf_cle, p['inf_cle']) and np.isclose(row.non_rec, p['non_rec']) # Arm 1 keeps the set's own cascade
    return p, float(post['rr_reinfection_rec'])


def run_cohort(sample, phi, lam, on, seed, n=n_cohort):
    """ One v4 agent cohort: the per-person outputs, with the recent and remote shares """
    p, pi = cohort_pars(sample, phi)
    tb = TBCohort(lam=lam, latent_reinf=on, pi=pi, **p)
    sim = ss.Sim(tb, solver='abm', n_agents=n, total_pop=n, start=2000, stop=2000 + years_q2, dt=ss.months(1), rand_seed=seed, verbose=0,
                 analyzers=CohortTracker())
    T = sc.timer()
    sim.run()
    an = sim.analyzers[0]
    out = an.tr.results()
    out.recent, out.remote = an.recent, an.remote
    out.update(sample=sample, phi=phi, lam=lam, on=on, seed=seed, time=T.toc(output=True))
    return out


def chain_curves(p, structure, lam, on, per_year=12):
    """ The exact cohort curves from one person's chain: A, B, C (ever progressed, ever asymptomatic, ever symptomatic) and latent shares """
    kw = dict(lam=lam, latent_reinf=on, cascade=True, tmax=years_q2, per_year=per_year)
    out = sc.objdict()
    t, X, idx = tbode4.cohort_curve(p, structure, **kw)
    out.t = t
    out.latent = X[:, idx.latent].sum(axis=1)
    if structure == 'two-state':
        out.recent, out.remote = X[:, idx.recent], X[:, idx.remote]
    _, X, idx = tbode4.cohort_curve(p, structure, absorbing=['NON', 'ASY'], **kw)
    out.A = X[:, idx.NON] + X[:, idx.ASY]
    _, X, idx = tbode4.cohort_curve(p, structure, absorbing=['ASY'], **kw)
    out.B = X[:, idx.ASY]
    _, X, idx = tbode4.cohort_curve(p, structure, absorbing=['SYM'], **kw)
    out.C = X[:, idx.SYM]
    return out


def q2_chains():
    """ Exact chains for every set, lam and arm: two-state at each phi, and tb_LAI_TPT's model solved exactly (k_asy = k_non = k) """
    post = pd.read_csv(samples_csv)
    chains = sc.objdict()
    for s in range(n_sets):
        row = post.iloc[s]
        exact = {k: float(row[k]) for k in cascade_keys}
        exact.update(inf_non=float(row.inf_non), inf_asy=float(row.inf_asy), k=float(row.k_asy), pi=float(row.rr_reinfection_rec), sym_dead=0.0)
        for lam, on in q2_arms():
            chains[f'p{s}_exact_lam{lam}_{onoff(on)}'] = chain_curves(exact, 'exact', lam, on)
            for phi in phis:
                p, pi = cohort_pars(s, phi)
                chains[f'p{s}_phi{phi}_lam{lam}_{onoff(on)}'] = chain_curves(dict(p, pi=pi, sym_dead=0.0), 'two-state', lam, on)
    return chains


def q2_arms():
    return [(lam, on) for lam in lams for on in ([True] if lam == 0 else [True, False])]


def onoff(on):
    return 'on' if on else 'off'


def q2():
    """ All v4 agent cohorts and chains; saved as it goes, and resumes """
    data = sc.load(fn_q2) if fn_q2.exists() else sc.objdict(runs=sc.objdict())
    post = pd.read_csv(samples_csv).head(n_sets)
    data.settings = dict(n_sets=n_sets, n_agents=n_cohort, seeds=seeds_q2, years=years_q2, lams=lams, phis=phis, dt='monthly')
    data.weights = post.sample_weight.values/post.sample_weight.sum()
    if 'chains' not in data:
        T = sc.timer()
        data.chains = q2_chains()
        sc.save(fn_q2, data)
        print(f'{len(data.chains)} chains in {T.toc(output=True):.0f} s', flush=True)
    jobs = [dict(sample=s, phi=phi, lam=lam, on=on, seed=seed) for lam, on in q2_arms() for phi in phis for s in range(n_sets) for seed in seeds_q2]
    key = lambda j: f"p{j['sample']}_phi{j['phi']}_lam{j['lam']}_{onoff(j['on'])}_s{j['seed']}"
    todo = [j for j in jobs if key(j) not in data.runs]
    for chunk in [todo[i:i + 48] for i in range(0, len(todo), 48)]:
        T = sc.timer()
        res = sc.parallelize(run_cohort, iterkwargs=chunk, ncpus=tbode4.ncpus)
        for j, r in zip(chunk, res):
            data.runs[key(j)] = r
        sc.save(fn_q2, data)
        print(f'{len(data.runs)}/{len(jobs)} runs ({T.toc(output=True):.0f} s for {len(chunk)})', flush=True)
    tbode4.save_versions(outdir, 'v4')
    (outdir / 'versions.json').replace(outdir / 'versions_v4.json')
    return data


def q2_perperson():
    """
    Seed 0 of each set's reinfection-on cohorts again, for each person's total time latent and number of latent spells (added to
    cohort.py after the main runs). v4's draws are keyed by seed, so these are the same runs as in q2_v4.obj: checked on curve B.
    """
    fn_pp = outdir / 'q2_v4_perperson.obj'
    saved = sc.load(fn_q2)
    jobs = [dict(sample=s, phi=phi, lam=lam, on=True, seed=0) for lam in lams for phi in phis for s in range(n_sets)]
    T = sc.timer()
    res = sc.parallelize(run_cohort, iterkwargs=jobs, ncpus=tbode4.ncpus)
    out = sc.objdict()
    for j, r in zip(jobs, res):
        key = f"p{j['sample']}_phi{j['phi']}_lam{j['lam']}_on_s0"
        assert np.array_equal(r.B, saved.runs[key].B), f'{key} does not reproduce the saved run'
        out[key] = sc.objdict(time_latent=r.hists.time_latent, time_latent_open=r.hists.time_latent_open, n_spells=r.n_spells, n=r.n, dt_year=r.dt_year)
    sc.save(fn_pp, out)
    print(f'{len(jobs)} runs in {T.toc(output=True):.0f} s; all reproduce the saved runs (curve B identical)', flush=True)
    return out


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'q1'
    if which == 'q2p':
        q2_perperson()
        sys.exit()
    if which == 'q1':
        q1()
    elif which == 'q1n':
        q1n()
    elif which == 'q2':
        q2()
    elif which == 'summary':
        summarise(sc.load(fn_q1))
    elif which == 'summary_n':
        summarise_n(sc.load(fn_q1n), sc.load(fn_q1))
    else:
        raise ValueError(f'Unknown part "{which}"')
