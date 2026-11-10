"""
Experiment 25: acceptance test of Starsim's network fix (starsim 6dcad9fd, in d83f4f87) with tbsim 619c296. See README.md.

Questions 1-2 repeat exp 18's q3, q1, q2 (network workflow); question 3 repeats exp 22's diagnose.py; question 4 repeats exp 22's
conversion with kappa unchanged and converted; question 5 repeats exp 21's homogeneous monthly arm. The code is copied from those
experiments, not imported, so that their outputs stay where they are (parallel workers re-import modules). Run from the repo root:

    PYTHONPATH=v4 experiments/25_starsim_network_fix/.venv/bin/python experiments/25_starsim_network_fix/run.py q3           # calibration on the network ODE
    PYTHONPATH=v4 experiments/25_starsim_network_fix/.venv/bin/python experiments/25_starsim_network_fix/run.py q1           # per transition
    PYTHONPATH=v4 experiments/25_starsim_network_fix/.venv/bin/python experiments/25_starsim_network_fix/run.py q2           # equilibrium
    PYTHONPATH=v4 experiments/25_starsim_network_fix/.venv/bin/python experiments/25_starsim_network_fix/run.py diagnose     # per step
    PYTHONPATH=v4 experiments/25_starsim_network_fix/.venv/bin/python experiments/25_starsim_network_fix/run.py conversions  # converted sets, predictions
    PYTHONPATH=v4 experiments/25_starsim_network_fix/.venv/bin/python experiments/25_starsim_network_fix/run.py convert      # agents for question 4
    PYTHONPATH=v4 experiments/25_starsim_network_fix/.venv/bin/python experiments/25_starsim_network_fix/run.py homog        # question 5
    PYTHONPATH=v4 experiments/25_starsim_network_fix/.venv/bin/python experiments/25_starsim_network_fix/run.py summary      # every summary from saved runs

Rerunning resumes: parts already saved are kept.
"""
import os
os.environ.setdefault('NUMBA_NUM_THREADS', '1') # Before numba is imported (by starsim), so every process agrees (tbode4 sets it too, but after starsim here)
import sys
import importlib.util
import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss
from scipy.optimize import brentq

outdir = sc.thispath(__file__)
root = outdir.parent.parent
sys.path.insert(0, str(root / 'v4'))
sys.path.insert(0, str(root / 'experiments' / '17_latent_reinfection'))
import common as c # noqa: E402 (Cliff's v4 helpers)
import run_equilibrium as req # noqa: E402
import tbode4 # noqa: E402
from tbode4 import fastode as fo # noqa: E402
import run_v4 as rv17 # noqa: E402 (exp 17: the two-state model's parameters)


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


exp21 = load_module('exp21_run', root / 'experiments' / '21_calibrated_comparison' / 'run.py') # Sets, homogeneous ODE runs, outputs (read only)
exp19 = exp21.exp19

K = 10 # Mean degree of ss.RandomNet (fresh contacts every step)
dts = sc.objdict(month=ss.months(1), week=ss.weeks(1))
target = 0.005 # Active (asymptomatic + symptomatic) prevalence, for q3
n_agents = 50_000
seed0 = 25000
fn = sc.objdict(q1=outdir / 'q1.obj', q2=outdir / 'q2.obj', q3=outdir / 'q3.obj', diag=outdir / 'diagnose.obj', conv=outdir / 'conversions.obj',
                runs=outdir / 'runs.obj', homog=outdir / 'homog.obj')
contact_flows = ['transmission', 'reinfect_cle', 'reinfect_rec', 'reinfect_trd', 'reinfect_remote']


# %% Questions 1-2: exp 18's sets on Starsim's network ODE

sets18 = ['A', 'B2'] # A: tbsim's defaults; B2: exp 17's two-state set B at phi = 0.01, every reinfection risk at pi
start18, stop18, reps18 = 1900, 2000, 40


def set_kw(name):
    """ The parameter set and extra tbsim.TB parameters of each named set (exp 18) """
    if name == 'A':
        return 'A', {}
    return 'B', {k: v for k, v in rv17.two_state_pars(0.01).items() if k != 'fit'}


def rr_of(name, flow):
    """ Each transition by contact's risk relative to a susceptible's """
    import tbsim
    parset, pars = set_kw(name)
    p = sc.mergedicts(c.parsets[parset], pars)
    dflt = tbsim.tb_pars()
    get = lambda k: float(p[k]) if k in p else float(dflt[k].raw if hasattr(dflt[k], 'raw') else dflt[k])
    return dict(transmission=1.0, reinfect_cle=get('rr_reinfection_cleared'), reinfect_rec=get('rr_reinfection_rec'),
                reinfect_trd=get('rr_reinfection_treat'), reinfect_remote=get('rr_reinfection_inf'))[flow]


def make_sim18(name, beta, solver, dt, network, init, n=None, seed=0, start=start18, stop=stop18):
    """ A sim of one of exp 18's sets: `beta` is per year (mass action) or per contact (network) """
    parset, pars = set_kw(name)
    tb = c.make_tb(parset, beta, init=init, model='tb', **pars)
    kw = dict(start=start, stop=stop, dt=dt, total_pop=n or 1e5, rand_seed=seed, verbose=0)
    if network:
        kw['networks'] = ss.RandomNet(n_contacts=ss.constant(K), dur=0)
    if n:
        kw['n_agents'] = n
    return ss.Sim(tb, solver=solver, **kw)


def ode_equilibrium(name, beta, dt, network, init, years=600, tol=1e-8, maxiter=20):
    """ The ODE's equilibrium shares at this step, run from `init` until they settle """
    shares = init
    for i in range(maxiter):
        sim = make_sim18(name, beta, 'ode', dt, network, shares, start=0, stop=years)
        sim.run()
        new = c.end_shares(sim)
        change = max(abs(new[k] - shares.get(k, 0)) for k in new)
        shares = new
        if change < tol:
            break
    return shares


def q3():
    """ The per-contact b for 500/100k on the network ODE at each step, against the matched b from the mass-action beta """
    data = sc.load(fn.q3) if fn.q3.exists() else sc.objdict()
    for name in sets18:
        parset, pars = set_kw(name)
        if f'{name}_mass' not in data:
            T = sc.timer()
            cal = c.calibrate(parset, model='tb', **pars)
            data[f'{name}_mass'] = sc.objdict(beta=cal.beta, equilibrium=cal.equilibrium, prevalence=cal.prevalence)
            sc.save(fn.q3, data)
            print(f'{name} mass action: beta = {cal.beta:.5f} ({T.toc(output=True):.0f} s)', flush=True)
        ma = data[f'{name}_mass']
        for dkey, dt in dts.items():
            key = f'{name}_net_{dkey}'
            if key in data:
                continue
            T = sc.timer()
            matched = tbode4.beta_edge_from_ode(ma.beta, K, dt.years)
            eq = {}
            def gap(b):
                eq[b] = ode_equilibrium(name, b, dt, True, ma.equilibrium)
                return c.prevalence(eq[b]) - target
            b = brentq(gap, 0.9*matched, 1.1*matched, xtol=1e-7*matched)
            shares = eq[b] if b in eq else ode_equilibrium(name, b, dt, True, ma.equilibrium)
            data[key] = sc.objdict(beta_edge=b, matched=matched, equilibrium=shares, prevalence=c.prevalence(shares))
            sc.save(fn.q3, data)
            print(f'{key}: per-contact b = {b:.5f} (matched {matched:.5f}, ratio {b/matched:.5f}) ({T.toc(output=True):.0f} s)', flush=True)
    tbode4.save_versions(outdir, 'v4 (exp 25 environment)')
    return data


def flows_over(sim, nsteps):
    """ Each transition by contact's total over the first nsteps steps, per person """
    r = sim.results.tb
    pop = float(sim.results.alive.values[0])
    return {f: float(np.sum(r[f].values[:nsteps]))/pop for f in contact_flows if f in r}


def q1_one(name, dkey, seed, n=1_000_000, nsteps=24):
    """ Agents over nsteps steps from the network ODE's equilibrium, and the network ODE from the same start """
    cal = sc.load(fn.q3)
    dt = dts[dkey]
    stop_ = start18 + nsteps*dt.years
    net = cal[f'{name}_net_{dkey}']
    sim = make_sim18(name, net.beta_edge, 'abm', dt, True, net.equilibrium, n=n, seed=seed, start=start18, stop=stop_)
    sim.run()
    out = sc.objdict(agents=flows_over(sim, nsteps), b=net.beta_edge)
    if seed == seed0:
        osim = make_sim18(name, net.beta_edge, 'ode', dt, True, net.equilibrium, start=start18, stop=stop_)
        osim.run()
        out.ode = flows_over(osim, nsteps)
    return out


def q1():
    data = sc.load(fn.q1) if fn.q1.exists() else sc.objdict()
    seeds = [seed0 + s for s in range(4)]
    for name in sets18:
        for dkey in dts:
            key = f'{name}_{dkey}'
            if key in data:
                continue
            T = sc.timer()
            res = sc.parallelize(q1_one, iterkwargs=[dict(name=name, dkey=dkey, seed=s) for s in seeds], ncpus=2) # 1M agents each: 2 at a time
            data[key] = sc.objdict(runs=res, ode=res[0].ode, b=res[0].b)
            sc.save(fn.q1, data)
            print(f'{key} ({T.toc(output=True):.0f} s)', flush=True)
    summarise_q1(data)
    return data


def summarise_q1(data=None):
    """ Agents ÷ network ODE for each transition by contact; predicted 1, with exp 18's prediction (the old ODE rule) alongside """
    data = sc.ifelse(data, sc.load(fn.q1))
    rows = []
    for key, d in data.items():
        name, dkey = key.split('_')
        D = dts[dkey].years
        for f, ode_val in d.ode.items():
            if ode_val <= 0:
                continue
            ag = np.array([r.agents[f] for r in d.runs])
            rr = rr_of(name, f)
            x = rr*d.b*D
            old = (-np.expm1(-x)/x)/(-np.expm1(-d.b*D)/(d.b*D)) # Exp 18: agents' factor over the susceptibles' factor, applied to every transition
            rows.append(dict(set=name, dt=dkey, flow=f, rr=rr, agents_over_ode=ag.mean()/ode_val, se=ag.std(ddof=1)/np.sqrt(len(ag))/ode_val,
                             predicted=1.0, exp18_rule=old, events_per_run=ag.mean()*1e6))
    df = pd.DataFrame(rows)
    df.to_csv(outdir / 'q1_summary.csv', index=False)
    with pd.option_context('display.width', 200, 'display.max_rows', 200):
        print(df.round(4).to_string(index=False))
    return df


def outputs18(sim):
    """ Yearly means: active prevalence, latent, cleared (all three no-longer-infected states), TB deaths """
    r = sim.results.tb
    alive = sim.results.alive.values
    get = lambda s: r[s].values if s in r else 0
    out = sc.objdict(year=np.array(sim.t.yearvec))
    out.prevalence = (get('asymptomatic') + get('symptomatic'))/alive
    out.latent = (get('infection') + get('infection_remote'))/alive
    out.cleared = (get('cleared') + get('recovered') + get('treated'))/alive
    out.tb_deaths = get('tb_death')/alive/sim.t.dt_year
    return req.yearly(out)


def q2_one(name, dkey, seed, solver='abm'):
    net = sc.load(fn.q3)[f'{name}_net_{dkey}']
    dt = dts[dkey]
    if solver == 'ode':
        sim = make_sim18(name, net.beta_edge, 'ode', dt, True, net.equilibrium)
    else:
        sim = make_sim18(name, net.beta_edge, 'abm', dt, True, net.equilibrium, n=n_agents, seed=seed)
    sim.run()
    return outputs18(sim)


def ratio(runs, ode, key, years=(1950, 2000)):
    sel = lambda o: (o.year >= years[0]) & (o.year < years[1])
    return np.array([o[key][sel(o)].mean() for o in runs])/ode[key][sel(ode)].mean()


def q2():
    data = sc.load(fn.q2) if fn.q2.exists() else sc.objdict(ode=sc.objdict(), abm=sc.objdict())
    for name in sets18:
        for dkey in dts:
            key = f'{name}_{dkey}'
            if key not in data.ode:
                data.ode[key] = q2_one(name, dkey, 0, solver='ode')
            if key in data.abm and len(data.abm[key]) >= reps18:
                continue
            T = sc.timer()
            data.abm[key] = sc.parallelize(q2_one, iterkwargs=[dict(name=name, dkey=dkey, seed=seed0 + 100 + s) for s in range(reps18)], ncpus=tbode4.ncpus)
            sc.save(fn.q2, data)
            r = ratio(data.abm[key], data.ode[key], 'prevalence')
            print(f'{key:10s} {reps18} reps in {T.toc(output=True):.0f} s; prevalence agents ÷ network ODE {r.mean():.3f} ± {r.std(ddof=1)/np.sqrt(reps18):.3f}', flush=True)
    summarise_q2(data)
    return data


def summarise_q2(data=None):
    data = sc.ifelse(data, sc.load(fn.q2))
    rows = []
    for key, runs in data.abm.items():
        name, dkey = key.split('_')
        for out in ['prevalence', 'latent', 'cleared', 'tb_deaths']:
            r = ratio(runs, data.ode[key], out)
            rows.append(dict(set=name, dt=dkey, output=out, reps=len(r), mean=r.mean(), se=r.std(ddof=1)/np.sqrt(len(r)),
                             p10=np.percentile(r, 10), p90=np.percentile(r, 90)))
    df = pd.DataFrame(rows)
    df.to_csv(outdir / 'q2_summary.csv', index=False)
    with pd.option_context('display.width', 200, 'display.max_rows', 200):
        print(df.round(4).to_string(index=False))
    return df


# %% Question 4: exp 22's conversion of exp 21's homogeneous sets

arms = sc.objdict(convk_month=('convk', 'month'), conv_month=('conv', 'month'), convk_week=('convk', 'week'), conv_week=('conv', 'week')) # Order fixes the seeds
reps22 = 20


def convert(pars, dt_years, which):
    """
    Per-contact network parameters for one homogeneous set (README): b and each rr' always; 'conv' also converts trans_asymp to
    kappa' = -ln(1 - kappa*q)/(b*D), 'convk' leaves it unchanged (exp 22's rule). Returns tbsim/fastode parameters, with beta per contact
    """
    D = dt_years
    b = tbode4.beta_edge_from_ode(pars['beta'], K, D)
    q = -np.expm1(-b*D)
    same = lambda v: -np.log1p(-v*q)/(b*D) # The relative rate whose per-step chance is v times a susceptible's
    out = dict(pars, beta=b)
    keys = ['rr_reinfection_cleared', 'rr_reinfection_rec', 'rr_reinfection_treat', 'rr_reinfection_inf'] + (['trans_asymp'] if which == 'conv' else [])
    for k in keys:
        out[k] = same(pars[k])
    return out


def full_state(x, p):
    """ Exp 21's saved equilibria predate fastode's treatment states (exp 23); with treatment off, those are 0 """
    return np.concatenate([x, np.zeros(fo.n_states(p) - len(x))])


def outputs_of(x, p):
    x = full_state(x, p)
    sm = fo.summary(x, p)
    return sc.objdict(prevalence=sm.prevalence, prop_asymptomatic=sm.prop_asymptomatic, inc_per_prev=sm.incidence/sm.prevalence, prev_to_inc=sm.prev_to_inc,
                      mortality=sm.mortality, mort_to_inc=sm.mort_to_inc, recent=x[fo.L], remote=x[fo.M], latent=x[fo.L] + x[fo.M], cleared=x[fo.C] + x[fo.R])


def make_sim22(s, netp, dkey, seed=0, solver='abm'):
    init, share = exp21.init_of(s.x, s.pars)
    rest = {k: v for k, v in netp.items() if k != 'beta'}
    tb = c.make_tb('A', netp['beta'], init=init, model='tb', **rest) # beta is per contact on the network
    kw = dict(start=exp21.start, stop=exp21.stop, dt=dts[dkey], rand_seed=seed, verbose=0, networks=ss.RandomNet(n_contacts=ss.constant(K), dur=0))
    if solver == 'abm':
        return ss.Sim(tb, solver='abm', n_agents=n_agents, total_pop=n_agents, analyzers=exp19.Incidence(share_flagged=share), **kw)
    return ss.Sim(tb, solver='ode', total_pop=1e5, **kw)


def conversions():
    """ Each set's converted parameters, the standalone network ODE's prediction (kappa inside the exponent), and Starsim's network ODE """
    sets = sc.load(exp21.fn.sets)
    conv = sc.objdict()
    rows = []
    for i, s in enumerate(sets):
        ref = outputs_of(s.x, fo.make_pars(**s.pars))
        for arm, (which, dkey) in arms.items():
            D = dts[dkey].years
            netp = convert(s.pars, D, which)
            p = fo.make_pars(**netp, net_contacts=K, net_step=D, net_source_exp=1)
            pred = outputs_of(fo.equilibrium(p, x0=full_state(s.x, p)), p)
            osim = make_sim22(s, netp, dkey, solver='ode')
            osim.run()
            ss_ode = exp21.yearly_outputs(osim, False)
            ss_prev = exp21.derived(ss_ode, 'prevalence')
            conv[f's{i}_{arm}'] = sc.objdict(pars=netp, predicted=pred, predicted_ratio={k: pred[k]/ref[k] for k in ref if ref[k]}, starsim_ode=ss_ode)
            rows.append(dict(set=i, arm=arm, beta=s.pars['beta'], b=netp['beta'], pi=s.pars['rr_reinfection_cleared'], pi_net=netp['rr_reinfection_cleared'],
                             kappa=s.pars['trans_asymp'], kappa_net=netp['trans_asymp'], pred_prevalence_ratio=pred.prevalence/ref.prevalence,
                             starsim_ode_over_pred=ss_prev/pred.prevalence))
            print(f's{i} {arm}: predicted prevalence ÷ homogeneous {pred.prevalence/ref.prevalence:.4f}; Starsim network ODE ÷ predicted {ss_prev/pred.prevalence:.4f}', flush=True)
    sc.save(fn.conv, conv)
    df = pd.DataFrame(rows)
    df.to_csv(outdir / 'conversions.csv', index=False)
    with pd.option_context('display.width', 200, 'display.max_rows', 100):
        print(df.round(4).to_string(index=False))
    return conv


def convert_one(i, arm, seed):
    sets = sc.load(exp21.fn.sets)
    conv = sc.load(fn.conv)
    sim = make_sim22(sets[i], conv[f's{i}_{arm}'].pars, arms[arm][1], seed)
    sim.run()
    return exp21.yearly_outputs(sim, True)


def run_convert(dkeys=None):
    """ Question 4's agents, for the arms at the given steps (e.g. ['month']; default all) """
    if not fn.conv.exists():
        conversions()
    sets = sc.load(exp21.fn.sets)
    data = sc.load(fn.runs) if fn.runs.exists() else sc.objdict(abm=sc.objdict())
    for i in range(len(sets)):
        for a, arm in enumerate(arms):
            if dkeys and arms[arm][1] not in dkeys:
                continue
            key = f's{i}_{arm}'
            if key in data.abm and len(data.abm[key]) >= reps22:
                continue
            T = sc.timer()
            data.abm[key] = sc.parallelize(convert_one, iterkwargs=[dict(i=i, arm=arm, seed=seed0 + 1000 + 1000*i + 100*a + s) for s in range(reps22)], ncpus=tbode4.ncpus)
            sc.save(fn.runs, data)
            print(f'{key}: {reps22} reps in {T.toc(output=True):.0f} s', flush=True)
    summarise_runs(fn.runs, 'conv')
    return data


# %% Question 5: exp 21's homogeneous monthly arm

def homog_one(i, seed):
    sets = sc.load(exp21.fn.sets)
    sim = exp21.make_sim(sets[i], 'abm', 'month', seed)
    sim.run()
    return exp21.yearly_outputs(sim, True)


def run_homog():
    sets = sc.load(exp21.fn.sets)
    data = sc.load(fn.homog) if fn.homog.exists() else sc.objdict(abm=sc.objdict())
    for i in range(len(sets)):
        key = f's{i}_homog_month'
        if key in data.abm and len(data.abm[key]) >= reps22:
            continue
        T = sc.timer()
        data.abm[key] = sc.parallelize(homog_one, iterkwargs=[dict(i=i, seed=seed0 + 20000 + 1000*i + s) for s in range(reps22)], ncpus=tbode4.ncpus)
        sc.save(fn.homog, data)
        print(f'{key}: {reps22} reps in {T.toc(output=True):.0f} s', flush=True)
    summarise_runs(fn.homog, 'homog')
    return data


def summarise_runs(path, label):
    """ Agents ÷ the calibrated homogeneous ODE (exp 21's v4 ODE runs, and the standalone ODE for first-entry incidence), by set and pooled """
    sets = sc.load(exp21.fn.sets)
    ode21 = sc.load(exp21.fn.runs).ode
    conv = sc.load(fn.conv) if fn.conv.exists() else {}
    data = sc.load(path)
    rows = []
    for key, runs in data.abm.items():
        i, arm = int(key.split('_')[0][1:]), key.split('_', 1)[1]
        s, dkey = sets[i], arm.split('_')[-1]
        for out in exp21.outputs_compared:
            ref = exp21.ode_value(s, ode21[f's{i}_{dkey}'], out)
            if not ref:
                continue
            r = np.array([exp21.derived(o, out) for o in runs])/ref
            pred = conv[key].predicted_ratio.get(out, np.nan) if key in conv else 1.0
            rows.append(dict(set=i, arm=arm, output=out, reps=len(r), mean=r.mean(), se=r.std(ddof=1)/np.sqrt(len(r)), p10=np.percentile(r, 10),
                             p90=np.percentile(r, 90), predicted=pred))
    df = pd.DataFrame(rows)
    df.to_csv(outdir / f'{label}_by_set.csv', index=False)
    g = df.groupby(['arm', 'output'])
    pooled = pd.DataFrame(dict(sets=g['mean'].size(), mean=g['mean'].mean(), se_between=g['mean'].std(ddof=1)/np.sqrt(g['mean'].size()),
                               min_set=g['mean'].min(), max_set=g['mean'].max(), predicted=g['predicted'].mean(),
                               max_abs_z=g.apply(lambda x: ((x['mean'] - 1)/x['se']).abs().max()),
                               max_abs_z_vs_pred=g.apply(lambda x: ((x['mean'] - x['predicted'])/x['se']).abs().max()))).reset_index()
    pooled.to_csv(outdir / f'{label}_pooled.csv', index=False)
    with pd.option_context('display.width', 220, 'display.max_rows', 300):
        print(pooled.round(4).to_string(index=False))
    return df, pooled


# %% Question 3: exp 22's per-step diagnostic

targets = sc.objdict(transmission=('susceptible', None), reinfect_cle=('cleared', 'rr_reinfection_cleared'), reinfect_rec=('recovered', 'rr_reinfection_rec'),
                     reinfect_trd=('treated', 'rr_reinfection_treat'), reinfect_remote=('infection_remote', 'rr_reinfection_inf'))


class Snapshot(ss.Intervention):
    """ The edges and states transmission is about to use, with the expected infections by class under each rule """
    def __init__(self, netp, homp, D, **kwargs):
        super().__init__(**kwargs)
        self.netp, self.homp, self.D = netp, homp, D
        self.rows = []
        self.pending = None
        return

    def step(self):
        sim = self.sim
        tb = sim.diseases.tb
        n = int(np.asarray(sim.people.uid).max()) + 1
        alive = np.asarray(sim.people.alive.uids)
        A, Y = np.zeros(n, bool), np.zeros(n, bool)
        A[np.asarray(tb.asymptomatic.uids)] = True
        Y[np.asarray(tb.symptomatic.uids)] = True
        e = sim.networks[0].edges
        p1, p2 = np.asarray(e.p1), np.asarray(e.p2)
        nA = np.bincount(p1[A[p2]], minlength=n) + np.bincount(p2[A[p1]], minlength=n) # Asymptomatic contacts of each agent
        nY = np.bincount(p1[Y[p2]], minlength=n) + np.bincount(p2[Y[p1]], minlength=n)
        deg = np.bincount(p1, minlength=n) + np.bincount(p2, minlength=n)
        b, D, N = self.netp['beta'], self.D, len(alive)
        kap = self.netp['trans_asymp']
        row = dict(ti=sim.ti, N=N, A=int(A.sum()), Y=int(Y.sum()), mean_deg=deg[alive].mean(), var_deg=deg[alive].var(), self_loops=int((p1 == p2).sum()))
        for tname, (state, rrkey) in targets.items():
            if not hasattr(tb, state):
                continue
            u = np.asarray(getattr(tb, state).uids)
            rr_n = self.netp[rrkey] if rrkey else 1.0
            row[f'{tname}_exp_sat'] = -np.expm1(-rr_n*b*D*(kap*nA[u] + nY[u])).sum() # kappa inside the exponent
            qY = -np.expm1(-rr_n*b*D)
            row[f'{tname}_exp_lin'] = (1 - (1 - qY)**nY[u]*(1 - kap*qY)**nA[u]).sum() # kappa scaling the chance (the old pins)
            row[f'{tname}_uids'] = u
        self.pending = row
        return


class Count(ss.Analyzer):
    """ The infections by class in the step, matched to its snapshot """
    def step(self):
        snap = self.sim.interventions[0]
        row = snap.pending
        if row is None:
            return
        tb = self.sim.diseases.tb
        for tname in targets:
            if f'{tname}_uids' not in row:
                continue
            u = np.asarray(getattr(tb, tname).uids)
            row[f'{tname}_obs'] = int(np.isin(u, row.pop(f'{tname}_uids')).sum()) # Those infected by this transition who were in the class at the snapshot
        snap.rows.append(row)
        snap.pending = None
        return


def diag_one(i, seed, which='conveq'):
    """ One run of the diagnostic, on exp 22's kappa' equilibrium-matched parameters for set i (monthly) """
    sets = sc.load(exp21.fn.sets)
    conv22 = sc.load(root / 'experiments' / '22_network_conversion' / 'conversions.obj')
    s, netp = sets[i], conv22[f's{i}_{which}_month'].pars
    sim = make_sim22(s, netp, 'month', seed)
    sim.pars.interventions = Snapshot(netp, s.pars, dts.month.years)
    sim.pars.analyzers = [sim.pars.analyzers, Count()]
    sim.run()
    df = pd.DataFrame(sim.interventions[0].rows)
    df['set'], df['seed'] = i, seed
    return df


def diagnose():
    kw = [dict(i=i, seed=seed0 + 30000 + s) for i in [3, 7] for s in range(2)]
    dfs = sc.parallelize(diag_one, iterkwargs=kw, ncpus=tbode4.ncpus)
    df = pd.concat(dfs)
    sc.save(fn.diag, df.drop(columns=[col for col in df.columns if col.endswith('_uids')], errors='ignore'))
    return report_diag(df)


def report_diag(df=None):
    df = sc.ifelse(df, sc.load(fn.diag) if fn.diag.exists() else None)
    d = df[df.ti > 0]
    print(f"Mean degree {d.mean_deg.mean():.4f} (variance {d.var_deg.mean():.3f}); self-loops per step {d.self_loops.mean():.1f}; N {d.N.mean():.0f}")
    rows = []
    for i, g in d.groupby('set'):
        for tname in targets:
            if f'{tname}_obs' not in g:
                continue
            obs, es, el = g[f'{tname}_obs'].sum(), g[f'{tname}_exp_sat'].sum(), g[f'{tname}_exp_lin'].sum()
            if not obs:
                continue
            rows.append(dict(set=i, transition=tname, observed=obs, obs_per_sat=obs/es, obs_per_lin=obs/el, se=np.sqrt(obs)/es))
    out = pd.DataFrame(rows)
    out.to_csv(outdir / 'diagnose_summary.csv', index=False)
    with pd.option_context('display.width', 200):
        print(out.round(4).to_string(index=False))
    return out


def summary():
    for name, path, func in [('q1', fn.q1, summarise_q1), ('q2', fn.q2, summarise_q2), ('diagnose', fn.diag, report_diag)]:
        if path.exists():
            print(f'\n## {name}')
            func()
    if fn.runs.exists():
        print('\n## conversion')
        summarise_runs(fn.runs, 'conv')
    if fn.homog.exists():
        print('\n## homogeneous')
        summarise_runs(fn.homog, 'homog')
    return


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'q3'
    if which == 'convert':
        run_convert(sys.argv[2:] or None) # e.g. `convert month`: only the monthly arms
    else:
        dict(q1=q1, q2=q2, q3=q3, diagnose=diagnose, conversions=conversions, homog=run_homog, summary=summary)[which]()
