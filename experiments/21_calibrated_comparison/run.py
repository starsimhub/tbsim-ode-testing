"""
Experiment 21, v4 environment: v4 agents against v4's ODE at 12 parameter sets drawn from exp 20's calibrated posterior. See README.md.
Run from the repo root:

    PYTHONPATH=v4 v4/.venv/bin/python experiments/21_calibrated_comparison/run.py          # draw the sets (once), run everything (resumes)
    PYTHONPATH=v4 v4/.venv/bin/python experiments/21_calibrated_comparison/run.py summary  # summaries and figures from the saved runs
"""
import sys
import importlib.util
import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss
import matplotlib.pyplot as plt

outdir = sc.thispath(__file__)
root = outdir.parent.parent
sys.path.insert(0, str(root / 'v4'))
import common as c # noqa: E402
import run_equilibrium as req # noqa: E402
import tbode4 # noqa: E402
from tbode4 import fastode as fo # noqa: E402


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


calib = load_module('exp20_calib', root / 'experiments' / '20_ode_calibration' / 'calib.py') # Priors, parameter mapping, targets
exp19 = load_module('exp19_run', root / 'experiments' / '19_v4_incidence' / 'run.py')       # The first-entry incidence analyzer

n_sets, draw_seed = 12, 21
dts = sc.objdict(month=ss.months(1), week=ss.weeks(1))
reps, seed0, n_agents = 20, 21000, 50_000
start, stop = 1900, 2000
fn = sc.objdict(sets=outdir / 'sets.obj', runs=outdir / 'runs.obj')
outputs_compared = ['prevalence', 'prop_asymptomatic', 'inc_per_prev', 'prev_to_inc', 'mortality', 'mort_to_inc', 'latent', 'recent', 'remote', 'cleared']


def draw_sets():
    """ 12 distinct draws from exp 20's posterior (weighted resampling without replacement), with their parameters and ODE equilibria """
    res = sc.load(root / 'experiments' / '20_ode_calibration' / 'real.obj')
    rng = np.random.default_rng(draw_seed)
    idx = rng.choice(len(res.weights), size=n_sets, replace=False, p=res.weights)
    sets = []
    for i in idx:
        theta = res.theta[i]
        pars = calib.to_pars(theta)
        p = fo.make_pars(**pars)
        x = fo.equilibrium(p)
        sets.append(sc.objdict(index=int(i), weight=float(res.weights[i]), theta=dict(zip(calib.names, map(float, theta))), pars=pars, x=x,
                               ode=fo.summary(x, p), flows=fo.flows(x, p), loglik=float(res.ll[i])))
    sc.save(fn.sets, sets)
    df = pd.DataFrame([s.theta for s in sets])
    post = pd.read_csv(root / 'experiments' / '20_ode_calibration' / 'real_posterior.csv').set_index('parameter')
    print('Drawn sets (rows) against the posterior 95% interval:')
    with pd.option_context('display.width', 250):
        print(df.round(3).to_string())
        print(pd.DataFrame({'min drawn': df.min(), 'max drawn': df.max(), 'post 2.5%': post.q025, 'post 97.5%': post.q975}).round(3).to_string())
    return sets


def init_of(x, pars):
    init = dict(zip(fo.state_names, x))
    non0, non1 = init.pop('non0'), init.pop('non1')
    init['non_infectious'] = non0 + non1
    if not pars.get('k_latent'):
        init.pop('infection_remote')
    return init, non1/(non0 + non1)


def make_sim(s, solver, dkey, seed=0):
    init, share = init_of(s.x, s.pars)
    rest = {k: v for k, v in s.pars.items() if k != 'beta'}
    tb = c.make_tb('A', s.pars['beta'], init=init, model='tb', **rest) # Set 'A' is empty: every parameter comes from the drawn set
    kw = dict(start=start, stop=stop, dt=dts[dkey], rand_seed=seed, verbose=0)
    if solver == 'abm':
        return ss.Sim(tb, solver='abm', n_agents=n_agents, total_pop=n_agents, analyzers=exp19.Incidence(share_flagged=share), **kw)
    return ss.Sim(tb, solver='ode', total_pop=1e5, **kw)


def yearly_outputs(sim, agents):
    r = sim.results.tb
    alive = sim.results.alive.values
    dt = sim.t.dt_year
    get = lambda k: r[k].values if k in r else 0
    asy, sym = get('asymptomatic'), get('symptomatic')
    out = sc.objdict(year=np.array(sim.t.yearvec), prevalence=(asy + sym)/alive, asymptomatic=asy/alive, mortality=get('tb_death')/alive/dt,
                     recent=get('infection')/alive, remote=get('infection_remote')/alive,
                     cleared=(get('cleared') + get('recovered') + get('treated'))/alive)
    if agents:
        out.incidence = sim.analyzers[0].first/alive/dt
    return req.yearly(out)


def run_one(i, dkey, seed):
    sets = sc.load(fn.sets)
    sim = make_sim(sets[i], 'abm', dkey, seed)
    sim.run()
    return yearly_outputs(sim, True)


def run_ode(i, dkey):
    sets = sc.load(fn.sets)
    sim = make_sim(sets[i], 'ode', dkey)
    sim.run()
    return yearly_outputs(sim, False)


def main():
    sets = sc.load(fn.sets) if fn.sets.exists() else draw_sets()
    data = sc.load(fn.runs) if fn.runs.exists() else sc.objdict(ode=sc.objdict(), abm=sc.objdict())
    for i in range(n_sets):
        for dkey in dts:
            key = f's{i}_{dkey}'
            if key not in data.ode:
                data.ode[key] = run_ode(i, dkey)
            if key in data.abm and len(data.abm[key]) >= reps:
                continue
            T = sc.timer()
            data.abm[key] = sc.parallelize(run_one, iterkwargs=[dict(i=i, dkey=dkey, seed=seed0 + 1000*i + s) for s in range(reps)], ncpus=tbode4.ncpus)
            sc.save(fn.runs, data)
            print(f'{key}: {reps} reps in {T.toc(output=True):.0f} s', flush=True)
    tbode4.save_versions(outdir, 'v4')
    summarise()
    return data


def window_mean(o, key):
    sel = (o.year >= 1950) & (o.year < 2000)
    return o[key][sel].mean()


def derived(o, key):
    """ An output's 1950–2000 mean; ratios are ratios of the means """
    m = lambda k: window_mean(o, k)
    if key == 'prop_asymptomatic':
        return m('asymptomatic')/m('prevalence')
    if key == 'inc_per_prev':
        return m('incidence')/m('prevalence')
    if key == 'prev_to_inc':
        return m('prevalence')/m('incidence')
    if key == 'mort_to_inc':
        return m('mortality')/m('incidence')
    if key == 'latent':
        return m('recent') + m('remote')
    return m(key)


def ode_value(s, o_ode, key):
    """ The ODE's value: v4's ODE run for state shares and mortality; the standalone ODE for anything that needs first-entry incidence """
    e = s.ode
    first = dict(inc_per_prev=e.incidence/e.prevalence, prev_to_inc=e.prev_to_inc, mort_to_inc=e.mort_to_inc)
    return first[key] if key in first else derived(o_ode, key)


def summarise():
    sets = sc.load(fn.sets)
    data = sc.load(fn.runs)
    rows = []
    for key, runs in data.abm.items():
        i, dkey = int(key.split('_')[0][1:]), key.split('_')[1]
        s, o_ode = sets[i], data.ode[key]
        check = abs(derived(o_ode, 'prevalence')/s.ode.prevalence - 1) # v4's ODE stays at the standalone equilibrium
        for out in outputs_compared:
            ref = ode_value(s, o_ode, out)
            if not ref:
                continue
            r = np.array([derived(o, out) for o in runs])/ref
            rows.append(dict(set=i, dt=dkey, output=out, ode=ref, reps=len(r), mean=r.mean(), se=r.std(ddof=1)/np.sqrt(len(r)), p10=np.percentile(r, 10),
                             p90=np.percentile(r, 90), ode_drift=check))
    df = pd.DataFrame(rows)
    df.to_csv(outdir / 'by_set.csv', index=False)
    pooled = (df.groupby(['dt', 'output']).agg(sets=('mean', 'size'), mean=('mean', 'mean'), se_between=('mean', lambda x: x.std(ddof=1)/np.sqrt(len(x))),
                                                 min_set=('mean', 'min'), max_set=('mean', 'max'), max_abs_z=('mean', lambda x: 0.0)).reset_index())
    for k, g in df.groupby(['dt', 'output']):
        z = ((g['mean'] - 1)/g['se']).abs().max()
        pooled.loc[(pooled.dt == k[0]) & (pooled.output == k[1]), 'max_abs_z'] = z
    pooled.to_csv(outdir / 'pooled.csv', index=False)
    print(f'v4 ODE vs standalone equilibrium, max drift in prevalence: {df.ode_drift.max():.1e}')
    with pd.option_context('display.width', 200, 'display.max_rows', 300):
        print(pooled.round(4).to_string(index=False))
    plot(df, sets, data)
    return df, pooled


def plot(df, sets, data):
    """ Agents ÷ ODE by set for each output, and the agents' target outputs against the target intervals """
    shown = ['prevalence', 'prop_asymptomatic', 'inc_per_prev', 'mortality', 'latent', 'cleared']
    colors = dict(month='#3b6ea8', week='#c4620b')
    fig, axes = plt.subplots(2, 3, figsize=(14, 7), sharex=True)
    for ax, out in zip(axes.flat, shown):
        for j, dkey in enumerate(dts):
            d = df[(df.output == out) & (df.dt == dkey)].sort_values('set')
            ax.errorbar(d.set + (j - 0.5)*0.25, d['mean'], yerr=2*d['se'], fmt='o', ms=4, capsize=2, color=colors[dkey], label=f'{dkey}ly dt')
        ax.axhline(1, color='0.4', lw=0.8)
        ax.set_title(out, fontsize=9)
        ax.set_xticks(range(n_sets))
    axes[0][0].set_ylabel('Agents ÷ ODE, 1950–2000')
    axes[1][0].set_ylabel('Agents ÷ ODE, 1950–2000')
    for ax in axes[1]:
        ax.set_xlabel('Parameter set (drawn from exp 20\'s posterior)')
    axes[0][0].legend(fontsize=8, frameon=False)
    fig.suptitle('Exp 21: v4 agents ÷ v4 ODE at 12 calibrated parameter sets (mean ± 2 SE, 20 reps; incidence = first entry per episode, per prevalent case)', fontsize=10)
    fig.tight_layout()
    fig.savefig(outdir / 'agents_vs_ode.png', dpi=150)

    # The targets: agents' and ODE's values per set, against the target intervals
    tg = dict(prevalence=('prevalence', 1e5), prop_asymptomatic=('prop_asymptomatic', 1), prev_to_inc=('prev_to_inc', 1))
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.6))
    for ax, (tkey, (out, sc_)) in zip(axes, tg.items()):
        t = calib.targets[tkey]
        if t['kind'] == 'binomial':
            from scipy.stats import beta as betad
            m, lo, hi = t['x']/t['n'], *betad.ppf([0.025, 0.975], t['x'] + 0.5, t['n'] - t['x'] + 0.5)
        else:
            m, lo, hi = t['mean'], t['mean']*np.exp(-1.96*t['sigma']), t['mean']*np.exp(1.96*t['sigma'])
        ax.axhspan(lo*sc_, hi*sc_, color='0.88', zorder=0)
        ax.axhline(m*sc_, color='0.3', lw=1)
        for i, s in enumerate(sets):
            ode = ode_value(s, data.ode[f's{i}_month'], out)*sc_
            ag = np.array([derived(o, out) for o in data.abm[f's{i}_month']])*sc_
            ax.plot(i, ode, 'x', color='0.1', ms=6, label='ODE' if i == 0 else None)
            ax.errorbar(i + 0.2, ag.mean(), yerr=2*ag.std(ddof=1)/np.sqrt(len(ag)), fmt='o', ms=4, color='#3b6ea8', capsize=2, label='Agents, monthly (mean ± 2 SE)' if i == 0 else None)
        ax.set_title(tkey, fontsize=9)
        ax.set_xticks(range(n_sets))
        ax.set_xlabel('Parameter set')
    axes[0].legend(fontsize=7, frameon=False)
    fig.suptitle('Exp 21: the targets at each calibrated set, ODE and agents, against the target mean and 95% interval (grey)', fontsize=10)
    fig.tight_layout()
    fig.savefig(outdir / 'targets.png', dpi=150)
    print('Saved agents_vs_ode.png, targets.png')
    return


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'summary':
        summarise()
    elif len(sys.argv) > 1 and sys.argv[1] == 'sets':
        draw_sets()
    else:
        main()
