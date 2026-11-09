"""
Experiment 22, v4 environment: exp 21's 12 calibrated parameter sets, converted for ss.RandomNet(n_contacts=10, dur=0), as agents on the
network against the calibrated homogeneous ODE. See README.md. Run from the repo root:

    PYTHONPATH=v4 v4/.venv/bin/python experiments/22_network_conversion/run.py          # conversions, predictions, agents (resumes)
    PYTHONPATH=v4 v4/.venv/bin/python experiments/22_network_conversion/run.py summary  # summaries and figures from the saved runs
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
import tbode4 # noqa: E402
from tbode4 import fastode as fo # noqa: E402


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


exp21 = load_module('exp21_run', root / 'experiments' / '21_calibrated_comparison' / 'run.py') # Sets, ODE runs, outputs, analyzer
exp19 = exp21.exp19
calib = exp21.calib

K = 10
dts = sc.objdict(month=ss.months(1), week=ss.weeks(1))
arms = sc.objdict(conv_month=('conv', 'month'), conveq_month=('conveq', 'month'), beta_month=('beta', 'month'), conveq_week=('conveq', 'week'),
                  convk_month=('convk', 'month'), convk_week=('convk', 'week')) # Order fixes the seeds: new arms go at the end
reps, seed0, n_agents = 20, 22000, 50_000
fn = sc.objdict(conv=outdir / 'conversions.obj', runs=outdir / 'runs.obj')
outputs_compared = exp21.outputs_compared


def convert(pars, dt_years, which='conv', x=None):
    """
    Per-contact network parameters for one homogeneous set (README). Returns tbsim/fastode parameters, with beta per contact.

    - 'convk': b and pi', with kappa unchanged. v4's agents scale a contact's chance linearly by the source's infectiousness, so this
      matches every combination of target and source per contact and per step: the conversion to use;
    - 'conv': b, pi' and kappa' = -ln(1 - kappa*q)/(b*D), on the assumption that the source's infectiousness saturates too (it doesn't, so
      transmission from asymptomatic sources is kappa'/kappa too low);
    - 'conveq': b and kappa' as 'conv', and pi' chosen so the total force of reinfection at the homogeneous equilibrium x matches under that
      same assumption;
    - 'beta': b only (CONVENTIONS.md 1).
    """
    from scipy.optimize import brentq
    D = dt_years
    b = tbode4.beta_edge_from_ode(pars['beta'], K, D)
    q = -np.expm1(-b*D) # A susceptible's chance per symptomatic contact per step
    out = dict(pars, beta=b)
    rr = ['rr_reinfection_cleared', 'rr_reinfection_rec', 'rr_reinfection_treat', 'rr_reinfection_inf']
    same = lambda v: -np.log1p(-v*q)/(b*D) # The relative rate whose per-step chance is v times a susceptible's
    if which == 'convk':
        for k in rr:
            out[k] = same(pars[k])
    if which in ('conv', 'conveq'):
        for k in rr + ['trans_asymp']:
            out[k] = same(pars[k])
    if which == 'conveq':
        pi, kappa, kap_n = pars['rr_reinfection_cleared'], pars['trans_asymp'], out['trans_asymp']
        A, Y = x[fo.A], x[fo.Y]
        target = pi*pars['beta']*(kappa*A + Y) # The homogeneous force of reinfection at the equilibrium
        net = lambda p_: (K/D)*(-np.expm1(-p_*b*D)*Y - np.expm1(-p_*kap_n*b*D)*A)
        pi_n = brentq(lambda p_: net(p_) - target, 0.5*pi, 2*pi, xtol=1e-14)
        for k in rr:
            out[k] = pi_n
    return out


def predict(s, net_pars, D):
    """ The network standalone ODE's equilibrium outputs at net_pars (a reference version of proposal 03's fix) """
    p = fo.make_pars(**net_pars, net_contacts=K, net_step=D)
    x = fo.equilibrium(p, x0=s.x)
    sm = fo.summary(x, p)
    return sc.objdict(prevalence=sm.prevalence, prop_asymptomatic=sm.prop_asymptomatic, inc_per_prev=sm.incidence/sm.prevalence, prev_to_inc=sm.prev_to_inc,
                      mortality=sm.mortality, mort_to_inc=sm.mort_to_inc, recent=x[fo.L], remote=x[fo.M], latent=x[fo.L] + x[fo.M],
                      cleared=x[fo.C] + x[fo.R])


def homog(s):
    """ The homogeneous ODE's equilibrium outputs (standalone ODE), the reference """
    sm, x = s.ode, s.x
    return sc.objdict(prevalence=sm.prevalence, prop_asymptomatic=sm.prop_asymptomatic, inc_per_prev=sm.incidence/sm.prevalence, prev_to_inc=sm.prev_to_inc,
                      mortality=sm.mortality, mort_to_inc=sm.mort_to_inc, recent=x[fo.L], remote=x[fo.M], latent=x[fo.L] + x[fo.M],
                      cleared=x[fo.C] + x[fo.R])


def conversions():
    sets = sc.load(exp21.fn.sets)
    conv = sc.objdict()
    rows = []
    for i, s in enumerate(sets):
        ref = homog(s)
        for arm, (which, dkey) in arms.items():
            D = dts[dkey].years
            netp = convert(s.pars, D, which, x=s.x)
            pred = predict(s, netp, D)
            conv[f's{i}_{arm}'] = sc.objdict(pars=netp, predicted=pred, predicted_ratio={k: pred[k]/ref[k] for k in ref if ref[k]})
            rows.append(dict(set=i, arm=arm, beta=s.pars['beta'], b=netp['beta'], pi=s.pars['rr_reinfection_cleared'], pi_net=netp['rr_reinfection_cleared'],
                             kappa=s.pars['trans_asymp'], kappa_net=netp['trans_asymp'], pred_prevalence_ratio=pred.prevalence/ref.prevalence))
    sc.save(fn.conv, conv)
    df = pd.DataFrame(rows)
    df.to_csv(outdir / 'conversions.csv', index=False)
    with pd.option_context('display.width', 200, 'display.max_rows', 100):
        print(df.round(4).to_string(index=False))
    return conv


def make_sim(s, netp, dkey, seed):
    init, share = exp21.init_of(s.x, s.pars)
    rest = {k: v for k, v in netp.items() if k != 'beta'}
    tb = c.make_tb('A', netp['beta'], init=init, model='tb', **rest) # beta is per contact on the network
    return ss.Sim(tb, solver='abm', n_agents=n_agents, total_pop=n_agents, start=exp21.start, stop=exp21.stop, dt=dts[dkey], rand_seed=seed, verbose=0,
                  networks=ss.RandomNet(n_contacts=ss.constant(K), dur=0), analyzers=exp19.Incidence(share_flagged=share))


def run_one(i, arm, seed):
    sets = sc.load(exp21.fn.sets)
    conv = sc.load(fn.conv)
    sim = make_sim(sets[i], conv[f's{i}_{arm}'].pars, arms[arm][1], seed)
    sim.run()
    return exp21.yearly_outputs(sim, True)


def main():
    if not fn.conv.exists():
        conversions()
    sets = sc.load(exp21.fn.sets)
    data = sc.load(fn.runs) if fn.runs.exists() else sc.objdict(abm=sc.objdict())
    for i in range(len(sets)):
        for a, arm in enumerate(arms):
            key = f's{i}_{arm}'
            if key in data.abm and len(data.abm[key]) >= reps:
                continue
            T = sc.timer()
            data.abm[key] = sc.parallelize(run_one, iterkwargs=[dict(i=i, arm=arm, seed=seed0 + 1000*i + 100*a + s) for s in range(reps)], ncpus=tbode4.ncpus)
            sc.save(fn.runs, data)
            print(f'{key}: {reps} reps in {T.toc(output=True):.0f} s', flush=True)
    tbode4.save_versions(outdir, 'v4')
    summarise()
    return data


def summarise():
    sets = sc.load(exp21.fn.sets)
    ode21 = sc.load(exp21.fn.runs).ode # v4's homogeneous ODE runs (exp 21), at each step
    conv = sc.load(fn.conv)
    data = sc.load(fn.runs)
    rows = []
    for key, runs in data.abm.items():
        i, arm = int(key.split('_')[0][1:]), key.split('_', 1)[1]
        s, dkey = sets[i], arms[arm][1]
        o_ode = ode21[f's{i}_{dkey}']
        for out in outputs_compared:
            ref = exp21.ode_value(s, o_ode, out)
            if not ref:
                continue
            r = np.array([exp21.derived(o, out) for o in runs])/ref
            rows.append(dict(set=i, arm=arm, output=out, reps=len(r), mean=r.mean(), se=r.std(ddof=1)/np.sqrt(len(r)), p10=np.percentile(r, 10),
                             p90=np.percentile(r, 90), predicted=conv[key].predicted_ratio.get(out, np.nan)))
    df = pd.DataFrame(rows)
    df.to_csv(outdir / 'by_set.csv', index=False)
    g = df.groupby(['arm', 'output'])
    pooled = pd.DataFrame(dict(sets=g['mean'].size(), mean=g['mean'].mean(), se_between=g['mean'].std(ddof=1)/np.sqrt(g['mean'].size()),
                               min_set=g['mean'].min(), max_set=g['mean'].max(), predicted=g['predicted'].mean(),
                               max_abs_z=g.apply(lambda x: ((x['mean'] - 1)/x['se']).abs().max()),
                               max_abs_z_vs_pred=g.apply(lambda x: ((x['mean'] - x['predicted'])/x['se']).abs().max()))).reset_index()
    pooled.to_csv(outdir / 'pooled.csv', index=False)
    with pd.option_context('display.width', 220, 'display.max_rows', 300):
        print(pooled.round(4).to_string(index=False))
    plot(df)
    return df, pooled


def plot(df):
    shown = ['prevalence', 'prop_asymptomatic', 'inc_per_prev', 'mortality', 'remote', 'cleared']
    colors = dict(convk_month='#1f4e8c', convk_week='#2e8b57', beta_month='#c4620b', conv_month='#b9b9b9', conveq_month='#8c8c8c', conveq_week='#5e5e5e')
    labels = dict(convk_month='Converted (b, π′; κ unchanged), monthly', convk_week='Converted (b, π′; κ unchanged), weekly', beta_month='beta only (our convention), monthly',
                  conv_month='b, π′ and κ′, monthly', conveq_month='b, κ′, π′ equilibrium-matched, monthly', conveq_week='b, κ′, π′ equilibrium-matched, weekly')
    fig, axes = plt.subplots(2, 3, figsize=(15, 7.5), sharex=True)
    for ax, out in zip(axes.flat, shown):
        for j, arm in enumerate(labels):
            d = df[(df.output == out) & (df.arm == arm)].sort_values('set')
            x = d.set + (j - (len(arms) - 1)/2)*0.13
            ax.errorbar(x, d['mean'], yerr=2*d['se'], fmt='o', ms=4, capsize=2, color=colors[arm], label=labels[arm])
            ax.plot(x, d['predicted'], '_', ms=10, mew=2, color=colors[arm])
        ax.axhline(1, color='0.4', lw=0.8)
        ax.set_title(out, fontsize=9)
        ax.set_xticks(range(12))
    for ax in axes[1]:
        ax.set_xlabel('Parameter set (exp 21)')
    for ax in axes[:, 0]:
        ax.set_ylabel('Agents on the network ÷ homogeneous ODE')
    axes[0][0].legend(fontsize=7, frameon=False)
    fig.suptitle('Exp 22: agents on RandomNet(10, dur=0) ÷ the calibrated homogeneous ODE (mean ± 2 SE, 20 reps; dashes: the network standalone ODE\'s prediction)', fontsize=10)
    fig.tight_layout()
    fig.savefig(outdir / 'agents_vs_ode.png', dpi=150)
    print('Saved agents_vs_ode.png')
    return


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    if which == 'conversions':
        conversions()
    elif which == 'summary':
        summarise()
    else:
        main()
