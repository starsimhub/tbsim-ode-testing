"""
Experiment 11: TBExact incidence counted from transitions (event-based) vs end-of-step snapshots, against the ODE. See README.md.
"""

import os
import sys
import numpy as np
import pandas as pd
import sciris as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import tbode

outdir = os.path.dirname(os.path.abspath(__file__))
dir10 = os.path.join(outdir, '..', '10_transmission_step')
sys.path.insert(0, dir10)
from step_map import settings  # noqa: E402 (exp 10's parameter settings: k0, A, B)

start, stop = 1900, 2000
dts = {'1 week': 7/365, '1 month': 1/12}
n_agents, n_reps = 25_000, 190
# (label, setting, dt, transmission, exp 10 arm with the same seeds or None): weekly first, so the longest runs start first
arms = [
    ('4: B, weekly, mass action linear', 'B', '1 week', 'linear', None),
    ('1: B, monthly, mass action linear', 'B', '1 month', 'linear', '5: B, mass action linear'),
    ('2: B, monthly, network', 'B', '1 month', 'network', '4: B, network'),
    ('3: A, monthly, mass action linear', 'A', '1 month', 'linear', '3: A, mass action linear'),
]
outputs = ['prev_active', 'incidence_py', 'incidence_event_py', 'INFECTION', 'tb_deaths_py']
scale = dict(prev_active=1e5, incidence_py=1e5, incidence_event_py=1e5, tb_deaths_py=1e5)


def run_abm(setting, dt, transmission, seed):
    pars = settings()[setting]
    ode = tbode.TBODE(**pars)
    x = ode.equilibrium()
    if transmission == 'network':
        kw = dict(tb_class=tbode.TBExact)
    else:
        kw = dict(tb_class=tbode.TBExactMassAction, network=False, beta_ode=pars['beta'], infection_prob=transmission)
    sim = tbode.make_abm(ode.pars, seed=seed, n_agents=n_agents, start=start, stop=stop, dt=dts[dt], infection_timing='within', **kw)
    sim.init()
    tbode.init_from_ode(sim, x, seed=seed, ode=ode)
    sim.run()
    return tbode.annualize(tbode.abm_outputs(sim))


def last50(oo, k):
    """Each rep's mean over the last 50 years, in display units."""
    return np.array([np.nanmean(o[k][o.time >= stop - 50]) for o in oo]) * scale.get(k, 1)


if __name__ == '__main__':
    T = sc.timer()
    iterkwargs = [dict(setting=a[1], dt=a[2], transmission=a[3], seed=s) for a in arms for s in range(n_reps)]
    outs = sc.parallelize(run_abm, iterkwargs=iterkwargs, ncpus=tbode.ncpus)
    abms = {a[0]: outs[i * n_reps:(i + 1) * n_reps] for i, a in enumerate(arms)}
    sc.save(os.path.join(outdir, 'outputs.obj'), dict(abms=abms, arms=arms))
    T.toc('ABM runs')

    # The counting must not change the dynamics: every rep's prevalence equals exp 10's with the same seed
    out10 = sc.load(os.path.join(dir10, 'outputs.obj'))['abms']
    for lab, _, _, _, ref in arms:
        if ref:
            same = [np.array_equal(a.prev_active, b.prev_active, equal_nan=True) for a, b in zip(abms[lab], out10[ref])]
            assert all(same), f'{lab}: {len(same) - sum(same)} reps differ from exp 10'
            print(f'Check passed: {lab}: all {len(same)} reps bit-identical to exp 10 ({ref})')

    sets = settings()
    rows = []
    for lab, setting, dt, transmission, _ in arms:
        ode = tbode.TBODE(**sets[setting])
        eq = ode.outputs(ode.equilibrium())
        oo = abms[lab]
        vals = {k: last50(oo, k) for k in outputs}
        for k in outputs:
            ode_key = 'incidence_py' if k == 'incidence_event_py' else k
            ode_val = (tbode.ode_outputs(ode.outputs(ode.equilibrium(), time=np.zeros(1)))[ode_key][0]) * scale.get(k, 1)
            v = vals[k] / ode_val
            rows.append(dict(arm=lab, setting=setting, dt=dt, transmission=transmission, output=k, n_reps=len(v),
                             ratio_mean=v.mean(), ratio_se=v.std(ddof=1) / np.sqrt(len(v))))
        r = vals['incidence_event_py'] / vals['incidence_py']  # Paired, same runs
        rows.append(dict(arm=lab, setting=setting, dt=dt, transmission=transmission, output='event / snapshot (paired)', n_reps=len(r),
                         ratio_mean=r.mean(), ratio_se=r.std(ddof=1) / np.sqrt(len(r))))
        r = (vals['incidence_event_py'] / vals['prev_active']) / (eq.incidence_py[0] / eq.prev_active[0])  # Incidence per prevalent case, ABM / ODE
        rows.append(dict(arm=lab, setting=setting, dt=dt, transmission=transmission, output='event incidence per prevalent case, ABM / ODE', n_reps=len(r),
                         ratio_mean=r.mean(), ratio_se=r.std(ddof=1) / np.sqrt(len(r))))
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(outdir, 'summary.csv'), index=False)

    fig, axs = plt.subplots(1, len(arms), figsize=(18, 4.5), sharey=False)
    for ax, (lab, setting, dt, transmission, _) in zip(axs, arms):
        ode = tbode.TBODE(**sets[setting])
        ode_inc = ode.outputs(ode.equilibrium()).incidence_py[0] * 1e5
        oo = abms[lab]
        for k, col, name in [('incidence_py', 'C3', 'Snapshot count'), ('incidence_event_py', 'C0', 'Event count (transitions)')]:
            yr = np.floor(oo[0].time + 1e-9)
            years = np.unique(yr)
            Y = np.vstack([[np.nanmean(np.asarray(o[k])[yr == y]) for y in years] for o in oo]) * 1e5
            kern = np.ones(5) / 5
            mean = np.convolve(Y.mean(0), kern, mode='same') / np.convolve(np.ones(len(years)), kern, mode='same')
            ax.plot(years, mean, color=col, lw=1.8, label=f'{name}: {last50(oo, k).mean():.0f}')
        ax.axhline(ode_inc, color='k', ls=':', label=f'ODE: {ode_inc:.0f}')
        ax.set_ylim(0, ode_inc * 1.2)
        ax.set_title(lab, fontsize=10)
        ax.set_xlabel('Year')
        ax.legend(fontsize=8, frameon=False, loc='lower left')
    axs[0].set_ylabel('Incidence per 100k per year')
    fig.suptitle(f'TBExact incidence: snapshot vs event counting, mean of {n_reps} reps (5-year moving average); legend: mean 1950–2000', fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, 'incidence_counting.png'), dpi=110)

    with pd.option_context('display.width', 250, 'display.float_format', '{:.4f}'.format):
        print(df.to_string(index=False))
    T.toc('Total')
