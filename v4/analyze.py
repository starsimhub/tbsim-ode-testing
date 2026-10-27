"""
Figures and summary tables for v4_findings.qmd, from the saved results (no runs)

Agents ÷ ODE is the mean over 1950–2000 of each run's yearly values, divided by the ODE's mean over the same years; the
tables give the mean ± SE over runs.
"""
import numpy as np
import pandas as pd
import sciris as sc
import matplotlib.pyplot as plt
import common as c
from run_growth import crossing

res = c.folder / 'results'
figs = c.folder / 'figs'
colors = sc.objdict(day='#2a78d6', week='#eb6834', month='#1baf7a', net='#6250d6', net_matched='#e87ba4') # Fixed categorical slots
keys = ['prevalence', 'infection', 'non_infectious', 'cleared', 'incidence', 'tb_deaths']
window = (1950, 2000)


def ratios(data):
    """ Agents ÷ ODE for each arm and output: mean, SE, and the 10th and 90th percentiles over runs """
    rows = []
    for key, runs in data.abm.items():
        arm = key.split('_')[0]
        ode = data.ode[arm]
        sel = lambda o: (o.year >= window[0]) & (o.year < window[1])
        m = lambda o, out: o[out][sel(o)].mean()
        for out in keys:
            vals = np.array([m(r, out) for r in runs])/m(ode, out)
            rows.append(dict(arm=key, output=out, reps=len(runs), mean=vals.mean(), se=vals.std(ddof=1)/np.sqrt(len(vals)), p10=np.percentile(vals, 10), p90=np.percentile(vals, 90)))
        # Incidence per prevalent case: whether incidence matches once prevalence is allowed for
        vals = np.array([m(r, 'incidence')/m(r, 'prevalence') for r in runs])/(m(ode, 'incidence')/m(ode, 'prevalence'))
        rows.append(dict(arm=key, output='incidence_per_case', reps=len(runs), mean=vals.mean(), se=vals.std(ddof=1)/np.sqrt(len(vals)), p10=np.percentile(vals, 10), p90=np.percentile(vals, 90)))
    return pd.DataFrame(rows)


def table(df, outputs, arms):
    """ A Markdown table of agents ÷ ODE (mean ± SE) by arm and output """
    head = '| Arm | Reps | ' + ' | '.join(outputs) + ' |\n|' + '---|'*(len(outputs) + 2) + '\n'
    lines = []
    for arm in arms:
        sub = df[df.arm == arm].set_index('output')
        lines.append(f'| {arm} | {sub.reps.iloc[0]} | ' + ' | '.join(f'{sub.loc[o, "mean"]:.3f} ± {sub.loc[o, "se"]:.3f}' for o in outputs) + ' |')
    return head + '\n'.join(lines)


def plot_runs(ax, data, arm, out, label_of, scale=1.0):
    ode = data.ode[arm]
    for key, runs in data.abm.items():
        if key.split('_')[0] != arm or key not in label_of:
            continue
        col = colors[key.split('_', 2)[2] if key.count('_') >= 2 else key.split('_')[1]]
        mat = np.array([r[out] for r in runs])*scale
        ax.fill_between(runs[0].year, *np.percentile(mat, [10, 90], axis=0), color=col, alpha=0.15, lw=0)
        ax.plot(runs[0].year, mat.mean(axis=0), color=col, lw=2, label=f'{label_of[key]} ({len(runs)} reps)')
    ax.plot(ode.year, ode[out]*scale, 'k--', lw=2, label='ODE')
    ax.set_ylim(bottom=0)
    ax.set_xlabel('Year')


def fig_equilibrium(data, arms_labels, fn, title):
    fig, axs = plt.subplots(len(arms_labels), 2, figsize=(11, 3.6*len(arms_labels)), squeeze=False)
    for row, (arm, label_of) in zip(axs, arms_labels.items()):
        plot_runs(row[0], data, arm, 'infection', label_of)
        row[0].set_title(f'{arm}: latent infection (share)')
        plot_runs(row[1], data, arm, 'prevalence', label_of, scale=1e5)
        row[1].set_title(f'{arm}: active prevalence (per 100k)')
        row[1].legend(frameon=False, fontsize=8, loc='lower left')
    fig.suptitle(title + '\nlines: mean of runs; bands: 10–90% of runs (yearly means)')
    sc.figlayout()
    fig.savefig(figs / fn, dpi=130)
    plt.close(fig)


def fig_incidence(data, label_of, fn):
    arms = sorted({k.split('_')[0] for k in label_of})
    fig, axs = plt.subplots(1, len(arms), figsize=(5*len(arms), 3.6), squeeze=False)
    for ax, arm in zip(axs[0], arms):
        plot_runs(ax, data, arm, 'incidence', {k: v for k, v in label_of.items() if k.startswith(arm + '_')}, scale=1e5)
        ax.set_title(f'{arm}: incidence (per 100k per year)')
        ax.legend(frameon=False, fontsize=8, loc='lower left')
    sc.figlayout()
    fig.savefig(figs / fn, dpi=130)
    plt.close(fig)


def growth_summary(data):
    rows = []
    for key, runs in data.abm.items():
        arm, dt = key.split('_')
        for level in [0.0025, 0.004]:
            ode = crossing(data.ode[arm], level)
            t = np.array([crossing(r, level) for r in runs])
            t = t[np.isfinite(t)]
            rows.append(dict(arm=arm, dt=dt, level=f'{1e5*level:.0f}/100k', ode=ode, agents=t.mean(), lag=t.mean() - ode, se=t.std(ddof=1)/np.sqrt(len(t)), reps=len(t)))
    return pd.DataFrame(rows)


def fig_growth(data, fn):
    arms = sorted({k.split('_')[0] for k in data.abm})
    fig, axs = plt.subplots(1, len(arms), figsize=(5.5*len(arms), 3.6), squeeze=False)
    for ax, arm in zip(axs[0], arms):
        for key, runs in data.abm.items():
            if not key.startswith(arm + '_'):
                continue
            dt = key.split('_')[1]
            mat = np.array([r.prevalence for r in runs])*1e5
            ax.plot(runs[0].year, mat.mean(axis=0), color=colors[dt], lw=2, label=f'Agents, 1 {dt} ({len(runs)} reps)')
        ax.plot(data.ode[arm].year, data.ode[arm].prevalence*1e5, 'k--', lw=2, label='ODE')
        ax.set_title(f'{arm}: active prevalence (per 100k), from a 0.1% seed')
        ax.set_xlabel('Year')
        ax.set_ylim(bottom=0)
        ax.legend(frameon=False, fontsize=8)
    sc.figlayout()
    fig.savefig(figs / fn, dpi=130)
    plt.close(fig)


if __name__ == '__main__':
    figs.mkdir(exist_ok=True)
    md = []
    base = sc.load(res / 'baseline.obj')
    prog = sc.load(res / 'progression.obj')
    df = pd.concat([ratios(base), ratios(prog)])
    df.to_csv(res / 'summary.csv', index=False)
    outs = ['prevalence', 'infection', 'cleared', 'incidence', 'incidence_per_case', 'tb_deaths']
    md.append('## Baseline\n\n' + table(df, outs, list(base.abm.keys())))
    md.append('## Progression\n\n' + table(df, outs, list(prog.abm.keys())))

    fig_equilibrium(base, sc.objdict(
        A0 = {'A0_day': 'Agents, 1 day', 'A0_week': 'Agents, 1 week', 'A0_month': 'Agents, 1 month'},
        B0 = {'B0_day': 'Agents, 1 day', 'B0_week': 'Agents, 1 week', 'B0_month': 'Agents, 1 month'},
    ), 'baseline.png', 'Baseline (constant progression), started at the ODE equilibrium: agents vs ODE')
    fig_equilibrium(base, sc.objdict(
        A0 = {'A0_month': 'Agents, 1 month, homogeneous', 'A0_month_net': 'Network, rate β/K', 'A0_month_net_matched': 'Network, matched per-step probability'},
    ), 'network.png', 'A0 at a monthly step: homogeneous mixing vs a random network (K = 10, new contacts every step)')
    fig_equilibrium(prog, sc.objdict(
        Ak = {'Ak_week': 'Agents, 1 week', 'Ak_month': 'Agents, 1 month'},
        Bk = {'Bk_week': 'Agents, 1 week', 'Bk_month': 'Agents, 1 month'},
        Bremote = {'Bremote_week': 'Agents, 1 week', 'Bremote_month': 'Agents, 1 month'},
    ), 'progression.png', 'Progression by time since infection: agents vs ODE')
    fig_incidence(base, {'A0_week': 'Agents, 1 week', 'A0_month': 'Agents, 1 month', 'B0_week': 'Agents, 1 week', 'B0_month': 'Agents, 1 month'}, 'incidence.png')

    if (res / 'growth.obj').exists():
        growth = sc.load(res / 'growth.obj')
        gs = growth_summary(growth)
        gs.to_csv(res / 'growth.csv', index=False)
        md.append('## Growth\n\n' + gs.to_markdown(index=False, floatfmt='.2f'))
        fig_growth(growth, 'growth.png')

    one = pd.read_csv(res / 'onestep.csv')
    one['diff'] = one.abm - one.exact
    one['ode_diff'] = one.ode - one.exact
    worst = one.loc[one.groupby(['parset', 'dt', 'start'])['diff'].apply(lambda d: d.abs().idxmax())]
    md.append('## One step (largest agents − exact per start state)\n\n' + worst[['parset', 'dt', 'start', 'state', 'exact', 'abm', 'diff', 'z', 'ode_diff']].to_markdown(index=False, floatfmt='.5f'))
    if (res / 'timing.csv').exists():
        md.append('## Timing\n\n' + pd.read_csv(res / 'timing.csv').to_markdown(index=False, floatfmt='.3f'))
    text = '\n\n'.join(md)
    (res / 'summary.md').write_text(text)
    print(text)
