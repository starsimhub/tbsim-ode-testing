"""
Experiment 17, v4 environment: summaries and figures from the saved results. See README.md.

    v4/.venv/bin/python experiments/17_latent_reinfection/compare.py q1    # question 1: agents ÷ ODE (q1.obj -> q1_agents_vs_ode.png)
    v4/.venv/bin/python experiments/17_latent_reinfection/compare.py q1n   # by population size (q1_n_summary.csv -> q1_n_population_size.png)
    v4/.venv/bin/python experiments/17_latent_reinfection/compare.py q2    # question 2 (lai_tpt.obj, q2_v4.obj -> q2_summary.csv, q2_*.png)
    v4/.venv/bin/python experiments/17_latent_reinfection/compare.py q2 plot   # question 2's figures only, from q2_summary.obj
    v4/.venv/bin/python experiments/17_latent_reinfection/compare.py latent    # latent spell durations (-> q2_latent_spells.obj, .png); add plot to redraw only
"""
import sys
import numpy as np
import pandas as pd
import sciris as sc
import matplotlib.pyplot as plt

outdir = sc.thispath(__file__)
sys.path.insert(0, str(outdir))
import run_v4 as rv # noqa: E402

outputs_q1 = sc.objdict(prevalence='Active prevalence', latent='Latent (recent + remote)', recent='Recent latent', remote='Remote latent',
                        cleared='Cleared', recovered='Recovered', tb_deaths='TB deaths')
colors = sc.objdict(month='#3b6ea8', week='#c4620b')


def plot_q1():
    """ Agents ÷ ODE over 1950–2000 (mean ± 2 SE, and each rep), by output, phi and dt """
    data = sc.load(rv.fn_q1)
    fig, axes = plt.subplots(1, len(outputs_q1), figsize=(2.1*len(outputs_q1), 3.4), sharey=True)
    for ax, (key, label) in zip(axes, outputs_q1.items()):
        for j, dt in enumerate(rv.dts):
            for i, phi in enumerate(rv.phis):
                akey = f'phi{phi}_{dt}'
                if akey not in data.abm:
                    continue
                r = rv.ratio(data.abm[akey], data.ode[akey], key)
                x = i + (j - 0.5)*0.3
                ax.plot(np.full(len(r), x), r, '.', color=colors[dt], alpha=0.25, ms=3)
                ax.errorbar(x, r.mean(), yerr=2*r.std(ddof=1)/np.sqrt(len(r)), fmt='o', color=colors[dt], ms=4, capsize=2,
                            label=f'{dt}ly dt' if (i == 0 and key == 'prevalence') else None)
        ax.axhline(1, color='0.4', lw=0.8)
        ax.set_xticks(range(len(rv.phis)), [str(p) for p in rv.phis])
        ax.set_xlabel('φ')
        ax.set_title(label, fontsize=9)
    axes[0].set_ylabel('Agents ÷ ODE, 1950–2000')
    axes[0].legend(fontsize=8, frameon=False)
    fig.suptitle('Exp 17, question 1: v4 agents ÷ v4 ODE with latent reinfection (set B; mean ± 2 SE, dots = reps)', fontsize=10)
    fig.tight_layout()
    fig.savefig(outdir / 'q1_agents_vs_ode.png', dpi=160)
    print('Saved q1_agents_vs_ode.png')


def plot_q1n():
    """ Agents ÷ ODE over 1950–2000 by population size (phi = 0.01, monthly dt; mean ± 2 SE), from q1_n_summary.csv """
    import pandas as pd
    df = pd.read_csv(outdir / 'q1_n_summary.csv')
    outs = sc.objdict(prevalence='Active prevalence', tb_deaths='TB deaths', susceptible='Susceptible', cleared='Cleared')
    fig, axes = plt.subplots(1, len(outs), figsize=(2.6*len(outs), 3.2))
    for ax, (key, label) in zip(axes, outs.items()):
        d = df[df.output == key]
        ax.errorbar(d.n_agents, d['mean'], yerr=2*d.se, fmt='o', color=colors.month, ms=4, capsize=2)
        ax.axhline(1, color='0.4', lw=0.8)
        ax.set_xscale('log')
        ax.set_xticks(d.n_agents, [f'{n/1000:g}k' for n in d.n_agents], fontsize=8)
        ax.minorticks_off()
        ax.set_xlabel('Agents')
        ax.set_title(label, fontsize=9)
    axes[0].set_ylabel('Agents ÷ ODE, 1950–2000')
    fig.suptitle('Exp 17, question 1 by population size: set B, φ = 0.01, monthly dt (mean ± 2 SE; 40–160 reps)', fontsize=10)
    fig.tight_layout()
    fig.savefig(outdir / 'q1_n_population_size.png', dpi=160)
    print('Saved q1_n_population_size.png')


# Question 2
fn_q2_lai = outdir / 'lai_tpt.obj'
fn_q2_sum = outdir / 'q2_summary.obj'
models_q2 = ['lai_tpt'] + [f'phi{phi}' for phi in rv.phis] # tb_LAI_TPT as run, then the two-state model at each phi
label_q2 = dict(lai_tpt='tb_LAI_TPT as run (v0.10, monthly dt)', **{f'phi{phi}': f'Two-state, φ = {phi}' for phi in rv.phis})


def run_stats(o):
    """ Scalar outputs of one cohort run (cohort.Tracker results); times in years """
    s = sc.objdict()
    v010 = o.t[0] > 0 # v0.10's first point is the end of month 1, so time 0 is added
    t = np.concatenate([[0], o.t]) if v010 else o.t
    pad = lambda x, x0: np.concatenate([[x0], x]) if v010 else x
    for key, x0 in [('A', 0.0), ('B', 0.0), ('C', 0.0), ('latent', 1.0)]:
        for y in [1, 10, 30]:
            s[f'{key}{y}'] = float(np.interp(y, t, pad(o[key], x0)))
    if 'B_fp' in o: # tb_LAI_TPT's FirstPassage: end-of-step snapshots
        s.B1_snapshot, s.B10_snapshot, s.B30_snapshot = (float(np.interp(y, t, pad(o.B_fp, 0.0))) for y in [1, 10, 30])
    n = o.n
    k = np.arange(len(o.resets_pre))
    s.resets_per_person = float(k @ o.resets_pre/n)
    s.share_reset = float(o.resets_pre[1:].sum()/n)
    sp = o.split
    s.prog_if_reset = sp.n_reset_prog/sp.n_reset if sp.n_reset else np.nan
    s.prog_if_none = sp.n_none_prog/sp.n_none
    h = o.hists.reset_to_first_prog
    steps = np.arange(len(h))
    s.reset_to_first_prog_mean = float(steps @ h/h.sum()*o.dt_year)
    s.reset_to_first_prog_within_1y = float(h[:13].sum()/h.sum())
    h = o.hists.spell_prog
    s.spell_prog_mean = float(np.arange(len(h)) @ h/h.sum()*o.dt_year)
    h = o.hists.spell_clear
    s.spell_clear_mean = float(np.arange(len(h)) @ h/h.sum()*o.dt_year)
    e = o.episodes
    s.episodes_per_person = float(np.arange(len(e)) @ e/n)
    s.share_2plus_episodes = float(e[2:].sum()/n)
    return s


def wquantile(x, w, q):
    """ Weighted quantile (linear in the cumulative weight at each value's midpoint) """
    i = np.argsort(x)
    x, w = np.asarray(x)[i], np.asarray(w)[i]
    cw = (np.cumsum(w) - 0.5*w)/w.sum()
    return float(np.interp(q, cw, x))


def summarise_q2():
    """ Per run, per set (mean and SE over seeds; on - off paired by seed) and across sets (weighted); chains alongside """
    lai = sc.load(fn_q2_lai)
    v4 = sc.load(rv.fn_q2)
    w = np.asarray(v4.weights)
    assert np.allclose(w, lai.weights)
    runs = {}
    for key, o in lai.runs.items(): # p{i}_lam{lam}_{on|off}_s{seed}
        p, lam, arm, seed = key.split('_')
        runs[('lai_tpt', int(p[1:]), float(lam[3:]), arm, int(seed[1:]))] = run_stats(o)
    for key, o in v4.runs.items(): # p{i}_phi{phi}_lam{lam}_{on|off}_s{seed}
        p, phi, lam, arm, seed = key.split('_')
        runs[(phi, int(p[1:]), float(lam[3:]), arm, int(seed[1:]))] = run_stats(o)
    stat_keys = list(next(iter(runs.values())).keys())
    seeds = sorted({k[4] for k in runs})
    rows = []
    for model in models_q2:
        for lam, on in rv.q2_arms():
            arms = [('on' if on else 'off', False)] + ([('on-off', True)] if lam > 0 and on else [])
            for arm, diff in arms:
                for stat in stat_keys:
                    per_set, se_set = [], []
                    for i in range(rv.n_sets):
                        get = lambda a, s: runs[(model, i, lam, a, s)].get(stat, np.nan)
                        vals = np.array([get('on', s) - get('off', s) if diff else get(arm, s) for s in seeds])
                        per_set.append(vals.mean())
                        se_set.append(vals.std(ddof=1)/np.sqrt(len(vals)))
                    per_set, se_set = np.array(per_set), np.array(se_set)
                    if np.all(np.isnan(per_set)):
                        continue
                    rows.append(dict(model=model, lam=lam, arm=arm, stat=stat, mean=float(w @ per_set), se=float(np.sqrt(w**2 @ se_set**2)),
                                     p25=wquantile(per_set, w, 0.25), p75=wquantile(per_set, w, 0.75), per_set=per_set, se_set=se_set))
    # Chains (keys p{i}_{phi... | exact}_lam{lam}_{on|off}): exact curves at 1, 10, 30 years (on, off, and on - off)
    chain_rows = []
    for model in ['exact'] + [f'phi{phi}' for phi in rv.phis]:
        for lam, on in rv.q2_arms():
            arms = ['on' if on else 'off'] + (['on-off'] if lam > 0 and on else [])
            for arm in arms:
                for curve in ['A', 'B', 'C', 'latent']:
                    for y in [1, 10, 30]:
                        def val(i, a):
                            ch = v4.chains[f'p{i}_{model}_lam{lam}_{a}']
                            return float(np.interp(y, ch.t, ch[curve]))
                        per_set = np.array([val(i, 'on') - val(i, 'off') if arm == 'on-off' else val(i, arm) for i in range(rv.n_sets)])
                        chain_rows.append(dict(model=model, lam=lam, arm=arm, stat=f'{curve}{y}', mean=float(w @ per_set),
                                               p25=wquantile(per_set, w, 0.25), p75=wquantile(per_set, w, 0.75), per_set=per_set))
    agents = pd.DataFrame(rows)
    chains = pd.DataFrame(chain_rows)
    agents.drop(columns=['per_set', 'se_set']).to_csv(outdir / 'q2_summary.csv', index=False)
    chains.drop(columns=['per_set']).to_csv(outdir / 'q2_chains.csv', index=False)
    curves = q2_mean_curves(lai, v4, w)
    hists = q2_mean_hists(lai, v4, w)
    out = sc.objdict(agents=agents, chains=chains, curves=curves, hists=hists, weights=w)
    sc.save(fn_q2_sum, out)
    return out


def q2_mean_curves(lai, v4, w):
    """ Weighted mean curves (over sets, after averaging seeds) on a monthly grid from 0 to 30 years, for agents and chains """
    grid = np.arange(0, 30*12 + 1)/12
    curves = sc.objdict()
    def interp(o, key):
        t = np.concatenate([[0], o.t]) if o.t[0] > 0 else o.t
        x = np.concatenate([[1.0 if key == 'latent' else 0.0], o[key]]) if o.t[0] > 0 else o[key]
        return np.interp(grid, t, x)
    for lam, on in rv.q2_arms():
        a = 'on' if on else 'off'
        for key in ['A', 'B', 'C', 'latent']:
            per_set = [np.mean([interp(lai.runs[f'p{i}_lam{lam}_{a}_s{s}'], key) for s in rv.seeds_q2], axis=0) for i in range(rv.n_sets)]
            curves[f'lai_tpt_lam{lam}_{a}_{key}'] = w @ np.array(per_set)
            per_set = [np.interp(grid, v4.chains[f'p{i}_exact_lam{lam}_{a}'].t, v4.chains[f'p{i}_exact_lam{lam}_{a}'][key]) for i in range(rv.n_sets)]
            curves[f'exact_lam{lam}_{a}_{key}'] = w @ np.array(per_set)
            for phi in rv.phis:
                per_set = [np.mean([interp(v4.runs[f'p{i}_phi{phi}_lam{lam}_{a}_s{s}'], key) for s in rv.seeds_q2], axis=0) for i in range(rv.n_sets)]
                curves[f'phi{phi}_lam{lam}_{a}_{key}'] = w @ np.array(per_set)
                per_set = [np.interp(grid, v4.chains[f'p{i}_phi{phi}_lam{lam}_{a}'].t, v4.chains[f'p{i}_phi{phi}_lam{lam}_{a}'][key]) for i in range(rv.n_sets)]
                curves[f'chain_phi{phi}_lam{lam}_{a}_{key}'] = w @ np.array(per_set)
    curves.t = grid
    return curves


def q2_mean_hists(lai, v4, w):
    """ Weighted mean distributions (each set's seeds pooled, normalised, then weighted), reinfection on: per-person outputs """
    hists = sc.objdict()
    def norm(h):
        return h/h.sum() if h.sum() else h
    def pooled(get, n_sets, n):
        out = np.zeros(n)
        for i in range(n_sets):
            h = np.zeros(n)
            for s in rv.seeds_q2:
                x = get(i, s)
                h[:min(n, len(x))] += x[:n]
            out += w[i]*norm(h)
        return out
    n = 30*12 + 2
    for lam in rv.lams:
        for name in ['reset_to_first_prog', 'spell_prog', 'spell_clear']:
            hists[f'lai_tpt_lam{lam}_{name}'] = pooled(lambda i, s: lai.runs[f'p{i}_lam{lam}_on_s{s}'].hists[name], rv.n_sets, n)
            for phi in rv.phis:
                hists[f'phi{phi}_lam{lam}_{name}'] = pooled(lambda i, s: v4.runs[f'p{i}_phi{phi}_lam{lam}_on_s{s}'].hists[name], rv.n_sets, n)
        for name in ['resets_pre', 'episodes']:
            m = 21
            hists[f'lai_tpt_lam{lam}_{name}'] = pooled(lambda i, s: lai.runs[f'p{i}_lam{lam}_on_s{s}'][name], rv.n_sets, m)
            for phi in rv.phis:
                hists[f'phi{phi}_lam{lam}_{name}'] = pooled(lambda i, s: v4.runs[f'p{i}_phi{phi}_lam{lam}_on_s{s}'][name], rv.n_sets, m)
    hists.dt_year = 1/12
    return hists


def print_q2(out):
    """ The headline tables """
    a, ch = out.agents, out.chains
    def get(df, model, lam, arm, stat, col='mean'):
        r = df[(df.model == model) & np.isclose(df.lam, lam) & (df.arm == arm) & (df.stat == stat)]
        return float(r[col].iloc[0]) if len(r) else np.nan
    print('\nWhat latent reinfection adds (on - off), share of the cohort; agents: weighted mean ± SE [25-75% across sets]; chain in brackets (for lai_tpt: its model solved exactly)')
    for lam in [0.10, 0.60]:
        for stat in ['A10', 'A30', 'B10', 'B30']:
            line = f'  lam {lam:.2f} {stat:4s}'
            for model in models_q2:
                cm = 'exact' if model == 'lai_tpt' else model
                line += (f' | {model}: {get(a, model, lam, "on-off", stat):.5f} ± {get(a, model, lam, "on-off", stat, "se"):.5f} '
                         f'[{get(a, model, lam, "on-off", stat, "p25"):.4f}, {get(a, model, lam, "on-off", stat, "p75"):.4f}] (chain {get(ch, cm, lam, "on-off", stat):.5f})')
            print(line)
    print('\nWith reinfection on: agents (weighted mean ± SE); two-state ÷ tb_LAI_TPT')
    for lam in rv.lams:
        for stat in ['A1', 'A10', 'A30', 'B1', 'B10', 'B30', 'C10', 'C30', 'latent1', 'latent10', 'latent30']:
            ref = get(a, 'lai_tpt', lam, 'on', stat)
            line = f'  lam {lam:.2f} {stat:9s} tb_LAI_TPT {ref:.4f} ± {get(a, "lai_tpt", lam, "on", stat, "se"):.4f} (solved exactly {get(ch, "exact", lam, "on", stat):.4f})'
            for phi in rv.phis:
                v = get(a, f'phi{phi}', lam, 'on', stat)
                line += f' | φ {phi}: {v:.4f} ({v/ref:.3f})' if ref > 0 else f' | φ {phi}: {v:.4f}'
            print(line)
    print('\nPer person, reinfection on (weighted mean ± SE)')
    for lam in rv.lams:
        for stat in ['resets_per_person', 'share_reset', 'prog_if_reset', 'prog_if_none', 'reset_to_first_prog_mean', 'reset_to_first_prog_within_1y',
                     'spell_prog_mean', 'spell_clear_mean', 'episodes_per_person', 'share_2plus_episodes']:
            line = f'  lam {lam:.2f} {stat:30s}'
            for model in models_q2:
                line += f' | {model}: {get(a, model, lam, "on", stat):.4f} ± {get(a, model, lam, "on", stat, "se"):.4f}'
            print(line)
    print('\nChain check: v4 agents ÷ two-state chain (weighted mean over sets), reinfection on and off')
    for phi in rv.phis:
        for lam, on in rv.q2_arms():
            arm = 'on' if on else 'off'
            vals = [get(a, f'phi{phi}', lam, arm, s)/get(ch, f'phi{phi}', lam, arm, s) for s in ['A1', 'A10', 'A30', 'B1', 'B10', 'B30', 'C10', 'C30', 'latent1', 'latent10']]
            print(f'  φ {phi} lam {lam:.2f} {arm:3s}: ' + ' '.join(f'{v:.3f}' for v in vals) + '  (A1 A10 A30 B1 B10 B30 C10 C30 latent1 latent10)')
    print('\ntb_LAI_TPT as run, lam = 0: B counted from transitions vs end-of-step snapshots (tb_LAI_TPT\'s FirstPassage)')
    for stat in ['B1', 'B10']:
        print(f'  {stat}: transitions {get(a, "lai_tpt", 0.0, "on", stat):.4f}, snapshots {get(a, "lai_tpt", 0.0, "on", stat + "_snapshot"):.4f}, solved exactly {get(ch, "exact", 0.0, "on", stat):.4f}')
    print('  What latent reinfection adds, by snapshots (on - off):')
    for lam in [0.10, 0.60]:
        print(f'    lam {lam:.2f}: ' + ', '.join(f'{s} {get(a, "lai_tpt", lam, "on-off", s + "_snapshot"):.5f} ± {get(a, "lai_tpt", lam, "on-off", s + "_snapshot", "se"):.5f}' for s in ['B10', 'B30']))
    return


colors_q2 = sc.objdict(lai_tpt='0.15', **{f'phi{phi}': c for phi, c in zip(rv.phis, ['#3b6ea8', '#c4620b', '#2e8b57'])})
colors_lam = {0.0: '0.5', 0.10: '#3b6ea8', 0.60: '#c4620b'}


def plot_q2(out):
    """ Question 2's figures: what latent reinfection adds; cumulative curves; per-person distributions """
    a, ch = out.agents, out.chains
    def row(df, model, lam, arm, stat):
        r = df[(df.model == model) & np.isclose(df.lam, lam) & (df.arm == arm) & (df.stat == stat)]
        return r.iloc[0] if len(r) else None

    # 1. What latent reinfection adds, against phi
    stats = sc.objdict(B10='Ever asymptomatic TB by 10 years', B30='... by 30 years', A30='Ever progressed by 30 years')
    fig, axes = plt.subplots(2, len(stats), figsize=(4*len(stats), 6.4))
    for r_, lam in enumerate([0.10, 0.60]):
        for ax, (stat, label) in zip(axes[r_], stats.items()):
            x = np.arange(len(rv.phis))
            ag = [row(a, f'phi{phi}', lam, 'on-off', stat) for phi in rv.phis]
            ax.errorbar(x, [g['mean'] for g in ag], yerr=[2*g['se'] for g in ag], fmt='o', color='#3b6ea8', capsize=3, label='Two-state, v4 agents')
            ax.plot(x, [row(ch, f'phi{phi}', lam, 'on-off', stat)['mean'] for phi in rv.phis], '-', color='#3b6ea8', lw=1, label='Two-state, chain')
            ref = row(a, 'lai_tpt', lam, 'on-off', stat)
            ax.axhspan(ref['mean'] - 2*ref['se'], ref['mean'] + 2*ref['se'], color='0.8', alpha=0.6)
            ax.axhline(ref['mean'], color='0.15', lw=1.2, label='tb_LAI_TPT as run (v0.10, monthly dt)')
            ax.axhline(row(ch, 'exact', lam, 'on-off', stat)['mean'], color='0.15', lw=1, ls='--', label="tb_LAI_TPT's model, solved exactly (continuous time)")
            ax.set_xticks(x, [str(p) for p in rv.phis])
            ax.set_xlabel('φ')
            ax.set_ylim(bottom=0)
            ax.set_title(f'{label}, λ = {lam:.2f}/y', fontsize=9)
        axes[r_][0].set_ylabel('Added by latent reinfection\n(on − off, share of cohort)')
    axes[0][0].legend(fontsize=7, frameon=False, loc='lower left')
    fig.suptitle('Exp 17, question 2: extra progression from reinfecting latent people (12 sets, weighted mean ± 2 SE; band: tb_LAI_TPT ± 2 SE)', fontsize=10)
    fig.tight_layout()
    fig.savefig(outdir / 'q2_reinfection_adds.png', dpi=160)

    # 2. Cumulative curves with reinfection on, and the latent share
    c = out.curves
    fig, axes = plt.subplots(2, len(rv.phis), figsize=(4*len(rv.phis), 6.4), sharey='row')
    for j, phi in enumerate(rv.phis):
        for lam in rv.lams:
            for r_, key in enumerate(['B', 'latent']):
                ax = axes[r_][j]
                ax.plot(c.t, c[f'phi{phi}_lam{lam}_on_{key}'], color=colors_lam[lam], lw=1.5, label=f'Two-state, λ = {lam:.2f}')
                ax.plot(c.t, c[f'lai_tpt_lam{lam}_on_{key}'], color=colors_lam[lam], lw=1.2, ls='--', label=f'tb_LAI_TPT as run (v0.10, monthly dt), λ = {lam:.2f}')
        axes[0][j].set_title(f'φ = {phi}', fontsize=9)
        axes[1][j].set_xlabel('Years since infection')
    axes[0][0].set_ylabel('Ever asymptomatic TB')
    axes[1][0].set_ylabel('Latent')
    axes[0][0].legend(fontsize=7, frameon=False)
    fig.suptitle('Exp 17, question 2: cohorts with latent reinfection on (12 sets, weighted mean; agents)', fontsize=10)
    fig.tight_layout()
    fig.savefig(outdir / 'q2_curves.png', dpi=160)

    # 3. Per-person distributions, reinfection on
    h = out.hists
    yrs = np.arange(30*12 + 2)*h.dt_year
    fig, axes = plt.subplots(2, 4, figsize=(15, 6.4))
    for r_, lam in enumerate([0.10, 0.60]):
        for model in models_q2:
            kw = dict(color=colors_q2[model], lw=1.8 if model == 'lai_tpt' else 1.2, ls='--' if model == 'lai_tpt' else '-', label=label_q2[model])
            axes[r_][0].plot(yrs, np.cumsum(h[f'{model}_lam{lam}_reset_to_first_prog']), **kw)
            axes[r_][1].plot(yrs, np.cumsum(h[f'{model}_lam{lam}_spell_prog']), **kw)
            off = models_q2.index(model)*0.2 - 0.3
            for col, name in [(2, 'resets_pre'), (3, 'episodes')]:
                y = h[f'{model}_lam{lam}_{name}'][:6]
                axes[r_][col].bar(np.arange(6) + off, y, width=0.2, color=colors_q2[model], label=label_q2[model])
        axes[r_][0].set_xlim(0, 10)
        axes[r_][1].set_xlim(0, 10)
        axes[r_][2].set_yscale('log')
        axes[r_][3].set_yscale('log')
        titles = ['Time from last infection or reset\nto first progression (CDF)', 'Latent spells ending in progression:\nduration (CDF)',
                  'Latent clock resets before\nfirst progression (share of people)', 'Episodes of TB per person\n(share of people)']
        for col, title in enumerate(titles):
            axes[r_][col].set_title(f'{title}, λ = {lam:.2f}', fontsize=8)
        axes[r_][0].set_xlabel('Years')
        axes[r_][1].set_xlabel('Years')
        axes[r_][2].set_xlabel('Resets')
        axes[r_][3].set_xlabel('Episodes')
    axes[0][0].legend(fontsize=7, frameon=False)
    fig.suptitle('Exp 17, question 2: per person, latent reinfection on (12 sets, weighted; tb_LAI_TPT resets = every latent reinfection, two-state = remote → recent)', fontsize=10)
    fig.tight_layout()
    fig.savefig(outdir / 'q2_per_person.png', dpi=160)
    print('Saved q2_reinfection_adds.png, q2_curves.png, q2_per_person.png')
    return


fn_latent = outdir / 'q2_latent_spells.obj'


def spell_curves(h_prog, h_clear, h_open):
    """
    Latent spell durations (counts by duration in steps, ending in progression, clearance, or still open at the end) as Kaplan-Meier
    curves: S (share of spells still latent at each duration) and the cumulative shares ended by progression and by clearance
    """
    n = max(len(h_prog), len(h_clear), len(h_open))
    pad = lambda h: np.pad(np.asarray(h, float), (0, n - len(h)))
    hp, hc, ho = pad(h_prog), pad(h_clear), pad(h_open)
    at_risk = np.cumsum((hp + hc + ho)[::-1])[::-1] # Spells lasting at least d steps
    # S[d] is the share still latent after d steps; spells can end at d = 0 (begun and ended within one step, possible in v4)
    S, cp, cc = np.ones(n + 1), np.zeros(n + 1), np.zeros(n + 1)
    s = 1.0
    for d in range(n):
        if at_risk[d] > 0:
            cp[d + 1] = cp[d] + s*hp[d]/at_risk[d]
            cc[d + 1] = cc[d] + s*hc[d]/at_risk[d]
            s *= 1 - (hp[d] + hc[d])/at_risk[d]
        else:
            cp[d + 1], cc[d + 1] = cp[d], cc[d]
        S[d + 1] = s
    # Aligned so that index d is after d steps: S[d] = share lasting more than d steps; by_prog[d], by_clear[d] = ended so by d steps
    return sc.objdict(S=S[1:], by_prog=cp[1:], by_clear=cc[1:])


def latent_spells():
    """ Weighted mean Kaplan-Meier curves of latent spell duration, reinfection on, for each model and lam (seeds pooled within a set) """
    lai = sc.load(fn_q2_lai)
    v4 = sc.load(rv.fn_q2)
    w = np.asarray(v4.weights)
    out = sc.objdict()
    for lam in rv.lams:
        for model in models_q2:
            per_set = []
            for i in range(rv.n_sets):
                runs = [lai.runs[f'p{i}_lam{lam}_on_s{s}'] if model == 'lai_tpt' else v4.runs[f'p{i}_{model}_lam{lam}_on_s{s}'] for s in rv.seeds_q2]
                n = max(len(r.hists.spell_prog) for r in runs) + 1
                tot = lambda name: sum(np.pad(r.hists[name], (0, n - len(r.hists[name]))) for r in runs)
                per_set.append(spell_curves(tot('spell_prog'), tot('spell_clear'), tot('spell_open')))
            m = min(len(c.S) for c in per_set)
            out[f'{model}_lam{lam}'] = sc.objdict({k: w @ np.array([c[k][:m] for c in per_set]) for k in ['S', 'by_prog', 'by_clear']})
            # Mean duration, up to 30 years (the area under S), per set and weighted
            out[f'{model}_lam{lam}'].mean_30y = float(w @ np.array([c.S[:m].sum()/12 for c in per_set]))
        # At lam = 0 the only spell is the first, so the chains' latent share is its exact survival curve
    for model in ['exact'] + [f'phi{phi}' for phi in rv.phis]:
        S = [np.interp(np.arange(30*12 + 1)/12, v4.chains[f'p{i}_{model}_lam0.0_on'].t, v4.chains[f'p{i}_{model}_lam0.0_on'].latent) for i in range(rv.n_sets)]
        out[f'chain_{model}_lam0.0'] = w @ np.array(S)
    out.dt_year = 1/12
    sc.save(fn_latent, out)
    return out


def plot_latent_spells(out=None):
    """ Latent spell durations: still latent, ended by progression, ended by clearance; one row per lam """
    out = sc.load(fn_latent) if out is None else out
    fig, axes = plt.subplots(len(rv.lams), 3, figsize=(13, 3.2*len(rv.lams)))
    for r_, lam in enumerate(rv.lams):
        for model in models_q2:
            c = out[f'{model}_lam{lam}']
            t = np.arange(len(c.S))*out.dt_year
            kw = dict(color=colors_q2[model], lw=1.8 if model == 'lai_tpt' else 1.2, ls='--' if model == 'lai_tpt' else '-',
                      label=f'{label_q2[model]} (mean {c.mean_30y:.2f} y)')
            axes[r_][0].plot(t, c.S, **kw)
            axes[r_][1].plot(t, c.by_prog, **kw)
            axes[r_][2].plot(t, c.by_clear, **kw)
        if lam == 0:
            tc = np.arange(30*12 + 1)/12
            axes[r_][0].plot(tc, out['chain_exact_lam0.0'], color='0.15', lw=1, ls=':', label="tb_LAI_TPT's model, solved exactly")
            for phi in rv.phis:
                axes[r_][0].plot(tc, out[f'chain_phi{phi}_lam0.0'], color=colors_q2[f'phi{phi}'], lw=0.8, ls=':')
        axes[r_][0].set_yscale('log')
        axes[r_][0].set_ylim(1e-4, 1.1)
        for col, title in enumerate(['Still latent (log scale)', 'Ended by progression (cumulative)', 'Ended by clearance (cumulative)']):
            axes[r_][col].set_xlim(0, 10)
            axes[r_][col].set_title(f'{title}, λ = {lam:.2f}/y', fontsize=9)
            axes[r_][col].set_xlabel('Years since the latent spell began')
        axes[r_][0].set_ylabel('Share of latent spells')
        axes[r_][0].legend(fontsize=7, frameon=False, loc='lower left')
    fig.suptitle('Exp 17, question 2: how long latent spells last, and how they end (reinfection on; Kaplan–Meier; 12 sets, weighted; '
                 'dotted at λ = 0: exact chains)', fontsize=10)
    fig.tight_layout()
    fig.savefig(outdir / 'q2_latent_spells.png', dpi=160)
    print('Saved q2_latent_spells.png')
    for lam in rv.lams:
        print(f'lam {lam:.2f}: mean spell (to 30 y) ' + ', '.join(f'{m} {out[f"{m}_lam{lam}"].mean_30y:.3f}' for m in models_q2)
              + '; still latent at 2 / 5 / 10 y: ' + '; '.join(f'{m} ' + '/'.join(f'{out[f"{m}_lam{lam}"].S[12*y]:.4f}' for y in [2, 5, 10]) for m in models_q2)
              + '; ended by progression: ' + ', '.join(f'{m} {out[f"{m}_lam{lam}"].by_prog[-1]:.4f}' for m in models_q2))
    return


fn_pp = outdir / 'q2_time_latent.obj'


def time_latent():
    """ Each person's total time latent over 30 years, and their number of latent spells: weighted mean distributions (seed 0 of each set) """
    lai = sc.load(outdir / 'lai_tpt_perperson.obj')
    v4 = sc.load(outdir / 'q2_v4_perperson.obj')
    w = sc.load(rv.fn_q2).weights
    out = sc.objdict()
    for lam in rv.lams:
        for model in models_q2:
            get = lambda i: lai[f'p{i}_lam{lam}_on_s0'] if model == 'lai_tpt' else v4[f'p{i}_{model}_lam{lam}_on_s0']
            m = min(len(get(i).time_latent) for i in range(rv.n_sets))
            d = sc.objdict()
            for name in ['time_latent', 'time_latent_open']:
                d[name] = w @ np.array([get(i)[name][:m]/get(i).n for i in range(rv.n_sets)])
            d.n_spells = w @ np.array([get(i).n_spells/get(i).n for i in range(rv.n_sets)])
            t = np.arange(m)*get(0).dt_year
            cdf = np.cumsum(d.time_latent)
            d.t = t
            d.quantiles = {q: float(t[np.searchsorted(cdf, q)]) for q in [0.25, 0.5, 0.75, 0.9, 0.99, 0.999]}
            d.mean = float(t @ d.time_latent)
            d.over = {y: float(d.time_latent[t > y].sum()) for y in [1, 2, 5, 10, 20]}
            out[f'{model}_lam{lam}'] = d
    sc.save(fn_pp, out)
    return out


def plot_time_latent(out=None):
    """
    Density plots (per year, linear axes), one column per lam: each person's total time latent over 30 years (top), and the duration of
    individual latent spells (bottom, from the Kaplan-Meier curves of latent_spells(), so spells still latent at 30 years are allowed for)
    """
    out = sc.load(fn_pp) if out is None else out
    spells = sc.load(fn_latent)
    xmax = 20 # The same x-axis in every panel
    fig, axes = plt.subplots(2, len(rv.lams), figsize=(4.6*len(rv.lams), 8.0))
    for j, lam in enumerate(rv.lams):
        for model in models_q2:
            kw = dict(color=colors_q2[model], lw=1.8 if model == 'lai_tpt' else 1.2, ls='--' if model == 'lai_tpt' else '-')
            d = out[f'{model}_lam{lam}']
            q = d.quantiles
            # Durations are whole months, each plotted at its own value; 0 (v4 only: begun and ended within a step) is merged into 1 month
            first = lambda x: np.concatenate([[x[0] + x[1]], x[2:]])
            axes[0][j].plot(d.t[1:], first(d.time_latent)*12,
                            label=f'{label_q2[model]}\n  median {q[0.5]:.2f} y, mean {d.mean:.2f} y, 90% {q[0.9]:.1f} y', **kw)
            S = spells[f'{model}_lam{lam}'].S
            ts = np.arange(len(S))*spells.dt_year
            axes[1][j].plot(ts[1:], first(-np.diff(np.concatenate([[1.0], S])))*12, label=f'{label_q2[model]} (mean {spells[f"{model}_lam{lam}"].mean_30y:.2f} y)', **kw)
        for i, (title, xlabel) in enumerate([('Total time latent per person, over 30 years', 'Total years latent'),
                                             ('Duration of each latent spell', 'Years from the start of the spell')]):
            ax = axes[i][j]
            ax.set_xlim(0, xmax)
            ax.set_ylim(bottom=0)
            ax.set_title(f'{title}, λ = {lam:.2f}/y', fontsize=9)
            ax.set_xlabel(xlabel)
            ax.legend(fontsize=6.5, frameon=False, loc='upper right')
        axes[0][j].set_ylabel('Density (share of people per year)')
        axes[1][j].set_ylabel('Density (share of spells per year)')
    fig.suptitle('Exp 17, question 2: time spent latent (reinfection on; 12 sets, weighted; top: 250k agents per set, seed 0; bottom: 4 seeds, Kaplan–Meier)', fontsize=10)
    note = ('A latent spell is one continuous period of latent infection: it starts with a new infection (at the start of the cohort, or after the person has cleared '
            'infection or recovered from non-infectious TB) and ends when they progress to TB or clear. Reinfection while latent does not end it: in the two-state model, '
            'recent → remote → reinfected → recent is one spell; in tb_LAI_TPT the clock resets but the spell continues. '
            'Top: each person\'s total over all their spells in 30 years (someone who clears and is infected again adds the second spell to the first). '
            'Bottom: each spell on its own, still-latent spells at 30 years allowed for (Kaplan–Meier). At λ = 0 no one is infected twice, so the two rows are the same. '
            'Durations are counted in whole months, each plotted at its own value. v4 dates a reinfection at the end of its month, so some of its spells last 0 months '
            '(1.8% at λ = 0.10, 3.3% at 0.60); these are plotted with the 1-month spells. v0.10 infects at the start of the month, so spells that start with '
            'a reinfection are recorded about half a month longer in v0.10 than in v4.')
    import textwrap
    fig.text(0.01, 0.005, textwrap.fill(note, 235), fontsize=7.5, ha='left', va='bottom')
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(outdir / 'q2_time_latent.png', dpi=160)
    print('Saved q2_time_latent.png')
    for lam in rv.lams:
        for model in models_q2:
            d = out[f'{model}_lam{lam}']
            print(f'lam {lam:.2f} {model:9s} mean {d.mean:.3f} y; quantiles ' + ' '.join(f'{k}: {v:.2f}' for k, v in d.quantiles.items())
                  + '; share over ' + ' '.join(f'{k} y: {v:.4f}' for k, v in d.over.items())
                  + f'; still latent at 30 y {d.time_latent_open.sum():.4f}; 2+ spells {d.n_spells[2:].sum():.4f}')
    return


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'q1'
    if which == 'perperson':
        plot_time_latent(sc.load(fn_pp) if (len(sys.argv) > 2 and sys.argv[2] == 'plot') else time_latent())
        sys.exit()
    if which == 'latent':
        plot_latent_spells(sc.load(fn_latent) if (len(sys.argv) > 2 and sys.argv[2] == 'plot') else latent_spells())
        sys.exit()
    if which == 'q1':
        plot_q1()
    elif which == 'q1n':
        plot_q1n()
    elif which == 'q2':
        out = sc.load(fn_q2_sum) if (len(sys.argv) > 2 and sys.argv[2] == 'plot') else summarise_q2()
        print_q2(out)
        plot_q2(out)
    else:
        raise ValueError(f'Unknown part "{which}"')
