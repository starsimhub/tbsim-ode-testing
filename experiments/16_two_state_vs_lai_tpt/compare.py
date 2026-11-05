"""
Exp 16: summaries and figures from the saved outputs of run_lai_tpt.py (tb_LAI_TPT as run, v0.10) and run_fit.py (two-state fits and the
exact reference, v4). See README.md.

v4 environment. Run from the repo root: v4/.venv/bin/python experiments/16_two_state_vs_lai_tpt/compare.py; add `plot` to redraw the
figures from summary.obj only
"""
import sys
import json
from pathlib import Path
import numpy as np
import pandas as pd
import sciris as sc
import matplotlib.pyplot as plt

outdir = Path(__file__).resolve().parent
ink, muted, surface = '#0b0b0b', '#52514e', '#fcfcfb'
curve_col = dict(A='#2a78d6', B='#eb6834', C='#1baf7a') # Categorical slots 1-3: the three progression curves
ramp = ['#86b6ef', '#5598e7', '#2a78d6', '#1c5cab', '#104281', '#0d366b'] # Sequential blue, light -> dark with phi
benchmarks = { # The data behind the targets (tb_LAI_TPT cohort_progression.py): placebo-arm cumulative incidence by year since infection
    'Ferebee 1970': dict(N0=7744, years=list(range(1, 11)), cum=[0.04584, 0.04817, 0.05165, 0.05424, 0.05669, 0.05824, 0.05992, 0.06095, 0.06147, 0.06185]),
    'Sutherland 1968': dict(N0=2672, years=list(range(1, 11)), cum=[0.05842, 0.08003, 0.08813, 0.09173, 0.09263, 0.09353, 0.09533, 0.09623, 0.09623, 0.09623]),
}
hazard_years = [0.5, 1, 2, 5, 10]


def wquantile(x, w, q):
    """ Weighted quantile (midpoint rule on the sorted cumulative weight) """
    x, w = np.asarray(x, float), np.asarray(w, float)
    o = np.argsort(x)
    cw = (np.cumsum(w[o]) - 0.5*w[o])/w.sum()
    return float(np.interp(q, cw, x[o]))


def wsummary(x, w):
    """ Weighted mean and weighted 25-75% range across parameter sets, as tb_LAI_TPT summarises them """
    return dict(mean=float(np.average(x, weights=w)), lo=wquantile(x, w, 0.25), hi=wquantile(x, w, 0.75))


def abm_curves(runs, i):
    """ tb_LAI_TPT as run, parameter set i: seed-mean curves on the grid 0, 1/12, ... (adding t = 0), and seed SEs of the scalars """
    r = [x for x in runs if x.param == i]
    t = np.concatenate([[0], r[0].t])
    def mean(key, start):
        return np.concatenate([[start], np.mean([x[key] for x in r], axis=0)])
    out = sc.objdict(t=t, A=mean('A', 0), B=mean('B', 0), C=mean('C', 0), latent=mean('latent', 1), exit_cle=mean('exit_cle', 0),
                     exit_non=mean('exit_non', 0), exit_asy=mean('exit_asy', 0))
    prog = np.sum([x.prog_from_latent for x in r], axis=0)
    at_risk = np.sum([x.latent_start for x in r], axis=0)
    dt = 1/12
    out.hazard_counts = (prog, at_risk) # Per month, pooled over seeds: agents leaving INFECTION to NON/ASY/SYM, and agents in INFECTION at the start
    out.share_from_latent = float(np.mean([x.share_from_latent for x in r]))
    out.share_from_latent_se = float(np.std([x.share_from_latent for x in r], ddof=1)/np.sqrt(len(r)))
    out.Y1_seeds = np.array([np.interp(1, x.t, x.B) for x in r])
    out.Ratio_seeds = np.array([np.interp(1, x.t, x.B)/np.interp(10, x.t, x.B) for x in r])
    out.dt = dt
    return out


def hazard_at(c, years, window=0.25):
    """ The progression hazard among people still latent at each time (per year); for the ABM, pooled over a window of +/- 3 months """
    if 'hazard_counts' in c:
        prog, at_risk = c.hazard_counts
        tm = c.t[1:] - c.dt/2 # Each step's midpoint
        out = []
        for yr in years:
            m = np.abs(tm - yr) <= window
            out.append(prog[m].sum()/at_risk[m].sum()/c.dt if at_risk[m].sum() else np.nan)
        return np.array(out)
    return np.interp(years, c.t, c.hazard)


def scalars(c):
    """ The scalar outputs of one cohort """
    at = lambda key, yr: float(np.interp(yr, c.t, c[key]))
    m50 = c.t <= 50 + 1e-9
    out = dict(Y1=at('B', 1), Ratio=at('B', 1)/at('B', 10), A50=at('A', 50), B50=at('B', 50), C50=at('C', 50),
               exit_cle=at('exit_cle', 50), exit_non=at('exit_non', 50), exit_asy=at('exit_asy', 50), latent_10y=at('latent', 10), latent_25y=at('latent', 25),
               time_latent=float(np.trapezoid(c.latent[m50], c.t[m50])), share_from_latent=float(c.share_from_latent))
    for yr, h in zip(hazard_years, hazard_at(c, hazard_years)):
        out[f'hazard_{yr:g}y'] = float(h)
    if 'recent' in c:
        for yr in [1, 5, 10]:
            out[f'recent_share_{yr}y'] = at('recent', yr)/at('latent', yr)
    return out


def build():
    """ Combine the saved outputs into per-set tables and summaries """
    lai = sc.load(outdir / 'lai_tpt.obj')
    cv = sc.load(outdir / 'curves.obj')
    fits = pd.read_csv(outdir / 'fits.csv')
    n_sets, phis, arms = cv['n_sets'], cv['phis'], cv['arms']
    w = np.asarray(lai.weights)
    models = {'as_run': {}, 'exact': {}}
    for i in range(n_sets):
        models['as_run'][i] = abm_curves(lai.runs, i)
        models['exact'][i] = cv['curves'][(i, 'exact')]
        for arm in arms:
            for phi in phis:
                models.setdefault((arm, phi), {})[i] = cv['curves'][(i, arm, phi)]

    rows = []
    for key, by_set in models.items():
        for i, c in by_set.items():
            s = scalars(c)
            ref = models['as_run'][i]
            win = (c.t >= 0.25) & (c.t <= 20)
            Bref = np.interp(c.t[win], ref.t, ref.B)
            s['maxgap_B_vs_as_run'] = float(np.max(np.abs(c.B[win]/Bref - 1)))
            label = key if isinstance(key, str) else f'arm{key[0]}_phi{key[1]:g}'
            rows.append(dict(model=label, arm=None if isinstance(key, str) else key[0], phi=None if isinstance(key, str) else key[1], param=i, weight=w[i], **s))
    by_set = pd.DataFrame(rows)
    ref = by_set[by_set.model == 'as_run'].set_index('param')
    cols = [c for c in by_set.columns if c not in ['model', 'arm', 'phi', 'param', 'weight']]
    summ = {}
    for model, d in by_set.groupby('model', sort=False):
        d = d.set_index('param').loc[range(n_sets)]
        summ[model] = {c: wsummary(d[c], w) for c in cols if d[c].notna().all()}
        summ[model]['ratio_to_as_run'] = {c: wsummary(d[c]/ref[c], w) for c in ['B50', 'A50', 'C50', 'time_latent', 'latent_10y', 'exit_cle', 'exit_non', 'exit_asy']}
    # tb_LAI_TPT as run: the seed SE of each set's progression observables
    summ['as_run']['seed_se'] = {f'param{i:02d}': dict(Y1=float(models['as_run'][i].Y1_seeds.std(ddof=1)/2), Ratio=float(models['as_run'][i].Ratio_seeds.std(ddof=1)/2),
                                                       share_from_latent=models['as_run'][i].share_from_latent_se) for i in range(n_sets)}

    # Question 1: posterior weight meeting the targets, by arm and phi; tb_LAI_TPT as run and the exact model as references
    q1 = {}
    for (arm, phi), d in fits.groupby(['arm', 'phi']):
        wt = d.sample_weight
        q1[f'arm{arm}_phi{phi:g}'] = dict(within1=float(wt[d.zmax <= 1].sum()), within2=float(wt[d.zmax <= 2].sum()),
                                          zmax_median=wquantile(d.zmax, wt, 0.5), zmax_q90=wquantile(d.zmax, wt, 0.9),
                                          z_Y1_median=wquantile(d.z_Y1, wt, 0.5), z_Ratio_median=wquantile(d.z_Ratio, wt, 0.5))
        if arm == 2:
            bounds = dict(inf_cle=(0.5, 1.5), non_rec=(0.14, 0.23), non_asy=(0.21, 0.28))
            for k, (lo, hi) in bounds.items():
                at_bound = np.isclose(d[k], lo) | np.isclose(d[k], hi)
                q1[f'arm{arm}_phi{phi:g}'][f'{k}_ratio_median'] = wquantile(d[k]/d[f'{k}_sample'], wt, 0.5)
                q1[f'arm{arm}_phi{phi:g}'][f'{k}_at_bound'] = float(wt[at_bound].sum())
    d0 = fits[(fits.arm == 1) & (fits.phi == 0)]
    for name, (a, b) in dict(as_run=('z_Y1_calib', 'z_Ratio_calib'), exact=('z_Y1_exact', 'z_Ratio_exact')).items():
        zm = np.maximum(d0[a].abs(), d0[b].abs())
        q1[name] = dict(within1=float(d0.sample_weight[zm <= 1].sum()), within2=float(d0.sample_weight[zm <= 2].sum()), zmax_median=wquantile(zm, d0.sample_weight, 0.5),
                        z_Y1_median=wquantile(d0[a], d0.sample_weight, 0.5), z_Ratio_median=wquantile(d0[b], d0.sample_weight, 0.5))
    for arm in arms:
        ok = fits[(fits.arm == arm) & (fits.zmax <= 1)].groupby('sample').phi.max()
        wt = fits.groupby('sample').sample_weight.first()
        q1[f'arm{arm}_largest_phi_within1'] = {f'{phi:g}': float(wt[ok.reindex(wt.index) == phi].sum()) for phi in phis}
    return dict(models=models, by_set=by_set, summary=summ, q1=q1, phis=phis, arms=arms, n_sets=n_sets, weights=w, v4_ode_check=cv['v4_ode_check'],
                lai_check=lai.check)


def style(ax):
    ax.set_facecolor(surface)
    for side in ['top', 'right']:
        ax.spines[side].set_visible(False)
    ax.grid(alpha=0.25, lw=0.6)


def wmean_curve(by_set, key, w, t):
    """ Weighted mean across parameter sets of one curve, on grid t """
    return np.average([np.interp(t, c.t, c[key]) for c in by_set.values()], axis=0, weights=w)


def plot(S):
    models, phis, arms, w, q1 = S['models'], S['phis'], S['arms'], S['weights'], S['q1']
    xs = np.arange(len(phis))
    labels = [f'{phi:g}' for phi in phis]

    # 1. Targets met
    fig, ax = plt.subplots(figsize=(7, 4.6), facecolor=surface)
    style(ax)
    for arm, ls, mk in [(1, '-', 'o'), (2, '--', 's')]:
        for lev, col in [('within1', ramp[3]), ('within2', ramp[0])]:
            ax.plot(xs, [100*q1[f'arm{arm}_phi{l}'][lev] for l in labels], color=col, ls=ls, marker=mk, ms=6, lw=1.8,
                    label=f'Arm {arm} ({"cascade fixed" if arm == 1 else "cascade freed"}), both |z| ≤ {lev[-1]}')
    for name, txt, col in [('as_run', 'tb_LAI_TPT as run', muted), ('exact', 'Exact r₀e^(−kτ)', '#c3c2b7')]:
        for lev, ls in [('within1', ':'), ('within2', '-.')]:
            ax.axhline(100*q1[name][lev], color=col, lw=1, ls=ls)
        ax.annotate(f'{txt}: {100*q1[name]["within1"]:.0f}% / {100*q1[name]["within2"]:.0f}% (1σ / 2σ)', (0, 100*q1[name]['within1']), textcoords='offset points',
                    xytext=(2, -11 if name == 'as_run' else 4), fontsize=7, color=col)
    ax.set_xticks(xs, labels)
    ax.set_ylim(0, 105)
    ax.set_xlabel('φ: remote progression ÷ recent progression')
    ax.set_ylabel('Posterior weight meeting both progression targets (%)')
    ax.set_title('Can the two-state model meet the targets? (345 samples)', fontsize=10, color=ink)
    ax.legend(frameon=False, fontsize=7, loc='lower left')
    fig.tight_layout()
    fig.savefig(outdir / 'targets_met.png', dpi=150, facecolor=surface)
    plt.close(fig)

    # 2. Progression by time since infection, one figure per arm, one panel per phi
    t = np.linspace(0, 20, 241)
    ymax = 100*max(wmean_curve(m, 'A', w, t).max() for m in models.values())*1.05 # One y range for every panel
    for arm in arms:
        fig, axes = plt.subplots(2, 3, figsize=(16, 9), facecolor=surface, sharey=True)
        for ax, phi in zip(axes.flat, phis):
            style(ax)
            for cv in ['A', 'B', 'C']:
                ax.plot(t, 100*wmean_curve(models['as_run'], cv, w, t), color=curve_col[cv], lw=2, ls='-')
                ax.plot(t, 100*wmean_curve(models[(arm, phi)], cv, w, t), color=curve_col[cv], lw=1.6, ls='--')
                ax.plot(t, 100*wmean_curve(models['exact'], cv, w, t), color=curve_col[cv], lw=1, ls=':')
            for (name, b), mk in zip(benchmarks.items(), ['o', 's']):
                cum = np.array(b['cum'])
                se = np.sqrt(cum*(1 - cum)/b['N0'])
                ax.errorbar(b['years'], 100*cum, yerr=196*se, fmt=mk, ms=4, color=muted, mfc=surface, capsize=2, lw=0.8, label=name if phi == phis[0] else None)
            ax.set_title(f'φ = {phi:g}', fontsize=10, color=ink)
            ax.set_xlim(0, 20)
            ax.set_ylim(0, ymax)
            ax.set_xlabel('Years since infection')
        for ax in axes[:, 0]:
            ax.set_ylabel('Cumulative incidence (%)')
        handles = [plt.Line2D([], [], color=curve_col[c], lw=2, label=l) for c, l in [('A', 'A: any progression (NON or ASY)'), ('B', 'B: at least asymptomatic'), ('C', 'C: symptomatic')]]
        handles += [plt.Line2D([], [], color=muted, lw=2, ls=s, label=l) for s, l in [('-', 'tb_LAI_TPT as run (v0.10, monthly)'), ('--', f'Two-state, arm {arm}'), (':', 'Exact r₀e^(−kτ) (reference)')]]
        handles += axes.flat[0].get_legend_handles_labels()[0]
        fig.legend(handles=handles, loc='lower center', ncol=4, frameon=False, fontsize=8)
        fig.suptitle(f'Progression by time since infection, no-reinfection cohort: weighted mean of the 12 face-validity sets (arm {arm}: '
                     f'{"cascade fixed" if arm == 1 else "inf_cle, non_rec, non_asy freed within priors"})', fontsize=11, color=ink)
        fig.tight_layout(rect=(0, 0.07, 1, 0.97))
        fig.savefig(outdir / f'progression_by_time_arm{arm}.png', dpi=150, facecolor=surface)
        plt.close(fig)

    # 3. Latent-infection outcomes (tb_state_outcomes bottom-right panel, for a cohort): first exit by 50 years, and still latent at 10 years
    summ = S['summary']
    cats = [('exit_cle', 'Cleared'), ('exit_non', 'To NON'), ('exit_asy', 'To ASY'), ('latent_10y', 'Still latent at 10 y')]
    fig, axes = plt.subplots(1, 2, figsize=(15, 4.8), facecolor=surface, sharey=True)
    for ax, arm in zip(axes, arms):
        style(ax)
        keys = [('as_run', 'tb_LAI_TPT as run', '#52514e'), ('exact', 'Exact (reference)', '#c3c2b7')] + [(f'arm{arm}_phi{phi:g}', f'Two-state, φ = {phi:g}', ramp[j]) for j, phi in enumerate(phis)]
        bw = 0.8/len(keys)
        for j, (k, lab, col) in enumerate(keys):
            m = np.array([summ[k][c]['mean'] for c, _ in cats])
            lo = np.array([summ[k][c]['lo'] for c, _ in cats])
            hi = np.array([summ[k][c]['hi'] for c, _ in cats])
            ax.bar(np.arange(len(cats)) - 0.4 + bw*(j + 0.5), 100*m, bw*0.92, color=col, label=lab,
                   yerr=100*np.vstack([np.clip(m - lo, 0, None), np.clip(hi - m, 0, None)]), capsize=2, error_kw=dict(lw=0.8))
        ax.set_xticks(np.arange(len(cats)), [l for _, l in cats])
        ax.set_yscale('log')
        ax.set_ylabel('% of infections (weighted mean; bars: weighted 25–75% across sets)')
        ax.set_title(f'Latent-infection outcomes, by 50 years (arm {arm})', fontsize=10, color=ink)
        ax.legend(frameon=False, fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(outdir / 'latent_outcomes.png', dpi=150, facecolor=surface)
    plt.close(fig)

    # 4. Latent prevalence, recent share and latent hazard for param00 (set B), both arms
    fig, axes = plt.subplots(2, 3, figsize=(16, 8.5), facecolor=surface)
    t = np.linspace(0, 20, 241)
    for row, arm in enumerate(arms):
        ax = axes[row, 0]
        style(ax)
        ax.plot(models['as_run'][0].t, 100*models['as_run'][0].latent, color=ink, lw=2, label='tb_LAI_TPT as run')
        ax.plot(models['exact'][0].t, 100*models['exact'][0].latent, color='#c3c2b7', lw=1.5, ls=':', label='Exact (reference)')
        for j, phi in enumerate(phis):
            c = models[(arm, phi)][0]
            ax.plot(c.t, 100*c.latent, color=ramp[j], lw=1.3, ls='--', label=f'Two-state, φ = {phi:g}')
        ax.set_xlim(0, 20)
        ax.set_yscale('log')
        ax.set_ylim(1e-6, 150)
        ax.set_ylabel('% of the cohort still latent')
        ax.set_title(f'Set B: latent prevalence (arm {arm})', fontsize=10, color=ink)
        ax.legend(frameon=False, fontsize=7)
        ax = axes[row, 1]
        style(ax)
        for j, phi in enumerate(phis):
            c = models[(arm, phi)][0]
            ax.plot(c.t[1:], 100*c.recent[1:]/c.latent[1:], color=ramp[j], lw=1.5, label=f'φ = {phi:g}')
        ax.set_xlim(0, 5)
        ax.set_ylabel('% of latent infections that are recent')
        ax.set_title(f'Set B: recent share of latent infection (arm {arm})', fontsize=10, color=ink)
        ax.legend(frameon=False, fontsize=7)
        ax = axes[row, 2]
        style(ax)
        hy = np.array([0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 7, 10])
        h = hazard_at(models['as_run'][0], hy)
        ax.plot(hy, np.where(h > 0, h, np.nan), 'o', color=ink, ms=5, label='tb_LAI_TPT as run (ABM, ±3 months)')
        ax.plot(models['exact'][0].t[1:], models['exact'][0].hazard[1:], color='#c3c2b7', lw=1.5, ls=':', label='Exact (reference)')
        for j, phi in enumerate(phis):
            c = models[(arm, phi)][0]
            ax.plot(c.t[1:], c.hazard[1:], color=ramp[j], lw=1.3, ls='--', label=f'φ = {phi:g}')
        ax.set_yscale('log')
        ax.set_ylim(1e-5, 2)
        ax.set_xlim(0, 10)
        ax.set_ylabel('Progression hazard among people still latent (per year)')
        ax.set_title(f'Set B: latent progression hazard (arm {arm})', fontsize=10, color=ink)
        ax.legend(frameon=False, fontsize=7)
    for ax in axes[-1]:
        ax.set_xlabel('Years since infection')
    fig.tight_layout()
    fig.savefig(outdir / 'latent_dynamics.png', dpi=150, facecolor=surface)
    plt.close(fig)

    # 5. Scalar outputs against phi, relative to tb_LAI_TPT as run
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6), facecolor=surface)
    for ax, (key, lab) in zip(axes[:2], [('B50', 'Cumulative progression to ASY by 50 years, B(50)'), ('time_latent', 'Expected time spent latent (to 50 years)')]):
        style(ax)
        for arm, ls, mk in [(1, '-', 'o'), (2, '--', 's')]:
            m = np.array([summ[f'arm{arm}_phi{l}']['ratio_to_as_run'][key]['mean'] for l in labels])
            lo = np.array([summ[f'arm{arm}_phi{l}']['ratio_to_as_run'][key]['lo'] for l in labels])
            hi = np.array([summ[f'arm{arm}_phi{l}']['ratio_to_as_run'][key]['hi'] for l in labels])
            ax.errorbar(xs + (0.08 if arm == 2 else -0.08), m, yerr=[m - lo, hi - m], color=ramp[3] if arm == 1 else ramp[1], ls=ls, marker=mk, ms=6, capsize=3, label=f'Two-state, arm {arm}')
        ax.axhline(summ['exact']['ratio_to_as_run'][key]['mean'], color='#c3c2b7', lw=1.2, ls=':', label='Exact (reference)')
        ax.axhline(1, color=muted, lw=0.8)
        ax.set_xticks(xs, labels)
        ax.set_xlabel('φ')
        ax.set_ylabel('÷ tb_LAI_TPT as run (weighted mean; 25–75%)')
        ax.set_title(lab, fontsize=10, color=ink)
        ax.legend(frameon=False, fontsize=7)
    ax = axes[2]
    style(ax)
    for arm, ls, mk in [(1, '-', 'o'), (2, '--', 's')]:
        m = np.array([summ[f'arm{arm}_phi{l}']['share_from_latent']['mean'] for l in labels])
        ax.plot(xs, 100*m, color=ramp[3] if arm == 1 else ramp[1], ls=ls, marker=mk, ms=6, label=f'Two-state, arm {arm}')
    ax.axhline(100*summ['as_run']['share_from_latent']['mean'], color=ink, lw=1.2, label='tb_LAI_TPT as run')
    ax.axhline(100*summ['exact']['share_from_latent']['mean'], color='#c3c2b7', lw=1.2, ls=':', label='Exact (reference)')
    ax.set_xticks(xs, labels)
    ax.set_xlabel('φ')
    ax.set_ylabel('% (weighted mean across the 12 sets)')
    ax.set_title('Share of years 1–10 progression from people still latent at year 1', fontsize=10, color=ink)
    ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(outdir / 'scalars.png', dpi=150, facecolor=surface)
    plt.close(fig)
    return


if __name__ == '__main__':
    if 'plot' in sys.argv[1:]:
        S = sc.load(outdir / 'summary.obj')
    else:
        S = build()
        sc.save(outdir / 'summary.obj', S)
        S['by_set'].to_csv(outdir / 'dynamics_by_set.csv', index=False)
        (outdir / 'summary.json').write_text(json.dumps(dict(q1=S['q1'], summary=S['summary'], v4_ode_check=S['v4_ode_check'], lai_check=S['lai_check']), indent=2, default=float) + '\n')
    plot(S)

    summ, q1 = S['summary'], S['q1']
    print('Question 1: posterior weight with both |z| <= 1 / <= 2')
    for k in ['as_run', 'exact'] + [f'arm{a}_phi{p:g}' for a in S['arms'] for p in S['phis']]:
        v = q1[k]
        extra = '' if not k.startswith('arm2') else '; cascade ÷ sample (median): ' + ', '.join(f'{c} {v[f"{c}_ratio_median"]:.3f} (at bound {v[f"{c}_at_bound"]:.2f})' for c in ['inf_cle', 'non_rec', 'non_asy'])
        print(f'  {k:16s} {v["within1"]:.3f} / {v["within2"]:.3f}; median z: Y1 {v["z_Y1_median"]:+.2f}, Ratio {v["z_Ratio_median"]:+.2f}{extra}')
    print('Question 2: weighted mean [25%, 75%] across the 12 sets')
    keys = ['Y1', 'Ratio', 'A50', 'B50', 'C50', 'exit_cle', 'exit_non', 'exit_asy', 'latent_10y', 'latent_25y', 'time_latent', 'share_from_latent', 'hazard_1y', 'hazard_5y', 'maxgap_B_vs_as_run']
    for k in ['as_run', 'exact'] + [f'arm{a}_phi{p:g}' for a in S['arms'] for p in S['phis']]:
        print(f'  {k}:')
        print('    ' + '; '.join(f'{c} {summ[k][c]["mean"]:.4g} [{summ[k][c]["lo"]:.3g}, {summ[k][c]["hi"]:.3g}]' for c in keys if c in summ[k]))
        if k != 'as_run':
            print('    ÷ as run: ' + '; '.join(f'{c} {v["mean"]:.3f} [{v["lo"]:.3f}, {v["hi"]:.3f}]' for c, v in summ[k]['ratio_to_as_run'].items()))
