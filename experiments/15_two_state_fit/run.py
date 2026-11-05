"""
Exp 15: can v4's two-state latent model, with progression from the remote state (a fraction phi of recent progression), fit tb_LAI_TPT's
progression targets within its priors, and does it reproduce the exact model's progression by time since infection? See README.md.

Deterministic: one person's continuous-time Markov chain per structure (tbode4.latent), for each of tb_LAI_TPT's 345 posterior samples.

Run from the repo root: v4/.venv/bin/python experiments/15_two_state_fit/run.py; add `plot` to redraw the figure from outputs.obj
"""
import sys
import json
from pathlib import Path
import numpy as np
import pandas as pd
import sciris as sc
import matplotlib.pyplot as plt
from scipy.linalg import expm
from scipy.optimize import least_squares
import tbode4

outdir = Path(__file__).resolve().parent
samples_csv = outdir.parents[2] / 'tb_LAI_TPT/results/hm_real_20260728/trajectory_selection/collapsed_samples_eps1.2.csv'

phis = [0.0, 0.003, 0.01, 0.03, 0.1, 0.3] # Remote progression as a fraction of recent progression
bounds = dict(inf_non=(0.04, 0.45), inf_asy=(0.01, 0.60), k_latent=(0.0, 12.0), inf_non_remote=0.23, inf_asy_remote=0.20) # Priors (user, 2026-11-04)
feasible_tol = 1e-4   # Both log-ratios to the sample's targets within this
lams = [0.12, 0.60]   # Force of infection, per year, for the reinfection scenario (as in exp 14)
mu = 1/70             # Background deaths in the reinfection scenario
tmax, per_year = 50, 48
gap_from = 0.25
gap_years = [2, 5, 10, 50]
hazard_years = [0, 0.5, 1, 2, 5, 10]

# The data behind the targets: placebo-arm cumulative incidence by year since infection (tb_LAI_TPT cohort_progression.py)
benchmarks = {
    'Ferebee 1970': dict(N0=7744, years=list(range(1, 11)), cum=[0.04584, 0.04817, 0.05165, 0.05424, 0.05669, 0.05824, 0.05992, 0.06095, 0.06147, 0.06185]),
    'Sutherland 1968': dict(N0=2672, years=list(range(1, 11)), cum=[0.05842, 0.08003, 0.08813, 0.09173, 0.09263, 0.09353, 0.09533, 0.09623, 0.09623, 0.09623]),
}


def later_progression(p, structure):
    """ Progression between years 1 and 10, and the share of it from people still latent at year 1 (no reinfection, no deaths) """
    Q, x0, idx = tbode4.latent_chain(p, structure)
    x1, P9 = x0 @ expm(Q), expm(Q*9)
    from_latent, from_non = x1[idx.latent] @ P9[idx.latent, idx.ASY], x1[idx.NON]*P9[idx.NON, idx.ASY]
    return from_latent + from_non, from_latent/(from_latent + from_non)


def dynamics(p, structure):
    """ Cohort outputs for one parameter set: curve B, latent hazard, later progression, lifetime progression, time latent, with reinfection """
    t, X, idx = tbode4.cohort_curve(p, structure, tmax=tmax, per_year=per_year)
    later, share = later_progression(p, structure)
    out = sc.objdict(B=X[:, idx.ASY], hazard=tbode4.latent_hazard(X, idx), later=later, share_latent=share, lifetime=X[-1, idx.ASY],
                     time_latent=np.trapezoid(X[:, idx.latent].sum(axis=1), t))
    for lam in lams:
        out[f'reinf_{lam:g}'] = tbode4.ever_asy(p, structure, years=[30], lam=lam, mu=mu)[0]
    return out


def two_state_pars(base, f_non, phi, x):
    """ Two-state parameters from the fitted (r0, k_latent), the sample's other rates, and phi """
    r0, k_latent = x
    p = dict(base, inf_non=f_non*r0, inf_asy=(1 - f_non)*r0, k_latent=k_latent)
    p.update(inf_non_remote=min(phi*p['inf_non'], bounds['inf_non_remote']), inf_asy_remote=min(phi*p['inf_asy'], bounds['inf_asy_remote']))
    return p


def fit(base, f_non, phi, target):
    """ Fit (r0, k_latent) so that ProgressionY1 and ProgressionRatio equal the target, within the priors; best of several starts """
    lo = max(bounds['inf_non'][0]/f_non, bounds['inf_asy'][0]/(1 - f_non))
    hi = min(bounds['inf_non'][1]/f_non, bounds['inf_asy'][1]/(1 - f_non))
    if lo > hi:
        return None, np.inf, None

    def resid(x):
        got = tbode4.progression_targets(two_state_pars(base, f_non, phi, x), 'two-state')
        return np.array([np.log(got['ProgressionY1']/target['ProgressionY1']), np.log(got['ProgressionRatio']/target['ProgressionRatio'])])

    r0_base = base['inf_non'] + base['inf_asy']
    starts = [(np.clip(r0_base, lo, hi), np.clip(0.95*base['k'], 0.1, 11.9)), (np.clip(r0_base, lo, hi), 11.0), (0.5*(lo + hi), 3.0)]
    best = None
    for x0 in starts:
        res = least_squares(resid, x0, bounds=([lo, bounds['k_latent'][0]], [hi, bounds['k_latent'][1]]), xtol=1e-12, ftol=1e-12, gtol=1e-12)
        if best is None or res.cost < best.cost:
            best = res
    return best.x, float(np.max(np.abs(best.fun))), best.fun


def zscores(vals):
    """ Log-space z-scores of the progression observables against tb_LAI_TPT's targets """
    return {f'z_{k}': float(np.log(vals[k]/tbode4.target_means[k])/tbode4.target_sigmas[k]) for k in tbode4.target_means}


def run_sample(i, row):
    """ The exact model and the two-state fits at every phi, for one posterior sample """
    base = dict(inf_cle=row.inf_cle, inf_non=row.inf_non, inf_asy=row.inf_asy, non_rec=row.non_rec, non_asy=row.non_asy, pi=row.rr_reinfection_rec, k=row.k_asy)
    f_non = base['inf_non']/(base['inf_non'] + base['inf_asy'])
    target = tbode4.progression_targets(base, 'exact')
    A = dynamics(base, 'exact')
    t = np.linspace(0, tmax, tmax*per_year + 1)
    win = t >= gap_from - 1e-9
    gi = {yr: int(round(yr*per_year)) for yr in sorted(set(gap_years + hazard_years))}
    rows, curves = [], dict(A=A)
    for phi in phis:
        x, misfit, fun = fit(base, f_non, phi, target)
        out = dict(sample=i, sample_weight=row.sample_weight, phi=phi, k=base['k'], r0_exact=base['inf_non'] + base['inf_asy'], inf_cle=base['inf_cle'],
                   Y1_exact=target['ProgressionY1'], Ratio_exact=target['ProgressionRatio'], **{f'{k}_exact': v for k, v in zscores(target).items()},
                   misfit=misfit, feasible=misfit < feasible_tol)
        if x is None:
            rows.append(out)
            continue
        p = two_state_pars(base, f_non, phi, x)
        got = tbode4.progression_targets(p, 'two-state')
        out.update(r0=x[0], k_latent=x[1], inf_non=p['inf_non'], inf_asy=p['inf_asy'], inf_non_remote=p['inf_non_remote'], inf_asy_remote=p['inf_asy_remote'],
                   misfit_Y1_sigma=fun[0]/tbode4.target_sigmas['ProgressionY1'], misfit_Ratio_sigma=fun[1]/tbode4.target_sigmas['ProgressionRatio'],
                   Y1=got['ProgressionY1'], Ratio=got['ProgressionRatio'], **zscores(got),
                   zmax=max(abs(v) for v in zscores(got).values()), zmax_exact=max(abs(v) for v in zscores(target).values()),
                   at_bound=bool(np.isclose(x[1], bounds['k_latent'][1]) or np.isclose(x[0]*f_non, bounds['inf_non'][1]) or np.isclose(x[0]*(1 - f_non), bounds['inf_asy'][1])))
        D = dynamics(p, 'two-state')
        rel = D.B/np.where(A.B > 0, A.B, np.nan) - 1
        out['maxgap'] = float(np.nanmax(np.abs(rel[win])))
        for yr in gap_years:
            out[f'gap_{yr}y'] = float(rel[gi[yr]])
        for yr in hazard_years:
            out[f'hazard_{yr:g}y'] = float(D.hazard[gi[yr]])
            out[f'hazard_exact_{yr:g}y'] = float(A.hazard[gi[yr]])
            out[f'hazard_ratio_{yr:g}y'] = float(D.hazard[gi[yr]]/A.hazard[gi[yr]])
        out.update(later_ratio=D.later/A.later, share_latent=D.share_latent, share_latent_exact=A.share_latent,
                   lifetime_ratio=D.lifetime/A.lifetime, time_latent_ratio=D.time_latent/A.time_latent,
                   **{f'reinf_{lam:g}_ratio': D[f'reinf_{lam:g}']/A[f'reinf_{lam:g}'] for lam in lams})
        rows.append(out)
        curves[phi] = D
    return rows, curves


def wquantile(x, w, q):
    """ Weighted quantiles (midpoint rule on the sorted cumulative weight) """
    x, w = np.asarray(x, float), np.asarray(w, float)
    ok = ~np.isnan(x)
    x, w = x[ok], w[ok]
    if not len(x):
        return np.nan
    order = np.argsort(x)
    x, w = x[order], w[order]
    cw = (np.cumsum(w) - 0.5*w)/w.sum()
    return float(np.interp(q, cw, x))


def summarise(df):
    """ For each phi: the posterior weight with a feasible fit, and weighted median [10%, 90%] and max of each output over feasible samples """
    cols = [c for c in df.columns if c not in ['sample', 'sample_weight', 'phi', 'feasible', 'at_bound']]
    out = {}
    for phi, d in df.groupby('phi'):
        w = d.sample_weight.values
        ok = d.feasible.values
        s = dict(weight_feasible=float(w[ok].sum()/w.sum()), n_feasible=int(ok.sum()), weight_at_bound=float(w[d.at_bound.fillna(False).values.astype(bool)].sum()/w.sum()))
        s['misfit_sigma_all'] = {q: wquantile(np.maximum(np.abs(d.misfit_Y1_sigma), np.abs(d.misfit_Ratio_sigma)), w, qq) for q, qq in [('median', 0.5), ('q90', 0.9)]}
        s['misfit_sigma_all']['max'] = float(np.nanmax(np.maximum(np.abs(d.misfit_Y1_sigma), np.abs(d.misfit_Ratio_sigma))))
        f = d[ok]
        s['stats'] = {c: dict(median=wquantile(f[c], f.sample_weight, 0.5), q10=wquantile(f[c], f.sample_weight, 0.1), q90=wquantile(f[c], f.sample_weight, 0.9),
                              min=float(f[c].min()) if len(f) else np.nan, max=float(f[c].max()) if len(f) else np.nan) for c in cols}
        out[f'{phi:g}'] = s
    largest = df[df.feasible].groupby('sample').phi.max().reindex(df['sample'].unique(), fill_value=np.nan)
    wt = df.groupby('sample').sample_weight.first()
    out['largest_feasible_phi'] = {f'{phi:g}': float(wt[largest == phi].sum()) for phi in phis}
    out['largest_feasible_phi']['none'] = float(wt[largest.isna()].sum())
    return out


def plot(df, curves, summ):
    """ Feasibility by phi; set B's curve B(t) and latent hazard at each phi; side effects by phi """
    ink, muted, surface = '#0b0b0b', '#52514e', '#fcfcfb'
    ramp = ['#86b6ef', '#5598e7', '#2a78d6', '#1c5cab', '#104281', '#0d366b'] # One hue, light -> dark with phi
    t = np.linspace(0, tmax, tmax*per_year + 1)
    fig, axes = plt.subplots(1, 4, figsize=(20, 4.8), facecolor=surface)
    for ax in axes:
        ax.set_facecolor(surface)
        for side in ['top', 'right']:
            ax.spines[side].set_visible(False)
        ax.grid(alpha=0.25, lw=0.6)
    labels = [f'{phi:g}' for phi in phis]
    xs = np.arange(len(phis))

    ax = axes[0]
    by = {phi: d for phi, d in df.groupby('phi')}
    within = [100*by[phi].sample_weight[by[phi].zmax < 2].sum() for phi in phis]
    exact_within = 100*by[0].sample_weight[by[0].zmax_exact < 2].sum()
    ax.plot(xs, within, color=ramp[3], lw=2, marker='o', ms=8, label='Best fit within 2σ of both target means')
    ax.plot(xs, [100*summ[l]['weight_feasible'] for l in labels], color=ramp[0], lw=2, marker='s', ms=7, label="Reproduces the sample's own Y1 and ratio")
    ax.axhline(exact_within, color=muted, lw=1, ls='--')
    ax.annotate(f'Exact model (posterior): {exact_within:.0f}% within 2σ', (len(phis) - 1, exact_within), textcoords='offset points', xytext=(0, -12), ha='right', fontsize=8, color=muted)
    for x, v in zip(xs, within):
        ax.annotate(f'{v:.0f}%', (x, v), textcoords='offset points', xytext=(0, 8), ha='center', fontsize=8, color=muted)
    ax.set_xticks(xs, labels)
    ax.set_ylim(-5, 112)
    ax.set_xlabel('φ: remote progression ÷ recent progression')
    ax.set_ylabel('Posterior weight (%)')
    ax.set_title('(a) Fit to the progression targets, within the priors', fontsize=10, color=ink)
    ax.legend(frameon=False, fontsize=7, loc='lower left')

    c0 = curves[0]
    ax = axes[1]
    ax.plot(t, 100*c0['A'].B, color=ink, lw=2, label='Exact r·e^(−kτ)')
    for phi, col in zip(phis, ramp):
        if phi in c0:
            ok = bool(df[(df['sample'] == 0) & (df.phi == phi)].feasible.iloc[0])
            ax.plot(t, 100*c0[phi].B, color=col, lw=1.5, ls='-' if ok else '--', label=f'Two-state, φ = {phi:g}' + ('' if ok else ' (best, misses)'))
    for (name, b), mk in zip(benchmarks.items(), ['o', 's']):
        cum = np.array(b['cum'])
        se = np.sqrt(cum*(1 - cum)/b['N0'])
        ax.errorbar(b['years'], 100*cum, yerr=196*se, fmt=mk, ms=5, color=muted, mfc=surface, capsize=2, lw=1, label=name)
    ax.set_xlim(0, 12)
    ax.set_xlabel('Years since infection')
    ax.set_ylabel('Ever reached asymptomatic TB (%)')
    ax.set_title('(b) Set B (row 0): progression by time since infection', fontsize=10, color=ink)
    ax.legend(frameon=False, fontsize=7, loc='lower right')

    ax = axes[2]
    win = (t > 0) & (t <= 10)
    floor = 1e-6 # Below this the hazards are numerical noise, and progression is nil
    ax.plot(t[win], np.maximum(c0['A'].hazard[win], floor), color=ink, lw=2, label='Exact r₀·e^(−kτ)')
    for phi, col in zip(phis, ramp):
        if phi in c0:
            ax.plot(t[win], np.maximum(c0[phi].hazard[win], floor), color=col, lw=1.5, label=f'φ = {phi:g}')
    ax.set_yscale('log')
    ax.set_ylim(floor/2, 1.5)
    ax.set_xlabel('Years since infection')
    ax.set_ylabel('Progression hazard among people still latent (per year)')
    ax.set_title('(b) Set B: latent progression hazard', fontsize=10, color=ink)
    ax.legend(frameon=False, fontsize=7)

    ax = axes[3]
    measures = [('lifetime_ratio', 'Lifetime progression', '#2a78d6'), ('reinf_0.12_ratio', 'Ever ASY by 30 y, with reinfection (λ = 0.12)', '#eb6834'),
                ('time_latent_ratio', 'Time spent latent', '#1baf7a')]
    for j, (col, lab, color) in enumerate(measures):
        off = (j - 1)*0.18
        med = np.array([wquantile(by[phi][col], by[phi].sample_weight, 0.5) for phi in phis])
        lo = np.array([wquantile(by[phi][col], by[phi].sample_weight, 0.1) for phi in phis])
        hi = np.array([wquantile(by[phi][col], by[phi].sample_weight, 0.9) for phi in phis])
        ax.errorbar(xs + off, med, yerr=[med - lo, hi - med], fmt='o', ms=6, color=color, capsize=3, lw=1.5, label=lab)
    ax.axhline(1, color=muted, lw=0.8)
    ax.set_xticks(xs, labels)
    ax.set_xlabel('φ: remote progression ÷ recent progression')
    ax.set_ylabel('Two-state ÷ exact (weighted median, 10–90%)')
    ax.set_title('Side effects: best fits, all samples', fontsize=10, color=ink)
    ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(outdir / 'two_state_fit.png', dpi=150, facecolor=surface)
    return


if __name__ == '__main__' and 'plot' in sys.argv[1:]:
    saved = sc.load(outdir / 'outputs.obj')
    plot(saved['df'], saved['curves'], saved['summary'])
elif __name__ == '__main__':
    post = pd.read_csv(samples_csv)
    assert np.allclose(post.k_asy, post.k_non), 'k_asy and k_non differ in some samples'
    results = sc.parallelize(run_sample, iterkwargs=[dict(i=i, row=row) for i, row in post.iterrows()], ncpus=tbode4.ncpus)
    df = pd.DataFrame([r for rows, _ in results for r in rows])
    curves = {i: c for i, (_, c) in enumerate(results)}
    df.to_csv(outdir / 'fits.csv', index=False)
    summ = summarise(df)
    summ['settings'] = dict(samples_csv=str(samples_csv), phis=phis, bounds=bounds, feasible_tol=feasible_tol, lams=lams, mu=mu, tmax=tmax, per_year=per_year)
    (outdir / 'summary.json').write_text(json.dumps(summ, indent=2, default=float) + '\n')
    sc.save(outdir / 'outputs.obj', dict(df=df, curves=curves, summary=summ, t=np.linspace(0, tmax, tmax*per_year + 1)))
    plot(df, curves, summ)
    tbode4.save_versions(outdir, 'v4')

    fmt = lambda s, c, m=1: f"{m*s['stats'][c]['median']:.5g} [{m*s['stats'][c]['q10']:.5g}, {m*s['stats'][c]['q90']:.5g}], max {m*s['stats'][c]['max']:.5g}"
    print('Largest feasible phi (posterior weight):', {k: round(v, 3) for k, v in summ['largest_feasible_phi'].items()})
    for phi in phis:
        s = summ[f'{phi:g}']
        d = df[df.phi == phi]
        w = d.sample_weight.values
        print(f'phi = {phi:g}: weight feasible {s["weight_feasible"]:.3f} ({s["n_feasible"]} samples); at a prior bound {s["weight_at_bound"]:.3f}; '
              f'misfit to the sample (sigma) median {s["misfit_sigma_all"]["median"]:.3g}, q90 {s["misfit_sigma_all"]["q90"]:.3g}, max {s["misfit_sigma_all"]["max"]:.3g}')
        print(f'    all samples, max|z| to target means: two-state median {wquantile(d.zmax, w, 0.5):.3g} [q90 {wquantile(d.zmax, w, 0.9):.3g}], '
              f'exact {wquantile(d.zmax_exact, w, 0.5):.3g} [q90 {wquantile(d.zmax_exact, w, 0.9):.3g}]; weight within 2 sigma: two-state {w[d.zmax < 2].sum():.3f}, exact {w[d.zmax_exact < 2].sum():.3f}')
        print(f'    all samples, lifetime ratio median {wquantile(d.lifetime_ratio, w, 0.5):.5g} [q90 {wquantile(d.lifetime_ratio, w, 0.9):.5g}]; '
              f'with reinfection (0.12) median {wquantile(d["reinf_0.12_ratio"], w, 0.5):.5g} [q90 {wquantile(d["reinf_0.12_ratio"], w, 0.9):.5g}]')
        if s['n_feasible']:
            print('    feasible samples:')
            for c in ['k_latent', 'r0', 'maxgap', 'gap_50y', 'hazard_exact_1y', 'hazard_1y', 'hazard_exact_5y', 'hazard_5y', 'share_latent_exact', 'share_latent',
                      'lifetime_ratio', 'time_latent_ratio', 'reinf_0.12_ratio', 'reinf_0.6_ratio']:
                print(f'      {c:20s} {fmt(s, c)}')
