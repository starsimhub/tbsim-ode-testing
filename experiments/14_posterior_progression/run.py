"""
Exp 14: is v4's two-state latent model close to front-loaded progression (r·exp(−k·τ)) across tb_LAI_TPT's posterior? See README.md.

Deterministic: for each posterior sample, each structure is a continuous-time Markov chain for one person, solved exactly with
expm_multiply (generalises proposals/02_latent_progression/compare_structures.py to any parameter set).

Run from the repo root: v4/.venv/bin/python experiments/14_posterior_progression/run.py (about 3 minutes); add `plot` to redraw the figure
from outputs.obj
"""
import sys
import json
from pathlib import Path
import numpy as np
import pandas as pd
import sciris as sc
import matplotlib.pyplot as plt
from scipy.linalg import expm
from scipy.optimize import brentq
from scipy.sparse.linalg import expm_multiply
from scipy.stats import poisson
import tbode4

outdir = Path(__file__).resolve().parent
samples_csv = outdir.parents[2] / 'tb_LAI_TPT/results/hm_real_20260728/trajectory_selection/collapsed_samples_eps1.2.csv'

mu = 1/70                        # Background deaths in the reinfection scenario, per year
lams = [0.12, 0.60]              # Force of infection, per year: set B at 500/100k, and 5x (as in the options note)
reinf_years = [10, 30]
tmax, per_year = 50, 48          # Cohort grid: every 1/48 year to 50 years
gap_from = 0.25                  # The largest gap is taken over 3 months to 50 years
report_years = [0.25, 1, 5, 10, 50]
trunc = 1e-12                    # Truncated Poisson weight in the exact classes
thresholds = [0.01, 0.02, 0.05]


def n_classes(r0, k, tol):
    """ The number of Poisson-mixture classes (0..M-1) whose truncated weight is below tol """
    return max(2, int(poisson.isf(tol, r0/k)) + 2)


def chain(p, structure, k_latent=None, lam=0.0, latent_reinf=True, deaths=False):
    """
    One person's generator Q, initial vector x0, and the index of ASY (absorbing: 'ever reached asymptomatic TB')

    States: the structure's latent states, then CLE, NON, REC, ASY, DEAD. New infections (and the cohort's start) enter the latent
    states with weights `entry`. `structure` is 'exact' (Poisson-mixture classes) or 'two-state' (recent, remote)
    """
    r0 = p['inf_non'] + p['inf_asy']
    f_non = p['inf_non']/r0
    if structure == 'exact':
        nL = n_classes(r0, p['k'], trunc)
        entry = poisson.pmf(np.arange(nL), r0/p['k'])
        entry[-1] += 1 - entry.sum()
        prog = np.arange(nL)*p['k'] # Class m progresses at m·k
    else:
        nL = 2
        entry = np.array([1.0, 0.0])
        prog = np.array([r0, 0.0])
    CLE, NON, REC, ASY, DEAD = range(nL, nL + 5)
    n = nL + 5
    Q = np.zeros((n, n))

    def add(i, j, rate):
        Q[i, j] += rate
        Q[i, i] -= rate

    for i in range(nL):
        add(i, CLE, p['inf_cle'])
        add(i, NON, prog[i]*f_non)
        add(i, ASY, prog[i]*(1 - f_non))
    if structure != 'exact':
        add(0, 1, k_latent) # Recent -> remote
    add(NON, ASY, p['non_asy'])
    add(NON, REC, p['non_rec'])
    if lam > 0:
        for src in [CLE, REC]: # New infections of people no longer infected
            for j in range(nL):
                add(src, j, p['pi']*lam*entry[j])
        if latent_reinf:
            if structure == 'exact': # Clock reset: reassigned across the classes with the Poisson weights
                for i in range(nL):
                    for j in range(nL):
                        if i != j:
                            add(i, j, p['pi']*lam*entry[j])
            else: # Remote -> recent; reinfecting a recent infection does nothing
                add(1, 0, p['pi']*lam)
    if deaths:
        for i in range(n):
            if i not in (ASY, DEAD):
                add(i, DEAD, mu)
    x0 = np.zeros(n)
    x0[:nL] = entry
    return Q, x0, ASY


def cohort(p, structure, k_latent=None):
    """ The share ever reaching asymptomatic TB on the grid 0, 1/48, ..., 50 years (no reinfection, no deaths) """
    Q, x0, ASY = chain(p, structure, k_latent)
    return expm_multiply(Q.T, x0, start=0, stop=tmax, num=tmax*per_year + 1, endpoint=True)[:, ASY]


def ever_asy(p, structure, years, **kwargs):
    """ The share ever reaching asymptomatic TB at the given years """
    Q, x0, ASY = chain(p, structure, **kwargs)
    return np.array([(x0 @ expm(Q*t))[ASY] for t in years])


def y1(p, k_latent):
    """ ProgressionY1 of the two-state model """
    Q, x0, ASY = chain(p, 'two-state', k_latent)
    return (x0 @ expm(Q))[ASY]


def wquantile(x, w, q):
    """ Weighted quantiles (midpoint rule on the sorted cumulative weight) """
    order = np.argsort(x)
    x, w = np.asarray(x)[order], np.asarray(w)[order]
    cw = (np.cumsum(w) - 0.5*w)/w.sum()
    return np.interp(q, cw, x)


def run_sample(row):
    """ All outputs for one posterior sample """
    p = dict(inf_cle=row.inf_cle, inf_non=row.inf_non, inf_asy=row.inf_asy, non_rec=row.non_rec, non_asy=row.non_asy, k=row.k_asy, pi=row.rr_reinfection_rec)
    r0 = p['inf_non'] + p['inf_asy']
    t = np.linspace(0, tmax, tmax*per_year + 1)
    idx = {yr: int(round(yr*per_year)) for yr in report_years}
    curves = dict(A=cohort(p, 'exact'), B=cohort(p, 'two-state', p['k']))
    target = curves['A'][idx[1]]
    k_refit = brentq(lambda kl: y1(p, kl) - target, 1e-6, 1e4, xtol=1e-12, rtol=1e-12)
    curves['C'] = cohort(p, 'two-state', k_refit)

    out = dict(k=p['k'], r0=r0, r0_k=r0/p['k'], pi=p['pi'], k_latent_refit=k_refit, refit_ratio=k_refit/p['k'],
               classes_1e6=n_classes(r0, p['k'], 1e-6), classes_1e10=n_classes(r0, p['k'], 1e-10))
    window = t >= gap_from - 1e-9
    for s in 'ABC':
        out[f'Y1_{s}'] = curves[s][idx[1]]
        out[f'Ratio_{s}'] = curves[s][idx[1]]/curves[s][idx[10]] # tb_LAI_TPT's ProgressionRatio: B(1)/B(10) (cohort_progression.py)

    # Where the exact model's progression between years 1 and 10 comes from: people still latent at year 1, or already non-infectious
    Q, x0, ASY = chain(p, 'exact')
    nL, NON = ASY - 3, ASY - 2
    x1, P9 = x0 @ expm(Q), expm(Q*9)
    from_latent, from_non = x1[:nL] @ P9[:nL, ASY], x1[NON]*P9[NON, ASY]
    out.update(B1_to_10_A=from_latent + from_non, share_from_latent_A=from_latent/(from_latent + from_non), latent_at_1y_A=x1[:nL].sum(),
               non_at_1y_A=x1[NON], latent_hazard_1y=r0*np.exp(-p['k']))
    for s in 'BC':
        rel = curves[s]/np.where(curves['A'] > 0, curves['A'], np.nan) - 1
        out[f'maxgap_{s}'] = np.nanmax(np.abs(rel[window]))
        out[f'maxgap_at_{s}'] = t[window][np.nanargmax(np.abs(rel[window]))]
        for yr in report_years:
            out[f'gap_{s}_{yr:g}y'] = rel[idx[yr]]

    for lam in lams:
        for s, kw in dict(A=dict(structure='exact'), C=dict(structure='two-state', k_latent=k_refit)).items():
            w = ever_asy(p, years=reinf_years, lam=lam, deaths=True, **kw)
            wo = ever_asy(p, years=reinf_years, lam=lam, deaths=True, latent_reinf=False, **kw)
            for yr, a in zip(reinf_years, w - wo):
                out[f'added_{s}_lam{lam:g}_{yr}y'] = a
        for yr in reinf_years:
            out[f'added_CminusA_lam{lam:g}_{yr}y'] = out[f'added_C_lam{lam:g}_{yr}y'] - out[f'added_A_lam{lam:g}_{yr}y']
    return out, curves


def summarise(df, w):
    """ Weighted median, 10-90% interval and max over samples, and the posterior weight above each gap threshold """
    cols = [c for c in df.columns if c not in ['rand_seed', 'sample_weight', 'n_draws']]
    summ = {c: dict(median=float(wquantile(df[c], w, 0.5)), q10=float(wquantile(df[c], w, 0.1)), q90=float(wquantile(df[c], w, 0.9)),
                    min=float(df[c].min()), max=float(df[c].max())) for c in cols}
    exceed = {s: {f'>{th:g}': float(w[df[f'maxgap_{s}'] > th].sum()/w.sum()) for th in thresholds} for s in 'BC'}
    return dict(n_samples=len(df), stats=summ, weight_above=exceed, row0=df.iloc[0].to_dict())


def plot(df, w, curves):
    """ Largest gap against r0/k; C's gap over time for every sample; added progression from reinfection, C - A """
    blue, orange, ink, muted = '#2a78d6', '#eb6834', '#0b0b0b', '#52514e'
    size = 8 + 400*w/w.max()
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8), facecolor='#fcfcfb')
    for ax in axes:
        ax.set_facecolor('#fcfcfb')
        for side in ['top', 'right']:
            ax.spines[side].set_visible(False)
        ax.grid(alpha=0.25, lw=0.6)

    ax = axes[0]
    ax.scatter(df.r0_k, 100*df.maxgap_B, s=size, color=orange, alpha=0.55, lw=0, label='B: k_latent = k')
    ax.scatter(df.r0_k, 100*df.maxgap_C, s=size, color=blue, alpha=0.7, lw=0, label='C: k_latent refitted')
    ax.axhline(1, color=muted, lw=1, ls='--')
    ax.scatter([df.r0_k.iloc[0]], [100*df.maxgap_C.iloc[0]], s=80, facecolor='none', edgecolor=ink, lw=1.2, label='Set B (row 0)')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('r₀/k')
    ax.set_ylabel('Largest |two-state ÷ exact − 1|, 3 months–50 years (%)')
    ax.set_title('Largest gap from the exact model (point area ∝ weight)', fontsize=10, color=ink)
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1]
    t = np.linspace(0, tmax, tmax*per_year + 1)
    win = t >= gap_from
    for i, c in enumerate(curves):
        rel = c['C'][win]/c['A'][win] - 1
        ax.plot(t[win], 100*rel, color=blue, lw=0.8, alpha=float(0.08 + 0.9*w[i]/w.max()))
    ax.axhline(0, color=muted, lw=0.8)
    ax.set_xscale('log')
    ax.set_xlabel('Years since infection')
    ax.set_ylabel('C ÷ exact − 1, ever reached asymptomatic TB (%)')
    ax.set_title('C against the exact model, every sample (opacity ∝ weight)', fontsize=10, color=ink)

    ax = axes[2]
    for lam, col, lab in [(lams[0], '#6da7ec', f'λ = {lams[0]:g}/y'), (lams[1], '#184f95', f'λ = {lams[1]:g}/y')]: # C at both: light and dark blue
        ax.scatter(df.r0_k, 1e4*df[f'added_CminusA_lam{lam:g}_30y'], s=size, color=col, alpha=0.6, lw=0, label=lab)
    ax.axhline(0, color=muted, lw=0.8)
    ax.set_xscale('log')
    ax.set_xlabel('r₀/k')
    ax.set_ylabel('Progression added by latent reinfection by 30 y,\nC − exact (per 10,000 latent)')
    ax.set_title('Reinfection rule: C against the clock reset', fontsize=10, color=ink)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / 'posterior_progression.png', dpi=150, facecolor=fig.get_facecolor())
    return


if __name__ == '__main__' and 'plot' in sys.argv[1:]: # Redraw the figure from the saved outputs: run.py plot
    saved = sc.load(outdir / 'outputs.obj')
    plot(saved['df'], saved['df'].sample_weight.values, saved['curves'])
elif __name__ == '__main__':
    post = pd.read_csv(samples_csv)
    assert np.allclose(post.k_asy, post.k_non), 'k_asy and k_non differ in some samples'
    rows, curves = [], []
    for i, row in post.iterrows():
        out, c = run_sample(row)
        rows.append(dict(sample=i, sample_weight=row.sample_weight, n_draws=row.n_draws, rand_seed=row.rand_seed, **out))
        curves.append(c)
    df = pd.DataFrame(rows)
    w = df.sample_weight.values
    df.to_csv(outdir / 'samples.csv', index=False)
    summ = summarise(df.drop(columns=['sample']), w)
    summ['settings'] = dict(samples_csv=str(samples_csv), mu=mu, lams=lams, reinf_years=reinf_years, tmax=tmax, per_year=per_year, gap_from=gap_from, trunc=trunc)
    (outdir / 'summary.json').write_text(json.dumps(summ, indent=2, default=float) + '\n')
    sc.save(outdir / 'outputs.obj', dict(df=df, curves=curves, t=np.linspace(0, tmax, tmax*per_year + 1)))
    plot(df, w, curves)
    tbode4.save_versions(outdir, 'v4')

    st = summ['stats']
    fmt = lambda c, s=1: f"{s*st[c]['median']:.3g} [{s*st[c]['q10']:.3g}, {s*st[c]['q90']:.3g}], max {s*st[c]['max']:.3g}"
    print(f'{len(df)} samples. Weighted median [10%, 90%], max over samples:')
    for s in 'BC':
        print(f'  Largest gap {s} (%):   {fmt(f"maxgap_{s}", 100)};  weight above 1/2/5%: ' + ', '.join(f'{v:.3f}' for v in summ['weight_above'][s].values()))
        for yr in report_years:
            print(f'    gap {s} at {yr:g} y (%): {fmt(f"gap_{s}_{yr:g}y", 100)}')
    print(f'  k_latent/k (C):     {fmt("refit_ratio")}')
    print(f'  Ratio A / C:        {fmt("Ratio_A")} / {fmt("Ratio_C")}')
    print(f'  Classes at 1e-6 / 1e-10: {fmt("classes_1e6")} / {fmt("classes_1e10")}')
    print(f'  Exact model, progression between years 1 and 10: {fmt("B1_to_10_A")}; share from people still latent at year 1: {fmt("share_from_latent_A")}')
    print(f'  Exact model at year 1: latent {fmt("latent_at_1y_A")}, non-infectious {fmt("non_at_1y_A")}; latent hazard r0·exp(−k) {fmt("latent_hazard_1y")}/y')
    for lam in lams:
        for yr in reinf_years:
            print(f'  Reinfection λ={lam:g}, {yr} y: added A {fmt(f"added_A_lam{lam:g}_{yr}y")}; C − A {fmt(f"added_CminusA_lam{lam:g}_{yr}y")}')
    r0 = df.iloc[0]
    print(f'Row 0 (set B): maxgap C {100*r0.maxgap_C:.3f}%, k_latent {r0.k_latent_refit:.4f}, C − A by 10 y at λ=0.12: {r0["added_CminusA_lam0.12_10y"]:.2e}')
