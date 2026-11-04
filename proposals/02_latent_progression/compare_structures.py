"""
Latent-progression structures compared on a cohort, with set B's natural history (tb_LAI_TPT's best fit, k = 5.846). Deterministic:
each structure is a small continuous-time Markov chain for one person, solved with expm. Numbers for options.md.

Structures:
- exact:      tb_LAI_TPT / v0.10's hazard r·exp(−k·τ) with clock reset on reinfection, as Poisson-mixture classes (exact, 12 classes)
- two-state:  v4's k_latent = k, no progression from the remote stage (Cliff's Bremote); reinfection moves remote → recent
- two-state, refit: the same, with k_latent chosen so that ProgressionY1 matches the exact model

Scenarios:
- cohort: everyone latent at t = 0, no reinfection or deaths, asymptomatic TB absorbing (tb_LAI_TPT's ProgressionY1 / ProgressionRatio)
- reinfection: the same cohort with background deaths (1/70 per year) and a constant force of infection lam: cleared and recovered
  people are reinfected at pi·lam (new infections), and latent people at pi·lam (the structure's reinfection rule). Run with and
  without latent reinfection, to separate what each structure's reinfection rule adds.

Run with either environment: .venv/bin/python proposals/02_latent_progression/compare_structures.py
"""
import json
import numpy as np
from pathlib import Path
from scipy.linalg import expm
from scipy.optimize import brentq
from scipy.stats import poisson

outdir = Path(__file__).resolve().parent
k = 5.845726750296858
B = dict(inf_cle=0.5365965386150265, inf_non=0.4103330692741019, inf_asy=0.2741560588304815, non_rec=0.1706711906227213,
         non_asy=0.2426759489829796, pi=0.2415635040650483)
lam_B = 24.79881027584826*(0.9407760328175468*0.002792 + 0.002208) # Set B's force of infection at Bk's equilibrium (exp 13 calibration), per year
mu = 1/70


def chain(structure, lam=0.0, latent_reinf=True, deaths=False, k_latent=k, n_classes=12):
    """
    The generator Q and initial vector of one person's chain. States: the structure's latent states, then CLE, NON, REC, ASY (absorbing,
    'ever reached asymptomatic TB'), DEAD (absorbing). New infections (and the cohort's start) enter the latent states with `entry`
    """
    p = B
    r0 = p['inf_non'] + p['inf_asy']
    f_non = p['inf_non']/r0
    if structure == 'exact':
        w = poisson.pmf(np.arange(n_classes), r0/k)
        w[-1] += 1 - w.sum()
        prog = np.arange(n_classes)*k # Class m progresses at m·k
        nL = n_classes
        entry = w
    else: # recent, remote
        nL = 2
        prog = np.array([r0, 0.0])
        entry = np.array([1.0, 0.0])
    CLE, NON, REC, ASY, DEAD = nL, nL + 1, nL + 2, nL + 3, nL + 4
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


def ever_asy(structure, years, **kwargs):
    Q, x0, ASY = chain(structure, **kwargs)
    return np.array([(x0 @ expm(Q*t))[ASY] for t in years])


def cohort_metrics(structure, **kwargs):
    years = [0.25, 0.5, 1, 2, 5, 10, 50]
    b = ever_asy(structure, years, **kwargs)
    out = {f'B({t:g}y)': float(v) for t, v in zip(years, b)}
    out['ProgressionY1'] = float(b[2])
    out['ProgressionRatio'] = float(b[2]/b[5])
    return out


if __name__ == '__main__':
    exact = cohort_metrics('exact')
    k_refit = brentq(lambda kl: cohort_metrics('two-state', k_latent=kl)['ProgressionY1'] - exact['ProgressionY1'], 1, 20)
    structures = {'exact': dict(structure='exact'), 'two-state': dict(structure='two-state'), 'two-state, refit': dict(structure='two-state', k_latent=k_refit)}

    res = dict(k=k, k_latent_refit=k_refit, lam_B=lam_B, targets=dict(ProgressionY1=0.0523, ProgressionRatio=0.663), cohort={}, reinfection={})
    print(f'Cohort (no reinfection, no deaths); k_latent refit to the exact ProgressionY1: {k_refit:.3f}')
    for lab, kw in structures.items():
        res['cohort'][lab] = cohort_metrics(**kw)
        print(f'  {lab:18s} ' + ', '.join(f'{key} {v:.4f}' for key, v in res['cohort'][lab].items()))

    years = [10, 30]
    for mult in [1, 5]:
        lam = lam_B*mult
        tag = f'{mult}x lam_B ({lam:.3f}/y)'
        res['reinfection'][tag] = {}
        print(f'Reinfection at {tag}, with background deaths: ever asymptomatic by {years} years')
        for lab, kw in structures.items():
            with_ = ever_asy(years=years, lam=lam, deaths=True, **kw)
            without = ever_asy(years=years, lam=lam, deaths=True, latent_reinf=False, **kw)
            res['reinfection'][tag][lab] = dict(with_latent_reinfection=with_.tolist(), without=without.tolist(), added=(with_ - without).tolist())
            print(f'  {lab:18s} with {with_.round(4).tolist()}, without {without.round(4).tolist()}, added by latent reinfection {(with_ - without).round(4).tolist()}')
    (outdir / 'compare_structures.json').write_text(json.dumps(res, indent=2) + '\n')
