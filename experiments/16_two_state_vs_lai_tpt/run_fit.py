"""
Exp 16, v4 side: fit the two-state latent model (remote progression = phi x recent) to tb_LAI_TPT's progression target means within its
priors, for all 345 posterior samples, phi on a grid, and two arms (1: cascade fixed; 2: inf_cle, non_rec, non_asy freed within their
priors). Then compute the cohort outputs for the 12 face-validity parameter sets, for the fits and for the exact r0·exp(−k·τ) reference,
and check the chain against v4's own ODE. See README.md.

v4 environment. Run from the repo root: v4/.venv/bin/python experiments/16_two_state_vs_lai_tpt/run_fit.py
"""
import os
from pathlib import Path
import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss
from scipy.linalg import expm
from scipy.optimize import least_squares
import tbode4

outdir = Path(__file__).resolve().parent
lai_dir = outdir.parents[2] / 'tb_LAI_TPT'
samples_csv = lai_dir / 'results/hm_real_20260728/trajectory_selection/collapsed_samples_eps1.2.csv'
calib_csv = lai_dir / 'results/hm_real_20260728/trajectory_selection/scored_points.csv' # Every scored point; posterior_samples.csv lacks 115 of the 345

phis = [0.0, 0.003, 0.01, 0.03, 0.1, 0.3]
arms = [1, 2]
n_sets = 12
prior = dict(inf_non=(0.04, 0.45), inf_asy=(0.01, 0.60), k_latent=(0.0, 12.0), inf_cle=(0.5, 1.5), non_rec=(0.14, 0.23), non_asy=(0.21, 0.28),
             inf_non_remote=0.23, inf_asy_remote=0.20) # tb_LAI_TPT's priors, carried over as agreed (user, 2026-11-04)
freed = ['inf_cle', 'non_rec', 'non_asy'] # Arm 2
penalty = 0.01 # Arm 2's tie-break: weight on each freed parameter's log-distance from the sample's value, per log prior width
tmax, per_year = 51, 12 # Monthly grid, matching the v0.10 cohort (index i + 1 here is the end of its month i + 1)
cascade_keys = ['inf_cle', 'non_rec', 'non_asy', 'asy_non', 'asy_sym', 'sym_asy']


def sample_pars(row):
    """ A posterior sample's natural history (rates per year), with sym_dead = 0 as in tb_LAI_TPT's cohort """
    p = {k: float(row[k]) for k in cascade_keys}
    p.update(inf_non=float(row.inf_non), inf_asy=float(row.inf_asy), k=float(row.k_asy), pi=float(row.rr_reinfection_rec), sym_dead=0.0)
    return p


def two_state(base, f_non, phi, r0, k_latent, **cascade):
    """ Two-state parameters: recent rates from r0 and the NON : ASY split, remote = phi x recent (capped by the priors) """
    p = dict(base, **cascade, inf_non=f_non*r0, inf_asy=(1 - f_non)*r0, k_latent=k_latent)
    p.update(inf_non_remote=min(phi*p['inf_non'], prior['inf_non_remote']), inf_asy_remote=min(phi*p['inf_asy'], prior['inf_asy_remote']))
    return p


def zscores(vals):
    """ Log-space z-scores of ProgressionY1 and ProgressionRatio against tb_LAI_TPT's target means """
    return np.array([np.log(vals[k]/tbode4.target_means[k])/tbode4.target_sigmas[k] for k in ['ProgressionY1', 'ProgressionRatio']])


def fit(base, phi, arm):
    """ Fit the two-state model to the target means within the priors (best of several starts); returns (pars, z, x) """
    f_non = base['inf_non']/(base['inf_non'] + base['inf_asy'])
    lo_r0 = max(prior['inf_non'][0]/f_non, prior['inf_asy'][0]/(1 - f_non))
    hi_r0 = min(prior['inf_non'][1]/f_non, prior['inf_asy'][1]/(1 - f_non))
    names = ['r0', 'k_latent'] + (freed if arm == 2 else [])
    lo = [lo_r0, prior['k_latent'][0]] + [prior[k][0] for k in freed if arm == 2]
    hi = [hi_r0, prior['k_latent'][1]] + [prior[k][1] for k in freed if arm == 2]

    def pars(x):
        kw = dict(zip(names, x))
        return two_state(base, f_non, phi, kw.pop('r0'), kw.pop('k_latent'), **kw)

    def resid(x):
        z = zscores(tbode4.progression_targets(pars(x), 'two-state'))
        if arm == 2:
            reg = [penalty*np.log(v/base[k])/np.log(prior[k][1]/prior[k][0]) for k, v in zip(freed, x[2:])]
            return np.concatenate([z, reg])
        return z

    r0_base = np.clip(base['inf_non'] + base['inf_asy'], lo_r0, hi_r0)
    tail = [np.clip(base[k], *prior[k]) for k in freed] if arm == 2 else []
    starts = [[r0_base, np.clip(0.95*base['k'], 0.1, 11.9)] + tail, [r0_base, 11.0] + tail, [0.5*(lo_r0 + hi_r0), 3.0] + tail]
    best = None
    for x0 in starts:
        res = least_squares(resid, x0, bounds=(lo, hi), xtol=1e-12, ftol=1e-12, gtol=1e-12)
        if best is None or res.cost < best.cost:
            best = res
    p = pars(best.x)
    return p, zscores(tbode4.progression_targets(p, 'two-state')), dict(zip(names, best.x))


def fit_sample(i, row):
    """ All fits for one posterior sample """
    base = sample_pars(row)
    z_exact = zscores(tbode4.progression_targets(base, 'exact'))
    out = []
    for arm in arms:
        for phi in phis:
            p, z, x = fit(base, phi, arm)
            rec = dict(sample=i, sample_weight=row.sample_weight, arm=arm, phi=phi, z_Y1=z[0], z_Ratio=z[1], zmax=np.abs(z).max(),
                       z_Y1_exact=z_exact[0], z_Ratio_exact=z_exact[1], **{k: p[k] for k in ['inf_non', 'inf_asy', 'inf_non_remote', 'inf_asy_remote', 'k_latent'] + freed},
                       **{f'{k}_sample': base[k] for k in freed}, k=base['k'])
            out.append(rec)
    return out


def later_split(p, structure):
    """ Progression to ASY between years 1 and 10 (share of the cohort), and the share of it from people still latent at year 1 """
    Q, x0, idx = tbode4.latent_chain(p, structure)
    x1, P9 = x0 @ expm(Q), expm(Q*9)
    from_latent, from_non = x1[idx.latent] @ P9[idx.latent, idx.ASY], x1[idx.NON]*P9[idx.NON, idx.ASY]
    return from_latent + from_non, from_latent/(from_latent + from_non)


def dynamics(p, structure):
    """ The cohort outputs for one parameter set and structure, on the monthly grid """
    t, A = tbode4.first_passage(p, structure, 'A', tmax=tmax, per_year=per_year)
    _, B = tbode4.first_passage(p, structure, 'B', tmax=tmax, per_year=per_year)
    _, C = tbode4.first_passage(p, structure, 'C', tmax=tmax, per_year=per_year)
    _, ex = tbode4.latent_exits(p, structure, tmax=tmax, per_year=per_year)
    _, X, idx = tbode4.cohort_curve(p, structure, tmax=tmax, per_year=per_year, cascade=True)
    latent = X[:, idx.latent].sum(axis=1)
    later, share = later_split(p, structure)
    out = sc.objdict(t=t, A=A, B=B, C=C, exit_cle=ex.cleared, exit_non=ex.non, exit_asy=ex.asy, latent=latent, hazard=tbode4.latent_hazard(X, idx),
                     later=later, share_from_latent=share)
    if structure == 'two-state':
        out.recent, out.remote = X[:, idx.recent], X[:, idx.remote]
    return out


def v4_ode_check(p):
    """ v4's own tbsim.TB(solver='ode') cohort against the chain: the largest difference in any state's share over 50 years """
    tb_kw = {k: p[k] for k in cascade_keys + ['inf_non', 'inf_asy', 'k_latent', 'inf_non_remote', 'inf_asy_remote', 'sym_dead']}
    sim = tbode4.make_sim(beta=0.0, solver='ode', dt=ss.months(1), start=2000, stop=2050, init=dict(infection=1.0), mu=0.0, replace_tb_deaths=True, **tb_kw)
    sim.run()
    r = sim.results.tb
    t = np.asarray(sim.t.yearvec, dtype=float) - 2000
    _, X, idx = tbode4.cohort_curve(p, 'two-state', tmax=50, per_year=12, cascade=True)
    tc = np.linspace(0, 50, 50*12 + 1)
    pairs = dict(infection=idx.recent, infection_remote=idx.remote, non_infectious=idx.NON, asymptomatic=idx.ASY, symptomatic=idx.SYM, cleared=idx.CLE, recovered=idx.REC)
    N = sum(np.asarray(r[s].values, dtype=float) for s in pairs)
    diff = {s: float(np.max(np.abs(np.asarray(r[s].values, dtype=float)/N - np.interp(t, tc, X[:, j])))) for s, j in pairs.items()}
    return diff


if __name__ == '__main__':
    post = pd.read_csv(samples_csv)
    assert np.allclose(post.k_asy, post.k_non)

    # Question 1: fits for all samples
    T = sc.timer()
    res = sc.parallelize(fit_sample, iterkwargs=[dict(i=i, row=row) for i, row in post.iterrows()], ncpus=tbode4.ncpus)
    fits = pd.DataFrame([r for rows in res for r in rows])
    T.toc('Fits')

    # tb_LAI_TPT as run, for question 1: the calibration's own saved progression values (matched by parameters)
    calib = pd.read_csv(calib_csv)
    calib['inf_cle'] = np.exp(calib.log_inf_cle)
    zc = []
    for i, row in post.iterrows():
        d = (calib.inf_non - row.inf_non).abs() + (calib.k_prog - row.k_asy).abs() + (calib.inf_cle - row.inf_cle).abs()
        j = d.idxmin()
        assert d[j] < 1e-9, f'sample {i} not found in {calib_csv.name}'
        zc.append(zscores(dict(ProgressionY1=calib.mod__ProgressionY1[j], ProgressionRatio=calib.mod__ProgressionRatio[j])))
    zc = np.array(zc)
    fits = fits.merge(pd.DataFrame(dict(sample=np.arange(len(post)), z_Y1_calib=zc[:, 0], z_Ratio_calib=zc[:, 1])), on='sample')
    fits.to_csv(outdir / 'fits.csv', index=False)

    # Question 2: cohort outputs for the 12 face-validity sets
    curves = {}
    for i in range(n_sets):
        base = sample_pars(post.iloc[i])
        curves[(i, 'exact')] = dynamics(base, 'exact')
        for arm in arms:
            for phi in phis:
                f = fits[(fits['sample'] == i) & (fits.arm == arm) & (fits.phi == phi)].iloc[0]
                p = dict(base, **{k: f[k] for k in ['inf_non', 'inf_asy', 'inf_non_remote', 'inf_asy_remote', 'k_latent'] + freed})
                curves[(i, arm, phi)] = dynamics(p, 'two-state')
                curves[(i, arm, phi)].pars = p

    # Check: the chain is v4's model (param00, arm 1, phi = 0.03)
    f = fits[(fits['sample'] == 0) & (fits.arm == 1) & (fits.phi == 0.03)].iloc[0]
    pchk = dict(sample_pars(post.iloc[0]), **{k: f[k] for k in ['inf_non', 'inf_asy', 'inf_non_remote', 'inf_asy_remote', 'k_latent']})
    check = v4_ode_check(pchk)
    print('v4 ODE vs chain, largest difference in each state\'s share (param00, arm 1, phi = 0.03):', {k: f'{v:.1e}' for k, v in check.items()})

    sc.save(outdir / 'curves.obj', dict(curves=curves, phis=phis, arms=arms, n_sets=n_sets, weights=post.head(n_sets).sample_weight.values, v4_ode_check=check,
                                        settings=dict(prior=prior, freed=freed, penalty=penalty, tmax=tmax, per_year=per_year)))
    tbode4.save_versions(outdir, 'v4')
    os.replace(outdir / 'versions.json', outdir / 'versions_v4.json')

    w = fits.sample_weight
    for arm in arms:
        for phi in phis:
            d = fits[(fits.arm == arm) & (fits.phi == phi)]
            print(f'arm {arm}, phi = {phi:<5g}: weight with both |z| <= 1: {d.sample_weight[d.zmax <= 1].sum():.3f}, <= 2: {d.sample_weight[d.zmax <= 2].sum():.3f}')
    d = fits[(fits.arm == 1) & (fits.phi == 0)]
    zcm = np.maximum(d.z_Y1_calib.abs(), d.z_Ratio_calib.abs())
    print(f'tb_LAI_TPT as run (calibration values): weight with both |z| <= 1: {d.sample_weight[zcm <= 1].sum():.3f}, <= 2: {d.sample_weight[zcm <= 2].sum():.3f}')
