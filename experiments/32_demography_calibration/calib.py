"""
Experiment 32: the model per draw, targets and IMIS for recalibrating the ODE with ages and India's demography. See README.md.

The model is tbode4.fastode_age.run_cohort() at monthly steps with the mid-step force of infection (CONVENTIONS.md 8). Each draw starts at
its untreated equilibrium without ages, settles for 100 years at 1950's demography (frozen), then runs 1950–2041 with WPP's demography and
exp 24's treatment scale-up. Priors, parameterisations, target means and the likelihood are exp 24's; every target built on the prevalence
survey uses 15+ prevalence, as tb_LAI_TPT computes them.
"""
import os
for _k in ['NUMBA_NUM_THREADS', 'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']:
    os.environ.setdefault(_k, '1') # One thread per process: the workers are single-threaded; before importing numpy and starsim
import importlib.util
import multiprocessing as mp
import numpy as np
import sciris as sc
from scipy.stats import multivariate_normal
import tbode4
from tbode4 import fastode as fo, fastode_age as fa

__all__ = ['priors', 'names', 'targets', 'out_keys', 'simulate', 'run_draw', 'loglik', 'evaluate', 'imis']

here = sc.thispath(__file__)
root = here.parents[1]
_spec = importlib.util.spec_from_file_location('calib24', root / 'experiments' / '24_treatment_calibration' / 'calib.py')
cb24 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cb24)

priors, names, d = cb24.priors, cb24.names, cb24.d
from_unit, to_unit, to_pars, loglik = cb24.from_unit, cb24.to_unit, cb24.to_pars, cb24.loglik
targets = cb24.targets # Means, sigmas and survey sizes unchanged; the model's prevalence outputs are 15+
extra = ['prevalence_all_2020', 'prevalence_15plus_1995', 'prevalence_15plus_2024', 'incidence_2024', 'notifications_2020', 'mortality_2024',
         'population_2020'] # Reported: per person (per year for flows); population relative to 1950
out_keys = list(targets) + extra + ['error']
h = 1/12 # Monthly steps (and age groups one month wide)
start, stop = 1950, 2041
settle_years = 100
ncpus = int(os.environ.get('EXP32_NCPUS', tbode4.ncpus)) # 160 on zebra; the laptop default otherwise


#%% The model per draw

def untreated_equilibrium(nat):
    """ The start without ages (fastode's untreated equilibrium, births replacing deaths at 1/70, as exp 27 and 31): shares, and the share of NON already asymptomatic """
    x = fo.equilibrium(fo.make_pars(**nat))
    if x is None:
        return None, None
    s = dict(zip(fo.names(fo.make_pars(**nat)), x))
    init = {k: s[k] for k in fo.state_names if k not in ['non0', 'non1']}
    init['non_infectious'] = s['non0'] + s['non1']
    return init, s['non1']/(s['non0'] + s['non1'])


def model_pars(nat, trend):
    """ fastode_age's parameters: the natural history with mu = 0 (WPP's deaths instead), and exp 24's schedules (sym_dead among them) """
    nat = {k: v for k, v in nat.items() if k not in ['mu', 'sym_dead']}
    return fo.make_pars(**nat, mu=0, **fa.schedules(trend))


def run_draw(theta):
    """ One draw's run, 1950–2041 (per-step outputs), or None if it has no endemic equilibrium without ages or loses TB in the frozen start """
    nat, trend = to_pars(theta)
    init, flagged = untreated_equilibrium(nat)
    if init is None:
        return None
    p = model_pars(nat, trend)
    st = fa.settle(p, h, init, flagged, years=settle_years, start=start)
    X = st.X
    if (X[fa.A] + X[fa.Y]).sum()/X[:st.layout.ns].sum() < 1e-9:
        return None
    table, flagged = fa.start_table(st)
    out = fa.run_cohort(p, h, start, stop, table=table, share_flagged=flagged)
    out.trend = trend
    out.drift = st.drift
    return out


def outputs(out):
    """ The targets (exp 24's, with 15+ prevalence wherever the survey enters) and reported outputs, from a run's per-step outputs """
    year = out.year
    at = lambda y: int(np.argmin(np.abs(year - y)))
    total = lambda k, y: out[k][(year >= y - 1e-9) & (year < y + 1 - 1e-9)].sum() # Per person over the calendar year (people at each step's start)
    prev15 = lambda y: out['prev_15+'][at(y)]
    inc = {y: total('first', y) for y in [2015, 2020, 2024]}
    deaths = {y: total('deaths', y) for y in [2015, 2024]}
    notif = (1 + out.trend.over_detection)*total('tx_start', 2020)
    return dict(
        prevalence        = prev15(2020),
        prop_asymptomatic = out['prop_asymptomatic_15+'][at(2020)],
        prev_to_inc       = prev15(2020)/inc[2020],
        mort_to_inc       = deaths[2024]/inc[2024],
        prev_to_notif     = prev15(2020)/notif if notif > 0 else np.nan,
        inc_ratio         = inc[2024]/inc[2015],
        mort_ratio        = deaths[2024]/deaths[2015],
        prev_2040_ratio   = prev15(2040)/prev15(2020),
        prevalence_all_2020    = out.prevalence[at(2020)],
        prevalence_15plus_1995 = prev15(1995),
        prevalence_15plus_2024 = prev15(2024),
        incidence_2024    = inc[2024],
        notifications_2020= notif,
        mortality_2024    = deaths[2024],
        population_2020   = out.alive[at(2020)],
    )


def progression(nat):
    """ The progression targets from tbode4.latent's exact chain, as exp 24 """
    chain = {k: nat[k] for k in ['inf_cle', 'inf_non', 'inf_asy', 'k_latent', 'inf_non_remote', 'inf_asy_remote', 'non_rec', 'non_asy']}
    return tbode4.progression_targets(chain, 'two-state')


def simulate(theta):
    """ The outputs at one parameter vector, as a dict over out_keys: NaN for every model output if there's no endemic state; 'error' 1 if the run raised """
    out = {k: np.nan for k in out_keys}
    out['error'] = 0.0
    nat, _ = to_pars(theta)
    try:
        r = run_draw(theta)
        if r is not None:
            out.update(outputs(r))
    except Exception as e: # Counted and reported in the summary, not silently dropped
        print(f'simulate error at {list(np.round(theta, 6))}: {e!r}', flush=True)
        out['error'] = 1.0
    out.update(progression(nat))
    return out


def trajectory(theta):
    """ One draw by calendar year (fastode_age.calendar(): stocks on 1 January, flows over the year, per 100,000), with notifications as
    reported ((1 + over_detection) × treatment starts); None if no endemic state """
    r = run_draw(theta)
    if r is None:
        return None
    cal = fa.calendar(r)
    cal['notifications'] *= 1 + r.trend.over_detection
    return cal


#%% Parallel evaluation: a persistent pool of single-threaded workers

_pool = None


def pool():
    """ A fork pool of ncpus workers, made once; the parent warms WPP's data (loaded on first use) before forking, so workers inherit it """
    global _pool
    if _pool is None:
        fa.Demography(fa.cohort_grid(h), cohort=True)
        _pool = mp.get_context('fork').Pool(ncpus)
    return _pool


def evaluate(Z, tgts=None):
    """ Outputs (array [point, out_keys]) and log-likelihoods at unit-cube points Z (rows) """
    theta = from_unit(Z)
    outs = pool().map(simulate, list(theta), chunksize=1) if len(theta) > 1 else [simulate(t) for t in theta]
    ll = np.array([loglik(o, tgts) for o in outs])
    return np.array([[o[k] for k in out_keys] for o in outs]), ll


def imis(tgts=None, n0=10_000, b=1_000, j=3_000, max_iter=300, seed=1, verbose=True):
    """
    Incremental mixture importance sampling (Raftery and Bao 2010) on the unit cube, where the prior is uniform: exp 24's algorithm and settings.

    Starts with n0 prior draws; each iteration adds b draws from a multivariate normal centred on the current highest-weight point, with the
    weighted covariance of its b nearest points (Mahalanobis distance under the prior's covariance); stops when the expected number of unique
    points in a resample of j reaches (1 − 1/e)·j. Importance weights are likelihood × prior ÷ the mixture of the prior and every component.

    Unlike exp 24's, the sum of the components' densities at each point is kept and updated as components and points are added, rather than
    recomputed for every point and component each iteration. At the end it is recomputed in full and the largest relative difference stored.
    """
    rng = np.random.default_rng(seed)
    Z = rng.random((n0, d))
    outs, ll = evaluate(Z, tgts)
    comps, dens = [], [] # Components (mean, covariance) and their frozen densities
    csum = np.zeros(n0) # Sum over components of each component's density at each point
    log = []
    prior_cov_inv = np.eye(d)*12.0 # Covariance of a uniform on [0, 1] is 1/12
    T = sc.timer()
    for it in range(max_iter + 1):
        n = len(Z)
        q = n0/n + (b/n)*csum
        inside = np.all((Z >= 0) & (Z <= 1), axis=1)
        logw = np.where(inside & np.isfinite(ll), ll - np.log(q), -np.inf)
        w = np.exp(logw - np.max(logw))
        w /= w.sum()
        ess = 1/np.sum(w**2)
        unique = np.sum(1 - (1 - w)**j)
        log.append(dict(iteration=it, n=n, ess=ess, expected_unique=unique, max_weight=w.max(), max_ll=float(np.max(ll)), time=T.toc(output=True)))
        if verbose:
            print(f'IMIS {it:3d}: n = {n:7,}, ESS = {ess:8.1f}, expected unique in {j} = {unique:7.1f} (target {(1 - np.exp(-1))*j:.0f}), '
                  f'max weight {w.max():.4f}, max ll {np.max(ll):.3f}, {T.toc(output=True):.0f} s', flush=True)
        if unique >= (1 - np.exp(-1))*j or it == max_iter:
            break
        i = int(np.argmax(w))
        diff = Z - Z[i]
        dist = np.einsum('ij,jk,ik->i', diff, prior_cov_inv, diff)
        near = np.argsort(dist)[:b]
        wn = w[near] + 1/n
        wn /= wn.sum()
        cov = np.cov(Z[near].T, aweights=wn) + 1e-10*np.eye(d)
        comps.append((Z[i].copy(), cov))
        dens.append(multivariate_normal(Z[i], cov, allow_singular=True))
        csum += dens[-1].pdf(Z)
        Znew = rng.multivariate_normal(Z[i], cov, size=b)
        new_outs = np.full((b, len(out_keys)), np.nan)
        new_ll = np.full(b, -np.inf)
        ok = np.all((Znew >= 0) & (Znew <= 1), axis=1)
        if ok.any():
            new_outs[ok], new_ll[ok] = evaluate(Znew[ok], tgts)
        Z = np.vstack([Z, Znew])
        outs = np.vstack([outs, new_outs])
        ll = np.concatenate([ll, new_ll])
        csum = np.concatenate([csum, np.sum([dn.pdf(Znew) for dn in dens], axis=0)])
    full = np.sum([dn.pdf(Z) for dn in dens], axis=0) if dens else np.zeros(len(Z))
    check = float(np.max(np.abs(full - csum)/np.maximum(full, 1e-300))) if dens else 0.0
    if verbose:
        print(f'IMIS mixture density, running against recomputed: largest relative difference {check:.1e}', flush=True)
    return sc.objdict(Z=Z, theta=from_unit(Z), outs=outs, out_keys=out_keys, ll=ll, weights=w, comps=comps, log=log, n0=n0, b=b, j=j, seed=seed,
                      mixture_check=check)
