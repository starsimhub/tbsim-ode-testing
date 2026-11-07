"""
Experiment 20: priors, model, likelihood and IMIS for calibrating the v4 ODE (two-state latent model, homogeneous mixing, no
treatment) to tb_LAI_TPT's six applicable targets. See README.md.

The model is tbode4.fastode (a standalone copy of v4's ODE, checked by v4/tests/test_fastode.py) at its endemic equilibrium, with
incidence counted as the first entry to asymptomatic TB per episode; the progression targets come from tbode4.latent's exact chain
(a no-reinfection latent cohort, as exps 14–17).
"""
import numpy as np
import sciris as sc
from scipy.stats import binom, multivariate_normal
import tbode4
from tbode4 import fastode as fo

__all__ = ['priors', 'names', 'targets', 'to_pars', 'simulate', 'loglik', 'imis']

# Priors: (low, high, log-uniform?). tb_LAI_TPT's box (run_history_matching.parameter_bounds), without its treatment, household and
# trend parameters, with its time-since-infection decay replaced by the two-state parameters (README)
priors = sc.objdict(
    beta        = (2.0, 100.0, True),
    trans_asymp = (0.5, 1.0, False),
    pi          = (0.2, 1.0, False),  # Every reinfection risk: cleared, recovered, treated, and remote latent -> recent
    inf_cle     = (0.5, 1.5, True),
    inf_non     = (0.04, 0.45, False),
    inf_asy     = (0.01, 0.60, True),
    k_latent    = (0.1, 12.0, False),
    phi         = (0.0, 0.1, False),  # Remote progression ÷ recent
    non_rec     = (0.14, 0.23, False),
    non_asy     = (0.21, 0.28, False),
    asy_non     = (0.70, 2.03, False),
    asy_sym     = (0.56, 0.99, False),
    sym_asy     = (0.46, 0.85, False),
    sym_dead    = (0.20, 0.50, False),
)
names = list(priors.keys())
d = len(names)
mu = 1/70 # Background mortality, with births replacing every death (tbsim's default)

# The targets (README): Binomial for the survey proportions, log-Normal for the rest (mean, log-space sigma)
targets = sc.objdict(
    prevalence        = dict(kind='binomial', x=77, n=24815),
    prop_asymptomatic = dict(kind='binomial', x=123, n=188),
    prev_to_inc       = dict(kind='lognormal', mean=310.296/197.482, sigma=0.136),
    ProgressionY1     = dict(kind='lognormal', mean=0.0523, sigma=0.094),
    ProgressionRatio  = dict(kind='lognormal', mean=0.663, sigma=0.071),
)
# Reported but not in the likelihood: without treatment the model can't reach it (coverage check; dropped by the user, 2026-11-06,
# to come back when treatment is added)
reported = sc.objdict(
    mort_to_inc       = dict(kind='lognormal', mean=21.073/187.31, sigma=0.127),
)


def from_unit(z):
    """ Unit-cube coordinates (rows) -> parameter values; log-uniform priors are uniform in the log """
    z = np.atleast_2d(z)
    out = np.empty_like(z, dtype=float)
    for j, (lo, hi, lg) in enumerate(priors.values()):
        out[:, j] = np.exp(np.log(lo) + z[:, j]*(np.log(hi) - np.log(lo))) if lg else lo + z[:, j]*(hi - lo)
    return out


def to_unit(theta):
    theta = np.atleast_2d(theta)
    out = np.empty_like(theta, dtype=float)
    for j, (lo, hi, lg) in enumerate(priors.values()):
        out[:, j] = (np.log(theta[:, j]) - np.log(lo))/(np.log(hi) - np.log(lo)) if lg else (theta[:, j] - lo)/(hi - lo)
    return out


def to_pars(theta):
    """ One parameter vector -> tbsim/fastode parameters (rates per year) """
    t = dict(zip(names, np.asarray(theta, dtype=float)))
    return dict(beta=t['beta'], trans_asymp=t['trans_asymp'], rr_reinfection_cleared=t['pi'], rr_reinfection_rec=t['pi'],
                rr_reinfection_treat=t['pi'], rr_reinfection_inf=t['pi'], inf_cle=t['inf_cle'], inf_non=t['inf_non'], inf_asy=t['inf_asy'],
                k_latent=t['k_latent'], inf_non_remote=t['phi']*t['inf_non'], inf_asy_remote=t['phi']*t['inf_asy'], non_rec=t['non_rec'],
                non_asy=t['non_asy'], asy_non=t['asy_non'], asy_sym=t['asy_sym'], sym_asy=t['sym_asy'], sym_dead=t['sym_dead'], mu=mu)


def simulate(theta):
    """ The six target outputs at one parameter vector (NaN for every output if there is no endemic equilibrium) """
    pars = to_pars(theta)
    p = fo.make_pars(**pars)
    out = sc.objdict({k: np.nan for k in list(targets) + list(reported)})
    x = fo.equilibrium(p)
    if x is not None:
        s = fo.summary(x, p)
        out.update(prevalence=s.prevalence, prop_asymptomatic=s.prop_asymptomatic, mort_to_inc=s.mort_to_inc, prev_to_inc=s.prev_to_inc)
    chain = dict(inf_cle=pars['inf_cle'], inf_non=pars['inf_non'], inf_asy=pars['inf_asy'], k_latent=pars['k_latent'], inf_non_remote=pars['inf_non_remote'],
                 inf_asy_remote=pars['inf_asy_remote'], non_rec=pars['non_rec'], non_asy=pars['non_asy'])
    out.update(tbode4.progression_targets(chain, 'two-state'))
    return out


def loglik(out, tgts=None):
    """ The log-likelihood of one set of outputs (-inf if any is missing) """
    tgts = tgts or targets
    total = 0.0
    for k, t in tgts.items():
        v = out[k]
        if not np.isfinite(v) or v <= 0:
            return -np.inf
        if t['kind'] == 'binomial':
            total += binom.logpmf(t['x'], t['n'], min(v, 1 - 1e-12))
        else:
            total += -0.5*((np.log(t['mean']) - np.log(v))/t['sigma'])**2 - np.log(t['sigma']*np.sqrt(2*np.pi)) # The data's log mean ~ Normal(log of the model value, sigma)
    return total


def evaluate(Z, tgts=None, ncpus=None):
    """ Outputs and log-likelihoods at unit-cube points Z (rows), in parallel """
    theta = from_unit(Z)
    outs = sc.parallelize(simulate, iterarg=list(theta), ncpus=ncpus or tbode4.ncpus) if len(theta) > 8 else [simulate(t) for t in theta]
    ll = np.array([loglik(o, tgts) for o in outs])
    return outs, ll


def imis(tgts=None, n0=10_000, b=1_000, j=3_000, max_iter=150, seed=1, verbose=True):
    """
    Incremental mixture importance sampling (Raftery and Bao 2010), on the unit cube where the prior is uniform.

    Starts with n0 prior draws; each iteration adds b draws from a multivariate normal centred on the current highest-weight point,
    with the weighted covariance of its b nearest points (Mahalanobis distance under the prior's covariance); stops when the expected
    number of unique points in a resample of j reaches (1 − 1/e)·j. Importance weights are likelihood × prior ÷ the mixture of the prior
    and every component so far.

    Returns an objdict: Z (unit-cube points), theta, outs, ll, weights, the components, and a per-iteration log
    """
    rng = np.random.default_rng(seed)
    Z = rng.random((n0, d))
    outs, ll = evaluate(Z, tgts)
    comps = [] # (mean, cov)
    log = []
    prior_cov_inv = np.eye(d)*12.0 # Covariance of a uniform on [0, 1] is 1/12
    T = sc.timer()
    for it in range(max_iter + 1):
        n = len(Z)
        # Mixture density: (n0/n)·prior + (b/n)·Σ components; the prior density on the cube is 1
        q = np.full(n, n0/n)
        for m, c in comps:
            q += (b/n)*multivariate_normal(m, c, allow_singular=True).pdf(Z)
        inside = np.all((Z >= 0) & (Z <= 1), axis=1)
        logw = np.where(inside & np.isfinite(ll), ll - np.log(q), -np.inf)
        w = np.exp(logw - np.max(logw))
        w /= w.sum()
        ess = 1/np.sum(w**2)
        unique = np.sum(1 - (1 - w)**j)
        log.append(dict(iteration=it, n=n, ess=ess, expected_unique=unique, max_weight=w.max(), time=T.toc(output=True)))
        if verbose:
            print(f'IMIS {it:3d}: n = {n:7,}, ESS = {ess:8.1f}, expected unique in {j} = {unique:7.1f} (target {(1 - np.exp(-1))*j:.0f}), max weight {w.max():.4f}', flush=True)
        if unique >= (1 - np.exp(-1))*j or it == max_iter:
            break
        # A new component at the highest-weight point
        i = int(np.argmax(w))
        diff = Z - Z[i]
        dist = np.einsum('ij,jk,ik->i', diff, prior_cov_inv, diff)
        near = np.argsort(dist)[:b]
        wn = w[near] + 1/n # Raftery and Bao: average the importance weights with equal weights
        wn /= wn.sum()
        cov = np.cov(Z[near].T, aweights=wn) + 1e-10*np.eye(d)
        comps.append((Z[i].copy(), cov))
        Znew = rng.multivariate_normal(Z[i], cov, size=b)
        new_outs = [sc.objdict({k: np.nan for k in list(targets) + list(reported)}) for _ in range(b)]
        new_ll = np.full(b, -np.inf)
        ok = np.all((Znew >= 0) & (Znew <= 1), axis=1) # Outside the prior: weight 0, no need to simulate
        if ok.any():
            o, l = evaluate(Znew[ok], tgts)
            for k, idx in enumerate(np.flatnonzero(ok)):
                new_outs[idx] = o[k]
            new_ll[ok] = l
        Z = np.vstack([Z, Znew])
        outs += new_outs
        ll = np.concatenate([ll, new_ll])
    return sc.objdict(Z=Z, theta=from_unit(Z), outs=outs, ll=ll, weights=w, comps=comps, log=log, n0=n0, b=b, j=j, seed=seed)
