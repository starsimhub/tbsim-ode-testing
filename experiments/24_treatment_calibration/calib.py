"""
Experiment 24: priors, model, likelihood and IMIS for calibrating the v4 ODE with treatment (exp 23's design, a fixed course), time-varying
sym_dead and over-detection of notifications, to tb_LAI_TPT's applicable targets. See README.md.

The model is tbode4.fastode: each draw starts at its untreated endemic equilibrium at the end of 1994 (everything is constant before
1995), then follows tb_LAI_TPT's treatment scale-up and sym_dead trajectory to 2040. Incidence is the first entry to asymptomatic TB
per episode plus relapses (CONVENTIONS.md 5). The progression targets come from tbode4.latent's exact chain.
"""
import numpy as np
import sciris as sc
from scipy.integrate import solve_ivp
from scipy.stats import binom, multivariate_normal
import tbode4
from tbode4 import fastode as fo

__all__ = ['priors', 'names', 'targets', 'to_pars', 'trajectories', 'simulate', 'loglik', 'imis']

# Priors: (low, high, log-uniform?). Exp 20's natural history (sym_dead as its 1995 value), plus tb_LAI_TPT's treatment, mortality-trend
# and over-detection parameters (run_history_matching.parameter_bounds), with its parameterisations that keep the trends monotone
priors = sc.objdict(
    beta             = (2.0, 100.0, True),
    trans_asymp      = (0.5, 1.0, False),
    pi               = (0.2, 1.0, False),  # Every reinfection risk: cleared, recovered, treated, and remote latent -> recent
    inf_cle          = (0.5, 1.5, True),
    inf_non          = (0.04, 0.45, False),
    inf_asy          = (0.01, 0.60, True),
    k_latent         = (0.1, 12.0, False),
    phi              = (0.0, 0.1, False),  # Remote progression ÷ recent
    non_rec          = (0.14, 0.23, False),
    non_asy          = (0.21, 0.28, False),
    asy_non          = (0.70, 2.03, False),
    asy_sym          = (0.56, 0.99, False),
    sym_asy          = (0.46, 0.85, False),
    sym_dead_1995    = (0.20, 0.50, False),  # TB death from symptomatic TB, up to 1995
    sym_dead_frac    = (0.0, 1.0, False),    # sym_dead_2024 = sym_dead_1995 − frac·(sym_dead_1995 − 0.05): a decrease
    tx_covg_2024     = (0.50, 0.95, False),  # Annual treatment coverage of symptomatic TB in 2024
    tx_covg_2010_frac= (0.0, 1.0, False),    # tx_covg_2010 = frac·tx_covg_2024
    efficacy_1995    = (0.30, 0.95, False),  # Chance a course cures (adherence 1), from 1995
    efficacy_frac    = (0.0, 1.0, False),    # efficacy_2024 = efficacy_1995 + frac·(0.95 − efficacy_1995): an increase
    p_relapse        = (0.02, 0.10, False),
    over_detection   = (0.0, 0.20, False),   # Reported notifications = (1 + over_detection) × treatment starts
)
names = list(priors.keys())
d = len(names)
mu = 1/70 # Background mortality, with births replacing every death (tbsim's default)
sym_dead_floor, efficacy_ceiling = 0.05, 0.95 # tb_LAI_TPT's SYM_DEAD_FLOOR, EFF_CEIL
t0 = 1994 + 364/365 # tb_LAI_TPT's first treatment knot, 1994-12-31: the model is constant before it
stop = 2041.0
tx_stages = 12 # The fixed 6-month course, as Erlang stages (CONVENTIONS.md 9; within 0.1% of prevalence, exp 23)

# The targets (README). Means and data sigmas are tb_LAI_TPT's (hm_common.calibration_sigmas, real_targets.csv)
targets = sc.objdict(
    prevalence        = dict(kind='binomial', x=77, n=24815),                                # 2020 survey, 15+ (point)
    prop_asymptomatic = dict(kind='binomial', x=123, n=188),                                 # 2020 survey
    prev_to_inc       = dict(kind='lognormal', mean=310.296/197.482, sigma=0.1356),          # Prevalence 2020 ÷ incidence 2020
    mort_to_inc       = dict(kind='lognormal', mean=21.073/187.31, sigma=0.1271),            # Mortality 2024 ÷ incidence 2024
    prev_to_notif     = dict(kind='lognormal', mean=310.296/116.666, sigma=0.1247),          # Prevalence 2020 ÷ notifications 2020
    inc_ratio         = dict(kind='lognormal', mean=187.31/237.197, sigma=0.1119),           # Incidence 2024 ÷ 2015
    mort_ratio        = dict(kind='lognormal', mean=21.073/32.03, sigma=0.1204),             # Mortality 2024 ÷ 2015
    prev_2040_ratio   = dict(kind='lognormal', mean=0.8064319874724902, sigma=0.1901),       # Prevalence 2040 ÷ 2020 (tb_LAI_TPT's anchor)
    ProgressionY1     = dict(kind='lognormal', mean=0.0523, sigma=0.094),
    ProgressionRatio  = dict(kind='lognormal', mean=0.663, sigma=0.071),
)
extra = ['incidence_2024', 'notifications_2020', 'mortality_2024', 'prevalence_2024', 'prevalence_1994', 'incidence_1994'] # Reported, per person (per year)


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
    """ One parameter vector -> the natural history (tbsim/fastode parameters, rates per year) and the treatment and trend values """
    t = dict(zip(names, np.asarray(theta, dtype=float)))
    nat = dict(beta=t['beta'], trans_asymp=t['trans_asymp'], rr_reinfection_cleared=t['pi'], rr_reinfection_rec=t['pi'],
               rr_reinfection_treat=t['pi'], rr_reinfection_inf=t['pi'], inf_cle=t['inf_cle'], inf_non=t['inf_non'], inf_asy=t['inf_asy'],
               k_latent=t['k_latent'], inf_non_remote=t['phi']*t['inf_non'], inf_asy_remote=t['phi']*t['inf_asy'], non_rec=t['non_rec'],
               non_asy=t['non_asy'], asy_non=t['asy_non'], asy_sym=t['asy_sym'], sym_asy=t['sym_asy'], sym_dead=t['sym_dead_1995'], mu=mu)
    trend = sc.objdict(
        sym_dead_1995 = t['sym_dead_1995'],
        sym_dead_2024 = t['sym_dead_1995'] - t['sym_dead_frac']*(t['sym_dead_1995'] - sym_dead_floor),
        tx_covg_2010  = t['tx_covg_2010_frac']*t['tx_covg_2024'],
        tx_covg_2024  = t['tx_covg_2024'],
        efficacy_1995 = t['efficacy_1995'],
        efficacy_2024 = t['efficacy_1995'] + t['efficacy_frac']*(efficacy_ceiling - t['efficacy_1995']),
        p_relapse     = t['p_relapse'],
        over_detection= t['over_detection'],
    )
    return nat, trend


def schedules(trend):
    """ Coverage, efficacy and sym_dead as functions of calendar time, as tb_LAI_TPT interpolates them (linear between knots, then constant) """
    coverage = lambda y: np.interp(y, [t0, 2010.0, 2024.0], [0.0, trend.tx_covg_2010, trend.tx_covg_2024], left=0.0)
    efficacy = lambda y: np.interp(y, [t0, 1995.0, 2024.0], [0.0, trend.efficacy_1995, trend.efficacy_2024], left=0.0)
    sym_dead = lambda y: np.interp(y, [1995.0, 2024.0], [trend.sym_dead_1995, trend.sym_dead_2024])
    return coverage, efficacy, sym_dead


def trajectories(theta, t_eval=None, rtol=1e-8, atol=1e-11):
    """
    One draw's run: the untreated equilibrium at t0, then the scale-up to `stop`. Returns the times (calendar years), the states, the
    cumulative incidence, treatment starts and TB deaths (per person), and fastode's parameters; None if there is no endemic equilibrium
    """
    nat, trend = to_pars(theta)
    p0 = fo.make_pars(**nat, tx_stages=tx_stages)
    x0 = fo.equilibrium(p0)
    if x0 is None:
        return None
    coverage, efficacy, sym_dead = schedules(trend)
    p = fo.make_pars(**{**nat, 'sym_dead': lambda t: sym_dead(t0 + t)}, tx_stages=tx_stages, tx_p_relapse=trend.p_relapse,
                     tx_rate=lambda t: float(-np.log1p(-coverage(t0 + t))), tx_success=lambda t: float(efficacy(t0 + t)))
    k = tx_stages
    fail = slice(fo.TX + k, fo.TX + 2*k)

    def rhs(t, z):
        x = z[:-3]
        dx = fo.rhs(t, x, p)
        sd = fo.value(p, 'sym_dead', t)
        inc = p['inf_asy']*x[fo.L] + p['inf_asy_remote']*x[fo.M] + p['non_asy']*x[fo.N0] + x[fo.RL]/p['dur_relapse']
        return np.concatenate([dx, [inc, fo.value(p, 'tx_rate', t)*x[fo.Y], sd*(x[fo.Y] + x[fail].sum())]])

    years = sc.ifelse(t_eval, np.array([t0, 1995.0, 2015.0, 2016.0, 2020.0, 2021.0, 2024.0, 2025.0, 2040.0, 2041.0]))
    sol = solve_ivp(rhs, (0, stop - t0), np.concatenate([x0, np.zeros(3)]), t_eval=np.asarray(years) - t0, method='LSODA', rtol=rtol, atol=atol)
    if not sol.success:
        return None
    Z = sol.y.T
    return sc.objdict(year=np.asarray(years), X=Z[:, :-3], cum_inc=Z[:, -3], cum_notif=Z[:, -2], cum_deaths=Z[:, -1], p=p, nat=nat, trend=trend)


def simulate(theta):
    """ The target outputs at one parameter vector (NaN for every model output if there is no endemic equilibrium) """
    out = sc.objdict({k: np.nan for k in list(targets) + extra})
    r = trajectories(theta)
    nat, trend = to_pars(theta)
    if r is not None:
        at = {float(y): i for i, y in enumerate(r.year)}
        prev = lambda y: r.X[at[y], fo.A] + r.X[at[y], fo.Y]
        annual = lambda cum, y: cum[at[y + 1]] - cum[at[y]]
        inc = {y: annual(r.cum_inc, y) for y in [2015.0, 2020.0, 2024.0]}
        deaths = {y: annual(r.cum_deaths, y) for y in [2015.0, 2024.0]}
        notif_2020 = (1 + trend.over_detection)*annual(r.cum_notif, 2020.0)
        x20 = r.X[at[2020.0]]
        out.update(
            prevalence        = prev(2020.0),
            prop_asymptomatic = x20[fo.A]/(x20[fo.A] + x20[fo.Y]),
            prev_to_inc       = prev(2020.0)/inc[2020.0],
            mort_to_inc       = deaths[2024.0]/inc[2024.0],
            prev_to_notif     = prev(2020.0)/notif_2020 if notif_2020 > 0 else np.nan,
            inc_ratio         = inc[2024.0]/inc[2015.0],
            mort_ratio        = deaths[2024.0]/deaths[2015.0],
            prev_2040_ratio   = prev(2040.0)/prev(2020.0),
            incidence_2024    = inc[2024.0],
            notifications_2020= notif_2020,
            mortality_2024    = deaths[2024.0],
            prevalence_2024   = prev(2024.0),
            prevalence_1994   = prev(t0),
            incidence_1994    = fo.flows(r.X[0], fo.make_pars(**nat, tx_stages=tx_stages)).incidence,
        )
    chain = dict(inf_cle=nat['inf_cle'], inf_non=nat['inf_non'], inf_asy=nat['inf_asy'], k_latent=nat['k_latent'], inf_non_remote=nat['inf_non_remote'],
                 inf_asy_remote=nat['inf_asy_remote'], non_rec=nat['non_rec'], non_asy=nat['non_asy'])
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


def imis(tgts=None, n0=10_000, b=1_000, j=3_000, max_iter=300, seed=1, verbose=True):
    """
    Incremental mixture importance sampling (Raftery and Bao 2010), on the unit cube where the prior is uniform; as exp 20's.

    Starts with n0 prior draws; each iteration adds b draws from a multivariate normal centred on the current highest-weight point,
    with the weighted covariance of its b nearest points (Mahalanobis distance under the prior's covariance); stops when the expected
    number of unique points in a resample of j reaches (1 − 1/e)·j. Importance weights are likelihood × prior ÷ the mixture of the prior
    and every component so far.
    """
    rng = np.random.default_rng(seed)
    Z = rng.random((n0, d))
    outs, ll = evaluate(Z, tgts)
    comps = []
    log = []
    prior_cov_inv = np.eye(d)*12.0 # Covariance of a uniform on [0, 1] is 1/12
    T = sc.timer()
    for it in range(max_iter + 1):
        n = len(Z)
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
            print(f'IMIS {it:3d}: n = {n:7,}, ESS = {ess:8.1f}, expected unique in {j} = {unique:7.1f} (target {(1 - np.exp(-1))*j:.0f}), max weight {w.max():.4f}, {T.toc(output=True):.0f} s', flush=True)
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
        Znew = rng.multivariate_normal(Z[i], cov, size=b)
        new_outs = [sc.objdict({k: np.nan for k in list(targets) + extra}) for _ in range(b)]
        new_ll = np.full(b, -np.inf)
        ok = np.all((Znew >= 0) & (Znew <= 1), axis=1)
        if ok.any():
            o, l = evaluate(Znew[ok], tgts)
            for kk, idx in enumerate(np.flatnonzero(ok)):
                new_outs[idx] = o[kk]
            new_ll[ok] = l
        Z = np.vstack([Z, Znew])
        outs += new_outs
        ll = np.concatenate([ll, new_ll])
    return sc.objdict(Z=Z, theta=from_unit(Z), outs=outs, ll=ll, weights=w, comps=comps, log=log, n0=n0, b=b, j=j, seed=seed)
