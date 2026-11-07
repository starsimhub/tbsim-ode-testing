"""
A standalone copy of tbsim v4's TB ODE (tbsim.TB as declared at tbsim dcfb522, run by Starsim's compartmental solver), for
calibration: the same equations as one function of the state, with the equilibrium found by root-finding (milliseconds, against
about a second through Starsim). Checked against v4's ODE in tests/test_fastode.py and exp 19; rerun that test after bumping the
v4 pins.

Scope: homogeneous mixing, no treatment (diagnosis rate theta = 0), births replacing every death (mu, as tbsim.TB's own),
optionally the two latent states (k_latent > 0). Shares of a constant population.

One state is split for counting only: non-infectious TB is NON0 (not yet asymptomatic in this episode) and NON1 (has been), with
identical rates, so incidence can be counted as the first entry to asymptomatic TB per episode (CONVENTIONS.md 5): progression
from latent infection to asymptomatic TB, plus NON0 -> asymptomatic. NON0 + NON1 is v4's `non_infectious`.
"""
import numpy as np
import sciris as sc
from scipy.integrate import solve_ivp
from scipy.optimize import root

__all__ = ['default_pars', 'make_pars', 'state_names', 'rhs', 'flows', 'run', 'equilibrium', 'summary']

# tbsim.tb_pars() at tbsim dcfb522 (rates per year), the parameters used here
default_pars = dict(
    beta=9.0, trans_asymp=0.82, rr_reinfection_cleared=1.00, rr_reinfection_rec=0.21, rr_reinfection_treat=3.15, rr_reinfection_inf=0.21,
    inf_cle=1.90, inf_non=0.16, inf_asy=0.06, k_latent=0.0, inf_non_remote=None, inf_asy_remote=0.0,
    non_rec=0.18, non_asy=0.25, asy_non=1.66, asy_sym=0.88, sym_asy=0.54, sym_dead=0.34, mu=1/70,
)

state_names = ['susceptible', 'infection', 'infection_remote', 'cleared', 'recovered', 'non0', 'non1', 'asymptomatic', 'symptomatic']
S, L, M, C, R, N0, N1, A, Y = range(len(state_names))


def make_pars(**kwargs):
    """ The parameters (rates per year), tbsim's defaults updated by kwargs; as tbsim.TB, inf_non_remote=None means inf_non """
    unknown = set(kwargs) - set(default_pars)
    if unknown:
        raise KeyError(f'Not parameters of this model: {sorted(unknown)}')
    p = sc.mergedicts(default_pars, kwargs)
    if p['inf_non_remote'] is None:
        p['inf_non_remote'] = p['inf_non']
    return {k: float(v) for k, v in p.items()}


def flows(x, p):
    """ Each transition's flow (per year, as a share of the population) at state x """
    lam = p['beta']*(p['trans_asymp']*x[A] + x[Y]) # Homogeneous mixing: exposure is the share of the population in each state
    N = x[N0] + x[N1]
    f = sc.objdict(
        transmission     = lam*x[S],
        reinfect_cle     = p['rr_reinfection_cleared']*lam*x[C],
        reinfect_rec     = p['rr_reinfection_rec']*lam*x[R],
        clearance        = p['inf_cle']*x[L],
        latent_non       = p['inf_non']*x[L],
        latent_asy       = p['inf_asy']*x[L],
        remote           = p['k_latent']*x[L],
        reinfect_remote  = p['rr_reinfection_inf']*lam*x[M],
        clearance_remote = p['inf_cle']*x[M],
        remote_non       = p['inf_non_remote']*x[M],
        remote_asy       = p['inf_asy_remote']*x[M],
        recovery         = p['non_rec']*N,
        worsening        = p['non_asy']*N,
        worsening_first  = p['non_asy']*x[N0], # The part of worsening that is a first entry to asymptomatic TB in the episode
        improvement      = p['asy_non']*x[A],
        symptoms         = p['asy_sym']*x[A],
        remission        = p['sym_asy']*x[Y],
        tb_death         = p['sym_dead']*x[Y],
    )
    f.incidence = f.latent_asy + f.remote_asy + f.worsening_first # First entry to asymptomatic TB per episode
    f.into_asy = f.latent_asy + f.remote_asy + f.worsening # Every entry from latent infection or non-infectious TB (not remission)
    return f


def rhs(t, x, p):
    """ dx/dt (shares of a constant population; births replace background and TB deaths) """
    f = flows(x, p)
    mu = p['mu']
    dx = np.empty_like(x)
    dx[S] = mu + f.tb_death - f.transmission - mu*x[S]
    dx[L] = f.transmission + f.reinfect_cle + f.reinfect_rec + f.reinfect_remote - f.clearance - f.latent_non - f.latent_asy - f.remote - mu*x[L]
    dx[M] = f.remote - f.reinfect_remote - f.clearance_remote - f.remote_non - f.remote_asy - mu*x[M]
    dx[C] = f.clearance + f.clearance_remote - f.reinfect_cle - mu*x[C]
    dx[R] = f.recovery - f.reinfect_rec - mu*x[R]
    out_non = p['non_rec'] + p['non_asy'] + mu
    dx[N0] = f.latent_non + f.remote_non - out_non*x[N0]
    dx[N1] = f.improvement - out_non*x[N1]
    dx[A] = f.latent_asy + f.remote_asy + f.worsening + f.remission - f.improvement - f.symptoms - mu*x[A]
    dx[Y] = f.symptoms - f.remission - f.tb_death - mu*x[Y]
    return dx


def initial(init=None):
    """ A state vector from a dict of shares by v4's state names (non_infectious goes to NON0); the rest susceptible """
    x = np.zeros(len(state_names))
    init = sc.ifelse(init, dict(symptomatic=0.01))
    for k, v in init.items():
        if k == 'non_infectious':
            x[N0] += v
        elif k in state_names:
            x[state_names.index(k)] += v
        elif v:
            raise KeyError(f'No state {k} in this model')
    x[S] += 1 - x.sum()
    return x


def run(p, init=None, years=100, t_eval=None, rtol=1e-10, atol=1e-13):
    """ The trajectory from init over `years` (solve_ivp, LSODA): (t, X (times × states)) """
    sol = solve_ivp(rhs, (0, years), initial(init) if not isinstance(init, np.ndarray) else init, args=(p,), method='LSODA',
                    t_eval=t_eval, rtol=rtol, atol=atol)
    return sol.t, sol.y.T


def equilibrium(p, x0=None, burn=300, tol=1e-12):
    """
    The endemic equilibrium: a burn-in run from x0 (or a 1% symptomatic seed), then root-finding on dx/dt = 0 with the shares
    summing to 1. Returns the state vector, or None if there is no endemic equilibrium (the burn-in dies out) or it fails
    """
    x = x0 if x0 is not None else initial()
    if burn:
        _, X = run(p, x, years=burn, t_eval=[burn], rtol=1e-8, atol=1e-11)
        x = X[-1]
    if x[A] + x[Y] < 1e-9:
        return None
    def g(z):
        d = rhs(0, z, p)
        d[S] = z.sum() - 1 # Replace one equation (dx sum to 0 by construction) by the constraint
        return d
    sol = root(g, x, method='hybr', tol=tol)
    if not sol.success or np.any(sol.x < -1e-12) or sol.x[A] + sol.x[Y] < 1e-9:
        return None
    return np.clip(sol.x, 0, None)


def summary(x, p):
    """ The calibration outputs at state x: prevalence (asymptomatic + symptomatic), share asymptomatic, incidence and mortality per person per year, and their ratios """
    f = flows(x, p)
    prev = x[A] + x[Y]
    return sc.objdict(prevalence=prev, prop_asymptomatic=x[A]/prev, incidence=f.incidence, into_asy=f.into_asy, mortality=f.tb_death,
                      mort_to_inc=f.tb_death/f.incidence, prev_to_inc=prev/f.incidence, latent=x[L] + x[M], non_infectious=x[N0] + x[N1])
