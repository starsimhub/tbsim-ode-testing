"""
A standalone copy of tbsim v4's TB ODE (tbsim.TB as declared at tbsim dcfb522 and unchanged at 619c296, run by Starsim's compartmental solver), for
calibration: the same equations as one function of the state, with the equilibrium found by root-finding (milliseconds, against
about a second through Starsim). Checked against v4's ODE in tests/test_fastode.py and exp 19; rerun that test after bumping the
v4 pins.

Scope: homogeneous mixing (the network option in force() is a reference only, not for calibration), births replacing every death
(mu, as tbsim.TB's own), optionally the two latent states (k_latent > 0). Shares of a constant population.

Treatment (optional, tx_rate > 0) is exp 23's design (CONVENTIONS.md 9), as declared in tbode4.TBTreat: people with symptomatic TB start
treatment at tx_rate, entering a course that will cure them (probability tx_success) or fail. A course lasts tx_dur, as tx_stages Erlang
stages (1: a rate; many: nearly fixed; TBTreat's exact fixed course is a conveyor of whole steps, which this copy approximates by stages).
Failing courses face TB death at sym_dead and return to symptomatic TB; a share tx_p_relapse of cures relapse to symptomatic TB at rate
1/dur_relapse, the rest go to `treated`. tx_rate, tx_success and sym_dead can be functions of time (years from the start of the run).

One state is split for counting only: non-infectious TB is NON0 (not yet asymptomatic in this episode) and NON1 (has been), with
identical rates, so incidence can be counted as the first entry to asymptomatic TB per episode (CONVENTIONS.md 5): progression
from latent infection to asymptomatic TB, plus NON0 -> asymptomatic, plus relapse after treatment. NON0 + NON1 is v4's `non_infectious`.

The state vector is the 9 states in `state_names`, then `treated`, `relapsing`, and the tx_stages stages of each course (see names()).
"""
import numpy as np
import sciris as sc
from scipy.integrate import solve_ivp
from scipy.optimize import root

__all__ = ['default_pars', 'make_pars', 'state_names', 'names', 'force', 'rhs', 'flows', 'run', 'equilibrium', 'summary']

# tbsim.tb_pars() at tbsim dcfb522 and 619c296 (rates per year), the parameters used here
default_pars = dict(
    beta=9.0, trans_asymp=0.82, rr_reinfection_cleared=1.00, rr_reinfection_rec=0.21, rr_reinfection_treat=3.15, rr_reinfection_inf=0.21,
    inf_cle=1.90, inf_non=0.16, inf_asy=0.06, k_latent=0.0, inf_non_remote=None, inf_asy_remote=0.0,
    non_rec=0.18, non_asy=0.25, asy_non=1.66, asy_sym=0.88, sym_asy=0.54, sym_dead=0.34, mu=1/70,
    net_contacts=0.0, net_step=0.0, net_source_exp=0.0, # A random network with fresh contacts each step (off: 0); see force()
    tx_rate=0.0, tx_success=1.0, tx_dur=0.5, tx_stages=1, tx_p_relapse=0.0, dur_relapse=1.5, # Treatment (off: tx_rate = 0); see the docstring
)

state_names = ['susceptible', 'infection', 'infection_remote', 'cleared', 'recovered', 'non0', 'non1', 'asymptomatic', 'symptomatic']
S, L, M, C, R, N0, N1, A, Y = range(len(state_names))
TR, RL, TX = 9, 10, 11 # Treated, relapsing, and the first stage of the courses (successes' stages, then failures')


def names(p):
    """ Every state's name, in the order of the state vector """
    k = int(p['tx_stages'])
    return state_names + ['treated', 'relapsing'] + [f'treatment_success_{i}' for i in range(k)] + [f'treatment_fail_{i}' for i in range(k)]


def n_states(p):
    return len(state_names) + 2 + 2*int(p['tx_stages'])


def value(p, key, t):
    """ A parameter's value at time t (years from the start), for those that can be functions of time """
    v = p[key]
    return v(t) if callable(v) else v


def make_pars(**kwargs):
    """ The parameters (rates per year), tbsim's defaults updated by kwargs; as tbsim.TB, inf_non_remote=None means inf_non """
    unknown = set(kwargs) - set(default_pars)
    if unknown:
        raise KeyError(f'Not parameters of this model: {sorted(unknown)}')
    p = sc.mergedicts(default_pars, kwargs)
    if p['inf_non_remote'] is None:
        p['inf_non_remote'] = p['inf_non']
    return {k: (v if callable(v) else int(v) if k == 'tx_stages' else float(v)) for k, v in p.items()}


def force(p, x, rr=1.0):
    """
    The force of infection on someone at relative risk rr (per year).

    Homogeneous mixing (net_contacts = 0): rr·beta·(trans_asymp·A + Y), exposure being the share of the population in each state.

    A random network with fresh contacts each step (net_contacts = K mean contacts, net_step = D years; beta is then per contact): each
    contact transmits at most once in a step. As v4's agents do it, the target's risk saturates at its own rate and the source's
    infectiousness scales the chance linearly: a symptomatic contact transmits with chance 1 − e^(−rr·beta·D), an asymptomatic one with
    trans_asymp times that, so the force is (K/D)·(1 − e^(−rr·beta·D))·(trans_asymp·A + Y). This is what agents on ss.RandomNet(dur=0) do
    (exp 22's diagnose.py, starsim 80c81f53), and what proposal 03 asks Starsim's ODE of a network to do (a reference version of that fix).

    With net_source_exp = 1, the source's infectiousness is inside the exponent too, as agents do it from starsim 6dcad9fd (exp 25):
    (K/D)·[(1 − e^(−rr·trans_asymp·beta·D))·A + (1 − e^(−rr·beta·D))·Y].

    The network option (net_contacts, net_step, net_source_exp) is a reference only, kept so exps 22 and 25 can be rerun: to calibrate on a
    network, use Starsim's ODE of the network (CONVENTIONS.md 1), which net_source_exp = 1 matches to within 0.4% (exp 25).
    """
    if not p['net_contacts']:
        return rr*p['beta']*(p['trans_asymp']*x[A] + x[Y])
    K, D, b = p['net_contacts'], p['net_step'], p['beta']
    if p['net_source_exp']:
        return (K/D)*(-np.expm1(-rr*p['trans_asymp']*b*D)*x[A] - np.expm1(-rr*b*D)*x[Y])
    return (K/D)*(-np.expm1(-rr*b*D))*(p['trans_asymp']*x[A] + x[Y])


def flows(x, p, t=0.0):
    """ Each transition's flow (per year, as a share of the population) at state x and time t """
    N = x[N0] + x[N1]
    k = int(p['tx_stages'])
    succ, fail = x[TX:TX + k], x[TX + k:TX + 2*k] # The stages of courses that will cure and that will fail
    exit_rate = k/p['tx_dur'] # Each stage's rate of moving on
    sym_dead = value(p, 'sym_dead', t)
    f = sc.objdict(
        transmission     = force(p, x)*x[S],
        reinfect_cle     = force(p, x, p['rr_reinfection_cleared'])*x[C],
        reinfect_rec     = force(p, x, p['rr_reinfection_rec'])*x[R],
        clearance        = p['inf_cle']*x[L],
        latent_non       = p['inf_non']*x[L],
        latent_asy       = p['inf_asy']*x[L],
        remote           = p['k_latent']*x[L],
        reinfect_remote  = force(p, x, p['rr_reinfection_inf'])*x[M],
        clearance_remote = p['inf_cle']*x[M],
        remote_non       = p['inf_non_remote']*x[M],
        remote_asy       = p['inf_asy_remote']*x[M],
        recovery         = p['non_rec']*N,
        worsening        = p['non_asy']*N,
        worsening_first  = p['non_asy']*x[N0], # The part of worsening that is a first entry to asymptomatic TB in the episode
        improvement      = p['asy_non']*x[A],
        symptoms         = p['asy_sym']*x[A],
        remission        = p['sym_asy']*x[Y],
        tb_death         = sym_dead*x[Y],
        reinfect_trd     = force(p, x, p['rr_reinfection_treat'])*x[TR],
        tx_start         = value(p, 'tx_rate', t)*x[Y], # Notifications
        tx_cure          = exit_rate*succ[-1],
        tx_fail          = exit_rate*fail[-1],
        tb_death_tx      = sym_dead*fail.sum(), # TB deaths of people on courses that will fail
        relapse          = x[RL]/p['dur_relapse'],
    )
    f.tb_deaths = f.tb_death + f.tb_death_tx # TB mortality
    f.incidence = f.latent_asy + f.remote_asy + f.worsening_first + f.relapse # First entry to asymptomatic TB per episode, and relapses
    f.into_asy = f.latent_asy + f.remote_asy + f.worsening # Every entry from latent infection or non-infectious TB (not remission)
    return f


def rhs(t, x, p):
    """ dx/dt (shares of a constant population; births replace background and TB deaths) """
    f = flows(x, p, t)
    mu = p['mu']
    dx = np.empty_like(x)
    dx[S] = mu + f.tb_deaths - f.transmission - mu*x[S]
    dx[L] = f.transmission + f.reinfect_cle + f.reinfect_rec + f.reinfect_trd + f.reinfect_remote - f.clearance - f.latent_non - f.latent_asy - f.remote - mu*x[L]
    dx[M] = f.remote - f.reinfect_remote - f.clearance_remote - f.remote_non - f.remote_asy - mu*x[M]
    dx[C] = f.clearance + f.clearance_remote - f.reinfect_cle - mu*x[C]
    dx[R] = f.recovery - f.reinfect_rec - mu*x[R]
    out_non = p['non_rec'] + p['non_asy'] + mu
    dx[N0] = f.latent_non + f.remote_non - out_non*x[N0]
    dx[N1] = f.improvement - out_non*x[N1]
    dx[A] = f.latent_asy + f.remote_asy + f.worsening + f.remission - f.improvement - f.symptoms - mu*x[A]
    dx[Y] = f.symptoms - f.remission - f.tb_death - f.tx_start + f.tx_fail + f.relapse - mu*x[Y]

    # Treatment: courses as Erlang stages; cures go to treated or relapsing, failures back to symptomatic TB
    k = int(p['tx_stages'])
    s = value(p, 'tx_success', t)
    rate = k/p['tx_dur']
    for first, inflow, death in [(TX, s*f.tx_start, 0.0), (TX + k, (1 - s)*f.tx_start, value(p, 'sym_dead', t))]:
        stages = x[first:first + k]
        into = np.concatenate([[inflow], rate*stages[:-1]])
        dx[first:first + k] = into - (rate + death + mu)*stages
    dx[TR] = (1 - p['tx_p_relapse'])*f.tx_cure - f.reinfect_trd - mu*x[TR]
    dx[RL] = p['tx_p_relapse']*f.tx_cure - f.relapse - mu*x[RL]
    return dx


def initial(init=None, p=None):
    """ A state vector from a dict of shares by v4's state names (non_infectious goes to NON0); the rest susceptible """
    p = sc.ifelse(p, default_pars)
    allnames = names(p)
    x = np.zeros(n_states(p))
    init = sc.ifelse(init, dict(symptomatic=0.01))
    for k, v in init.items():
        if k == 'non_infectious':
            x[N0] += v
        elif k in allnames:
            x[allnames.index(k)] += v
        elif v:
            raise KeyError(f'No state {k} in this model')
    x[S] += 1 - x.sum()
    return x


def run(p, init=None, years=100, t_eval=None, rtol=1e-10, atol=1e-13):
    """ The trajectory from init over `years` (solve_ivp, LSODA): (t, X (times × states)) """
    sol = solve_ivp(rhs, (0, years), initial(init, p) if not isinstance(init, np.ndarray) else init, args=(p,), method='LSODA',
                    t_eval=t_eval, rtol=rtol, atol=atol)
    return sol.t, sol.y.T


def equilibrium(p, x0=None, burn=300, tol=1e-12):
    """
    The endemic equilibrium: a burn-in run from x0 (or a 1% symptomatic seed), then root-finding on dx/dt = 0 with the shares
    summing to 1. Returns the state vector, or None if there is no endemic equilibrium (the burn-in dies out) or it fails
    """
    x = x0 if x0 is not None else initial(p=p)
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
    if not sol.success: # The root-finder stalled: burn in longer, and retry from there
        _, X = run(p, x, years=2000, t_eval=[2000], rtol=1e-10, atol=1e-13)
        x = X[-1]
        sol = root(g, x, method='hybr', tol=tol)
        if not sol.success and np.max(np.abs(rhs(0, x, p))) < 1e-12: # Accept the burn-in's end if it is already at rest
            sol = sc.objdict(success=True, x=x)
    if not sol.success or np.any(sol.x < -1e-12) or sol.x[A] + sol.x[Y] < 1e-9:
        return None
    return np.clip(sol.x, 0, None)


def summary(x, p, t=0.0):
    """
    The calibration outputs at state x: prevalence (asymptomatic + symptomatic), share asymptomatic, incidence (with relapses) and TB
    mortality (with deaths on treatment) per person per year, their ratios, and with treatment, notifications and the share on treatment
    """
    f = flows(x, p, t)
    prev = x[A] + x[Y]
    on_tx = x[TX:].sum()
    return sc.objdict(prevalence=prev, prop_asymptomatic=x[A]/prev, incidence=f.incidence, into_asy=f.into_asy, mortality=f.tb_deaths,
                      mort_to_inc=f.tb_deaths/f.incidence, prev_to_inc=prev/f.incidence, latent=x[L] + x[M], non_infectious=x[N0] + x[N1],
                      notifications=f.tx_start, notif_to_inc=f.tx_start/f.incidence, on_treatment=on_tx, relapse=f.relapse, mortality_tx=f.tb_death_tx)
