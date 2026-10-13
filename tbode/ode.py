"""
ODE version of TBsim, extended feature by feature to match the agent-based tbsim.TB.

Starts as a generalisation of tbsim.compartmental.TB_ODE: same compartments, equations and default
parameters, with options for each ABM feature as it is added. With all options at their defaults it
reproduces TB_ODE exactly.
"""

import numpy as np
import sciris as sc
from scipy.integrate import odeint, solve_ivp
from scipy.optimize import brentq, fsolve
from scipy.stats import poisson
from tbsim.compartmental.lshtm_ode import default_pars as lshtm_pars

__all__ = ['states', 'default_pars', 'latent_classes', 'TBODE', 'calibrate_beta', 'ode_cohort_progression']

# ODE compartments with constant progression (k_asy = k_non = 0). CLEARED, RECOVERED and TREATED together correspond
# to the ABM's single CLEARED state. With k > 0, INFECTION is split into latent classes (see TBODE.states)
states = ['SUSCEPTIBLE', 'INFECTION', 'CLEARED', 'RECOVERED', 'NON_INFECTIOUS',
          'ASYMPTOMATIC', 'SYMPTOMATIC', 'TREATMENT', 'TREATED']

# Auxiliary (bookkeeping-only) compartment appended after the real ones: NON_INFECTIOUS agents who have not yet
# been ASYMPTOMATIC in this episode. It is a subset of NON_INFECTIOUS, does not count towards N and does not
# feed back on the dynamics; it lets incidence be counted as in tb_LAI_TPT's IncidenceTracker (first ASYMPTOMATIC
# entry per episode: INFECTION -> ASYMPTOMATIC, plus NON_INFECTIOUS -> ASYMPTOMATIC only if not yet ASYMPTOMATIC)
aux_states = ['NON0']

default_pars = sc.mergedicts(lshtm_pars, dict(
    replace_tb_deaths = True,  # True: TB deaths re-enter SUSCEPTIBLE so N is constant (as TB_ODE). False: births = mu*N only, as in the ABM with ss.Births/ss.Deaths
    init_sym = 0.01,           # Initial fraction in SYMPTOMATIC (rest SUSCEPTIBLE); TB_ODE uses 1e3/1e5
    k_asy = 0.0,               # Decline of inf_asy with time since infection tau: inf_asy*exp(-k_asy*tau), as in tbsim.TB. Must equal k_non
    k_non = 0.0,               # Decline of inf_non with time since infection, as k_asy
    latent_tol = 1e-10,        # With k > 0: the last latent class also holds all higher classes, whose Poisson weight is below this
))


def latent_classes(inf_non, inf_asy, k, tol=1e-10):
    """
    Poisson-mixture representation of INFECTION-exit hazards inf_non*exp(-k*tau) and inf_asy*exp(-k*tau).

    The time to progression has survival exp(-(r0/k)*(1 - exp(-k*tau))), with r0 = inf_non + inf_asy, which equals
    sum_m Poisson(m; r0/k) * exp(-m*k*tau): a mixture of exponentials. So an agent assigned at infection to class m
    with probability Poisson(m; r0/k), and progressing at constant rate m*k, has exactly the same progression time.
    Class 0 never progresses. The destination is NON_INFECTIOUS or ASYMPTOMATIC in the ratio inf_non : inf_asy.
    Clearance (and death) act identically on every class, so the mixture stays exact with competing risks.

    Returns (weights, rate_non, rate_asy): per-class assignment probabilities and progression rates, truncated where
    the remaining Poisson weight is below ``tol`` (added to the last class).
    """
    r0 = inf_non + inf_asy
    lam = r0 / k
    M = int(poisson.isf(tol, lam)) + 1 if lam > 0 else 0
    m = np.arange(M + 1)
    weights = poisson.pmf(m, lam)
    weights[-1] += 1 - weights.sum()
    frac_non = inf_non / r0 if r0 > 0 else 0
    return weights, m * k * frac_non, m * k * (1 - frac_non)


class TBODE(sc.prettyobj):
    """
    Compartmental TB model integrated with scipy.odeint.

    With k_asy = k_non = k > 0 (progression declining with time since infection, as in tbsim.TB), INFECTION is
    split into latent classes INFECTION_0, INFECTION_1, ...: see ``latent_classes``. This is exact, not an
    approximation, and ``outputs`` sums the classes back into INFECTION.

    Args:
        kwargs: overrides for any entry in ``default_pars``

    **Example**::

        ode = tbode.TBODE(replace_tb_deaths=False)
        res = ode.run(start=1900, stop=2100, dt=1/12)
        eq = ode.equilibrium()
    """
    def __init__(self, **kwargs):
        self.pars = sc.objdict(sc.dcp(default_pars))
        unknown = set(kwargs) - set(self.pars)
        if unknown:
            raise KeyError(f'Unknown parameter(s): {unknown}')
        self.pars.update(kwargs)
        self.results = None

        p = self.pars
        if p.k_asy != p.k_non:
            raise NotImplementedError('TBODE requires k_asy = k_non')
        if p.k_asy > 0:
            self.latent_weights, self.rate_non, self.rate_asy = latent_classes(p.inf_non, p.inf_asy, p.k_asy, p.latent_tol)
            self.latent = [f'INFECTION_{m}' for m in range(len(self.latent_weights))]
        else:
            self.latent_weights, self.rate_non, self.rate_asy = np.ones(1), np.array([p.inf_non]), np.array([p.inf_asy])
            self.latent = ['INFECTION']
        self.states = ['SUSCEPTIBLE'] + self.latent + states[2:]  # Real compartments; aux_states follow
        self.n_states = len(self.states)
        self.n_latent = len(self.latent)
        return

    def idx(self, name):
        """Index of a compartment in the state vector."""
        return (self.states + aux_states).index(name)

    def derivs(self, y, t):
        """Right-hand side of the ODE system (real compartments followed by auxiliary ones)."""
        p = self.pars
        nL = self.n_latent
        S = y[0]
        L = y[1:1 + nL]
        CLE, REC, NON, ASY, SYM, TX, TXD, NON0 = y[1 + nL:]
        N = y[:self.n_states].sum()

        foi = (p.beta / N) * (p.trans_asymp * ASY + SYM)
        tb_deaths = p.sym_dead * SYM
        births = p.mu * N + (tb_deaths if p.replace_tb_deaths else 0)
        to_non = (self.rate_non * L).sum()
        to_asy = (self.rate_asy * L).sum()

        dS    = births - foi * S - p.mu * S
        dL    = self.latent_weights * foi * (S + p.rr_reinfection_cleared * CLE + p.rr_reinfection_rec * REC + p.rr_reinfection_treat * TXD) \
                - (p.inf_cle + self.rate_non + self.rate_asy + p.mu) * L
        dCLE  = p.inf_cle * L.sum() - foi * p.rr_reinfection_cleared * CLE - p.mu * CLE
        dREC  = p.non_rec * NON - foi * p.rr_reinfection_rec * REC - p.mu * REC
        dNON  = to_non + p.asy_non * ASY - (p.non_rec + p.non_asy + p.mu) * NON
        dASY  = to_asy + p.non_asy * NON + p.sym_asy * SYM - (p.asy_non + p.asy_sym + p.mu) * ASY
        dSYM  = p.asy_sym * ASY - p.sym_asy * SYM - p.theta * SYM + p.phi * TX - tb_deaths - p.mu * SYM
        dTX   = p.theta * SYM - p.phi * TX - p.delta * TX - p.mu * TX
        dTXD  = p.delta * TX - foi * p.rr_reinfection_treat * TXD - p.mu * TXD
        dNON0 = to_non - (p.non_rec + p.non_asy + p.mu) * NON0
        return np.concatenate([[dS], dL, [dCLE, dREC, dNON, dASY, dSYM, dTX, dTXD, dNON0]])

    def init_conditions(self):
        """Initial compartment sizes."""
        p = self.pars
        y0 = np.zeros(self.n_states + len(aux_states))
        y0[self.idx('SUSCEPTIBLE')] = p.N * (1 - p.init_sym)
        y0[self.idx('SYMPTOMATIC')] = p.N * p.init_sym
        return y0

    def outputs(self, y, time=None):
        """Derived outputs from a state array (n_times x n_states+aux, or a single state vector)."""
        p = self.pars
        y = np.atleast_2d(y)
        res = sc.objdict()
        if time is not None:
            res.time = time
        for i, s in enumerate(self.states + aux_states):
            res[s] = y[:, i]
        L = y[:, 1:1 + self.n_latent]
        res.INFECTION = L.sum(axis=1)
        res.N = y[:, :self.n_states].sum(axis=1)
        res.prev_active   = (res.ASYMPTOMATIC + res.SYMPTOMATIC) / res.N                       # point prevalence, all ages
        res.frac_asymp    = res.ASYMPTOMATIC / (res.ASYMPTOMATIC + res.SYMPTOMATIC)             # of prevalent ASY+SYM
        res.incidence_py  = (L @ self.rate_asy + p.non_asy * res.NON0) / res.N                  # first ASY entry per episode, per person-year
        res.tb_deaths_py  = p.sym_dead * res.SYMPTOMATIC / res.N                                # TB deaths per person-year
        res.tb_deaths_ppy = res.tb_deaths_py  # alias used by compare.py
        return res

    def run(self, start=1500, stop=2020, dt=1.0, y0=None):
        """
        Integrate from ``start`` to ``stop`` (years), reporting every ``dt`` years, from ``y0`` (default:
        ``init_conditions()``). Returns an objdict of arrays (see ``outputs``).
        """
        times = np.arange(start, stop + dt / 2, dt)
        y0 = self.init_conditions() if y0 is None else y0
        out = odeint(lambda y, t: self.derivs(y, t), y0, times, rtol=1e-8, atol=1e-6)
        self.results = self.outputs(out, time=times)
        return self.results

    def _norm_derivs(self, x):
        """Derivatives of the population *proportions* (N=1), removing the growth/decline of N."""
        f = self.derivs(x, 0)
        return f - x * f[:self.n_states].sum()

    def equilibrium(self, burn=2000):
        """
        Endemic equilibrium as population proportions, returned as a state vector (sums to 1 over the real
        compartments). Integrates the proportion system for ``burn`` years from a 1% SYMPTOMATIC seed
        (regardless of ``init_sym``), then polishes with a root solve.
        """
        x0 = np.zeros(self.n_states + len(aux_states))
        x0[self.idx('SUSCEPTIBLE')] = 0.99
        x0[self.idx('SYMPTOMATIC')] = 0.01
        x = odeint(lambda x, t: self._norm_derivs(x), x0, [0, burn], rtol=1e-10, atol=1e-12, mxstep=100_000)[-1]
        if x[self.idx('ASYMPTOMATIC')] + x[self.idx('SYMPTOMATIC')] > 1e-12:
            x = fsolve(self._norm_derivs, x, xtol=1e-13)
        return x

    def latent_age_rates(self, x):
        """
        At equilibrium ``x``, the rate of the exponential distribution of time since infection within each latent
        class: inflow / stock (constant inflow and a constant per-capita exit rate, so the age within a class is
        exponential). Used to give ABM agents started at the equilibrium their time since infection.
        """
        p = self.pars
        nL = self.n_latent
        CLE, REC = x[self.idx('CLEARED')], x[self.idx('RECOVERED')]
        ASY, SYM, TXD = x[self.idx('ASYMPTOMATIC')], x[self.idx('SYMPTOMATIC')], x[self.idx('TREATED')]
        N = x[:self.n_states].sum()
        foi = (p.beta / N) * (p.trans_asymp * ASY + SYM)
        inflow = self.latent_weights * foi * (x[0] + p.rr_reinfection_cleared * CLE + p.rr_reinfection_rec * REC + p.rr_reinfection_treat * TXD)
        return inflow / x[1:1 + nL]


def calibrate_beta(target_prev, bounds=(1, 50), **pars):
    """
    ODE beta whose equilibrium active (ASY+SYM) prevalence, all ages, equals ``target_prev``.

    Args:
        target_prev (float): target prevalence as a fraction (e.g. 500e-5)
        bounds (tuple): bracketing interval for beta
        pars: other TBODE parameters, held fixed
    """
    def gap(beta):
        ode = TBODE(beta=beta, **pars)
        return ode.outputs(ode.equilibrium()).prev_active[0] - target_prev
    return brentq(gap, *bounds, xtol=1e-6)


def ode_cohort_progression(pars, years=(1, 10), k_asy=0.0, k_non=0.0):
    """
    Fraction of a synchronised, no-reinfection latent cohort that has ever reached ASYMPTOMATIC by each of ``years``.

    Mirrors tb_LAI_TPT's cohort_progression (curve B): everyone starts in INFECTION at t=0, no transmission, no
    background or TB death. ``k_asy``/``k_non`` give the exp(-k*tau) decline in the INFECTION-exit hazards.
    Returns {'B': array, 'ProgressionY1': B(1), 'ProgressionRatio': B(1)/B(10)} when years == (1, 10).
    """
    p = sc.objdict(pars)

    def f(t, y):
        I, N, B = y
        r_non, r_asy = p.inf_non * np.exp(-k_non * t), p.inf_asy * np.exp(-k_asy * t)
        return [-(p.inf_cle + r_non + r_asy) * I, r_non * I - (p.non_rec + p.non_asy) * N, r_asy * I + p.non_asy * N]

    sol = solve_ivp(f, [0, max(years)], [1, 0, 0], t_eval=list(years), rtol=1e-10, atol=1e-12)
    B = sol.y[2]
    out = sc.objdict(years=np.array(years), B=B)
    if tuple(years) == (1, 10):
        out.ProgressionY1 = B[0]
        out.ProgressionRatio = B[0] / B[1]
    return out
