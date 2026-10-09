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
from tbsim.compartmental.lshtm_ode import default_pars as lshtm_pars

__all__ = ['states', 'default_pars', 'TBODE', 'calibrate_beta', 'ode_cohort_progression']

# ODE compartments. CLEARED, RECOVERED and TREATED together correspond to the ABM's single CLEARED state
states = ['SUSCEPTIBLE', 'INFECTION', 'CLEARED', 'RECOVERED', 'NON_INFECTIOUS',
          'ASYMPTOMATIC', 'SYMPTOMATIC', 'TREATMENT', 'TREATED']
n_states = len(states)

# Auxiliary (bookkeeping-only) compartment appended after the real ones: NON_INFECTIOUS agents who have not yet
# been ASYMPTOMATIC in this episode. It is a subset of NON_INFECTIOUS, does not count towards N and does not
# feed back on the dynamics; it lets incidence be counted as in tb_LAI_TPT's IncidenceTracker (first ASYMPTOMATIC
# entry per episode: INFECTION -> ASYMPTOMATIC, plus NON_INFECTIOUS -> ASYMPTOMATIC only if not yet ASYMPTOMATIC)
aux_states = ['NON0']

default_pars = sc.mergedicts(lshtm_pars, dict(
    replace_tb_deaths = True,  # True: TB deaths re-enter SUSCEPTIBLE so N is constant (as TB_ODE). False: births = mu*N only, as in the ABM with ss.Births/ss.Deaths
    init_sym = 0.01,           # Initial fraction in SYMPTOMATIC (rest SUSCEPTIBLE); TB_ODE uses 1e3/1e5
))


class TBODE(sc.prettyobj):
    """
    Compartmental TB model integrated with scipy.odeint.

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
        return

    def derivs(self, y, t):
        """Right-hand side of the ODE system (real compartments followed by auxiliary ones)."""
        p = self.pars
        S, INF, CLE, REC, NON, ASY, SYM, TX, TXD, NON0 = y
        N = y[:n_states].sum()

        foi = (p.beta / N) * (p.trans_asymp * ASY + SYM)
        tb_deaths = p.sym_dead * SYM
        births = p.mu * N + (tb_deaths if p.replace_tb_deaths else 0)

        dS    = births - foi * S - p.mu * S
        dINF  = foi * (S + p.rr_reinfection_cleared * CLE + p.rr_reinfection_rec * REC + p.rr_reinfection_treat * TXD) \
                - (p.inf_cle + p.inf_non + p.inf_asy + p.mu) * INF
        dCLE  = p.inf_cle * INF - foi * p.rr_reinfection_cleared * CLE - p.mu * CLE
        dREC  = p.non_rec * NON - foi * p.rr_reinfection_rec * REC - p.mu * REC
        dNON  = p.inf_non * INF + p.asy_non * ASY - (p.non_rec + p.non_asy + p.mu) * NON
        dASY  = p.inf_asy * INF + p.non_asy * NON + p.sym_asy * SYM - (p.asy_non + p.asy_sym + p.mu) * ASY
        dSYM  = p.asy_sym * ASY - p.sym_asy * SYM - p.theta * SYM + p.phi * TX - tb_deaths - p.mu * SYM
        dTX   = p.theta * SYM - p.phi * TX - p.delta * TX - p.mu * TX
        dTXD  = p.delta * TX - foi * p.rr_reinfection_treat * TXD - p.mu * TXD
        dNON0 = p.inf_non * INF - (p.non_rec + p.non_asy + p.mu) * NON0
        return np.array([dS, dINF, dCLE, dREC, dNON, dASY, dSYM, dTX, dTXD, dNON0])

    def init_conditions(self):
        """Initial compartment sizes."""
        p = self.pars
        y0 = np.zeros(n_states + len(aux_states))
        y0[states.index('SUSCEPTIBLE')] = p.N * (1 - p.init_sym)
        y0[states.index('SYMPTOMATIC')] = p.N * p.init_sym
        return y0

    def outputs(self, y, time=None):
        """Derived outputs from a state array (n_times x n_states+aux, or a single state vector)."""
        p = self.pars
        y = np.atleast_2d(y)
        res = sc.objdict()
        if time is not None:
            res.time = time
        for i, s in enumerate(states + aux_states):
            res[s] = y[:, i]
        res.N = y[:, :n_states].sum(axis=1)
        res.prev_active   = (res.ASYMPTOMATIC + res.SYMPTOMATIC) / res.N                       # point prevalence, all ages
        res.frac_asymp    = res.ASYMPTOMATIC / (res.ASYMPTOMATIC + res.SYMPTOMATIC)             # of prevalent ASY+SYM
        res.incidence_py  = (p.inf_asy * res.INFECTION + p.non_asy * res.NON0) / res.N          # first ASY entry per episode, per person-year
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
        return f - x * f[:n_states].sum()

    def equilibrium(self, burn=2000):
        """
        Endemic equilibrium as population proportions, returned as a state vector (sums to 1 over the real
        compartments). Integrates the proportion system for ``burn`` years from a 1% SYMPTOMATIC seed
        (regardless of ``init_sym``), then polishes with a root solve.
        """
        x0 = np.zeros(n_states + len(aux_states))
        x0[states.index('SUSCEPTIBLE')] = 0.99
        x0[states.index('SYMPTOMATIC')] = 0.01
        x = odeint(lambda x, t: self._norm_derivs(x), x0, [0, burn], rtol=1e-10, atol=1e-12, mxstep=100_000)[-1]
        if x[states.index('ASYMPTOMATIC')] + x[states.index('SYMPTOMATIC')] > 1e-12:
            x = fsolve(self._norm_derivs, x, xtol=1e-13)
        return x


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
