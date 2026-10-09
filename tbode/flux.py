"""
Per-step flux checks: compare the ABM's realised per-step transitions (from FluxCounts) with the exact
continuous-time model over one step, expm(Q*dt).
"""

import numpy as np
import pandas as pd
import sciris as sc
from scipy.linalg import expm
from tbsim import TBS
from .ode import TBODE
from .abm import FluxCounts, make_abm, init_from_ode

__all__ = ['nh_states', 'nh_generator', 'run_flux', 'flux_tables']

# Natural-history states over which per-step transitions are compared
nh_states = [TBS.INFECTION, TBS.CLEARED, TBS.NON_INFECTIOUS, TBS.ASYMPTOMATIC, TBS.SYMPTOMATIC, TBS.DEAD]


def nh_generator(p):
    """Continuous-time generator (per year) over nh_states for the natural history (no infection), incl. background death mu."""
    idx = {s: i for i, s in enumerate(nh_states)}
    Q = np.zeros((len(nh_states), len(nh_states)))
    rates = [(TBS.INFECTION, TBS.CLEARED, p.inf_cle), (TBS.INFECTION, TBS.NON_INFECTIOUS, p.inf_non), (TBS.INFECTION, TBS.ASYMPTOMATIC, p.inf_asy),
             (TBS.NON_INFECTIOUS, TBS.CLEARED, p.non_rec), (TBS.NON_INFECTIOUS, TBS.ASYMPTOMATIC, p.non_asy),
             (TBS.ASYMPTOMATIC, TBS.NON_INFECTIOUS, p.asy_non), (TBS.ASYMPTOMATIC, TBS.SYMPTOMATIC, p.asy_sym),
             (TBS.SYMPTOMATIC, TBS.ASYMPTOMATIC, p.sym_asy), (TBS.SYMPTOMATIC, TBS.DEAD, p.sym_dead)]
    rates += [(s, TBS.DEAD, p.mu) for s in nh_states if s != TBS.DEAD]
    for a, b, r in rates:
        Q[idx[a], idx[b]] += r
    Q[np.diag_indices_from(Q)] = -Q.sum(axis=1)
    return Q


def run_flux(beta, dt, seed, n_agents=50_000, start=1900, stop=1910, tb_class=None, ode_kw=None):
    """
    Run one ABM replicate from the ODE equilibrium with a FluxCounts analyzer; return its summed flows and stocks
    and the expected infections among SUSCEPTIBLE / CLEARED agents.
    """
    ode_kw = sc.mergedicts(dict(replace_tb_deaths=False, init_sym=0.0), ode_kw)
    ode = TBODE(beta=beta, **ode_kw)
    flux = FluxCounts(beta_ode=beta, name='flux')
    sim = make_abm(ode.pars, seed=seed, n_agents=n_agents, start=start, stop=stop, dt=dt, analyzers=flux, tb_class=tb_class)
    sim.init()
    init_from_ode(sim, ode.equilibrium(), seed=seed)
    f = sim.analyzers.flux
    f.snapshot()
    sim.run()
    return dict(flows=f.flows.sum(axis=0), stock=f.stock.sum(axis=0), exp_S=f.exp_inf_S.sum(), exp_C=f.exp_inf_C.sum())


def flux_tables(outs, pars, dt, label):
    """
    Summarise run_flux outputs (one list of replicates at one dt) into natural-history and transmission tables.

    Natural-history rows give, for each flow A->B, the ABM per-step probability (count / stock of A) as a rate
    (divided by dt), the exact continuous-time value expm(Q*dt)[A, B]/dt, their ratio, and the binomial standard
    error of the ratio. 'A stays' rows give per-step probabilities instead of rates. Transmission rows compare realised with ODE-expected infections from SUSCEPTIBLE and CLEARED.
    """
    flows = sum(o['flows'] for o in outs)
    stock = sum(o['stock'] for o in outs)
    P_exact = expm(nh_generator(sc.objdict(pars)) * dt)
    rows = []
    for ia, a in enumerate(nh_states[:-1]):
        if a == TBS.CLEARED:
            continue
        for ib, b in enumerate(nh_states):
            p_ex = P_exact[ia, ib]
            n = int(flows[a, b])
            if a != b and p_ex < 1e-9 and n == 0:
                continue
            p_abm = n / stock[a]
            name = f'{a.name} stays' if a == b else f'{a.name}->{b.name}'
            se = np.sqrt(p_ex * (1 - p_ex) / stock[a]) / p_ex if p_ex > 0 else np.nan
            per = 1 if a == b else dt  # 'stays' rows are reported as per-step probabilities, flows as rates per year
            rows.append(dict(dt=label, flow=name, n=n, abm_rate=p_abm / per, exact_rate=p_ex / per,
                             ratio=p_abm / p_ex if p_ex > 0 else np.nan, se_ratio=se))

    dead = int(TBS.DEAD)
    inf_S = flows[TBS.SUSCEPTIBLE].sum() - flows[TBS.SUSCEPTIBLE, TBS.SUSCEPTIBLE] - flows[TBS.SUSCEPTIBLE, dead]
    inf_C = flows[TBS.CLEARED].sum() - flows[TBS.CLEARED, TBS.CLEARED] - flows[TBS.CLEARED, dead]
    exp_S = sum(o['exp_S'] for o in outs)
    exp_C = sum(o['exp_C'] for o in outs)
    trans = dict(dt=label, infections_S=int(inf_S), expected_S=exp_S, ratio_S=inf_S / exp_S,
                 infections_C=int(inf_C), expected_C=exp_C, ratio_C=inf_C / exp_C,
                 se_S=1 / np.sqrt(max(inf_S, 1)), se_C=1 / np.sqrt(max(inf_C, 1)))
    return pd.DataFrame(rows), trans
