"""
Matched runs of tbsim v4's TB model: the same declared model as an ODE or as agents, with the settings of the v0.10 experiments.
"""
import numpy as np
import sciris as sc
import starsim as ss
import tbsim

__all__ = ['states', 'state_map', 'beta_edge_from_ode', 'make_sim', 'shares', 'odepars_to_tb']

# The model's states, and the matching compartment of the v0.10 ODE (tbode.TBODE)
states = ['susceptible', 'infection', 'cleared', 'recovered', 'non_infectious', 'asymptomatic', 'symptomatic', 'treatment', 'treated']
state_map = {s: s.upper() for s in states}

# TBODE parameters that tbsim.TB doesn't take, or takes differently
_not_tb = {'N', 'mu', 'replace_tb_deaths', 'init_sym', 'k_asy', 'k_non', 'latent_tol', 'rr_reinfection_inf', 'rr_reinfection_non'}


def odepars_to_tb(pars):
    """ The tbsim.TB keyword arguments for a dict of TBODE parameters (rates per year), leaving out those it handles otherwise (see make_sim) """
    return {k: v for k, v in pars.items() if k not in _not_tb}


def beta_edge_from_ode(beta_ode, n_contacts, dt):
    """
    The per-contact beta (per year) that gives agents on ss.RandomNet the ODE's infections per step, as in v0.10 (tbode.beta_edge_from_ode)

    On ss.RandomNet(n_contacts=K) (dur=0, fresh contacts each step), K is the mean degree (checked in exp 12 for Starsim v4). A
    contact with an infectious person transmits in a step with probability q = 1 − exp(−beta_edge·dt), so the step's hazard is about
    K·q·w, for w the infectious share (weighted by trans_asymp). Setting it to the ODE's beta_ode·w·dt gives q = beta_ode·dt/K.
    v4's own mapping, beta_ode/K (as its ODE approximation of the network assumes), gives about b·dt/2 fewer infections (4% at monthly dt).

    Args:
        beta_ode (float): the ODE's beta (per year)
        n_contacts (int): K
        dt (float): the step, in years
    """
    q = beta_ode*dt/n_contacts
    if q >= 1:
        raise ValueError(f'beta_ode={beta_ode} is too large for n_contacts={n_contacts} at dt={dt}')
    return -np.log1p(-q)/dt


def make_sim(beta, solver='ode', dt=ss.months(1), start=1900, stop=2000, init=None, n_agents=None, total_pop=1e5, n_contacts=None, edge_beta='step',
             mu=1/70, replace_tb_deaths=False, seed=0, verbose=0, analyzers=None, **tb_kw):
    """
    A sim of tbsim.TB, set up as the v0.10 experiments were

    Args:
        beta (float): the ODE's beta, per year (homogeneous mixing); on a network, each contact gets beta/n_contacts
        solver (str): 'ode' or 'abm' (or any Starsim solver)
        dt (ss.dur): the step
        start, stop (int): years
        init (dict): the share of people in each state at the start (default: tbsim's, 5% latent)
        n_agents (int): agents, for solver='abm'
        total_pop (float): the population (as an ODE, or the people the agents stand for)
        n_contacts (int): with a number, an ss.RandomNet with this mean degree (fresh contacts each step); None: homogeneous mixing
        edge_beta (str): on a network, the beta of each contact: 'step' (beta_edge_from_ode, matching the ODE's infections per step, as in v0.10)
            or 'rate' (beta/n_contacts, v4's own mapping)
        mu (float): background mortality, per year, with births at the same per-person rate
        replace_tb_deaths (bool): True: births also replace each TB death, so N is constant (tbsim.TB's own mu, as the reference ODE);
            False: births only replace background deaths (ss.birth and ss.death), as in the v0.10 experiments (TBODE(replace_tb_deaths=False))
        seed (int): the random seed
        analyzers (list): analyzers for the sim
        tb_kw (dict): other parameters of tbsim.TB (rates per year)
    """
    kw = dict(start=start, stop=stop, dt=dt, solver=solver, rand_seed=seed, verbose=verbose)
    if analyzers is not None:
        kw['analyzers'] = analyzers
    if n_agents is not None:
        kw['n_agents'] = n_agents
    kw['total_pop'] = total_pop
    if n_contacts:
        kw['networks'] = ss.RandomNet(n_contacts=n_contacts)
        if edge_beta == 'step':
            beta = beta_edge_from_ode(beta, n_contacts, dt.years)
        elif edge_beta == 'rate':
            beta = beta/n_contacts
        else:
            raise ValueError(f"edge_beta must be 'step' or 'rate', not {edge_beta}")
    if replace_tb_deaths:
        tb = tbsim.TB(beta=beta, mu=mu, init=init, **tb_kw)
    else:
        tb = tbsim.TB(beta=beta, mu=0, init=init, **tb_kw)
        kw['demographics'] = [ss.birth(rate=ss.peryear(mu)), ss.death(rate=ss.peryear(mu))]
    return ss.Sim(tb, **kw)


def shares(sim):
    """
    The share of people in each state over time, and the outputs compared with the ODE (shares of everyone alive)

    Returns an objdict of arrays: `year`, `N`, each state (lower case), `prev_active` (asymptomatic + symptomatic),
    `frac_asymp`, and the flows per person per year (each flow at a time point is the number during the step that ends there)
    """
    r = sim.results.tb
    out = sc.objdict(year=np.asarray(sim.t.yearvec, dtype=float))
    N = sum(np.asarray(r[s].values, dtype=float) for s in states)
    out.N = N
    for s in states:
        out[s] = np.asarray(r[s].values, dtype=float)/N
    out.prev_active = out.asymptomatic + out.symptomatic
    with np.errstate(divide='ignore', invalid='ignore'):
        out.frac_asymp = out.asymptomatic/out.prev_active
    dt = sim.t.dt.years if hasattr(sim.t.dt, 'years') else float(sim.t.dt)
    N_start = np.concatenate([[np.nan], N[:-1]]) # People at the start of each step
    for key in ['transmission', 'latent_asy', 'worsening', 'tb_death']:
        out[f'{key}_py'] = np.asarray(r[key].values, dtype=float)/N_start/dt
    return out
