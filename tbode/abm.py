"""
Build agent-based tbsim.TB runs matched to a TBODE configuration.
"""

import numpy as np
import sciris as sc
import starsim as ss
import tbsim
from tbsim import TBS
from .ode import states as ode_states, n_states

__all__ = ['beta_edge_from_ode', 'StateCounts', 'FluxCounts', 'make_abm', 'init_from_ode', 'abm_cohort_progression']


def beta_edge_from_ode(beta_ode, n_contacts, dt):
    """
    Per-edge ABM beta (per year) giving the same per-step force of infection as the ODE's mass-action beta.

    On ss.RandomNet(n_contacts=K, dur=0), K is the mean number of contacts per agent (~K/2 edges per agent,
    with transmission both ways along each edge; checked empirically for starsim 3.6.1). Each contact with an
    infectious agent transmits with probability q = 1 - exp(-beta_edge*dt), so the per-step hazard is ~K*q*w,
    where w is the (trans_asymp-weighted) infectious fraction. The ODE hazard over dt is beta_ode*w*dt, so
    q = beta_ode*dt/K.

    Args:
        beta_ode (float): ODE beta (per year)
        n_contacts (int): RandomNet n_contacts (K)
        dt (float): ABM timestep in years
    """
    q = beta_ode * dt / n_contacts
    if q >= 1:
        raise ValueError(f'beta_ode={beta_ode} is too large for n_contacts={n_contacts} at dt={dt}')
    return -np.log(1 - q) / dt


class StateCounts(ss.Analyzer):
    """
    Record, each timestep, the number of living agents in each TBS state and the number of incident episodes.

    Incidence follows tb_LAI_TPT's IncidenceTracker (sym_only=False): an episode is counted when an agent first
    enters ASYMPTOMATIC or SYMPTOMATIC after being in SUSCEPTIBLE, INFECTION or CLEARED (so ASYMPTOMATIC <->
    NON_INFECTIOUS cycling within an episode is not double counted).
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.define_states(ss.BoolState('been_active', default=False))  # entered ASY/SYM during the current episode
        return

    def init_results(self):
        super().init_results()
        self.define_results(
            ss.Result('n_alive', dtype=int, scale=False),
            ss.Result('new_incident', dtype=int, scale=False),
            *[ss.Result(s.name, dtype=int, scale=False) for s in TBS],
        )
        return

    def step(self):
        tb = self.sim.get_tb()
        alive = self.sim.people.alive
        st = tb.state
        active = st.isin((TBS.ASYMPTOMATIC, TBS.SYMPTOMATIC))
        self.results.new_incident[self.ti] = np.count_nonzero(active & ~self.been_active & alive)
        self.been_active[active] = True
        self.been_active[st.isin((TBS.SUSCEPTIBLE, TBS.INFECTION, TBS.CLEARED))] = False

        counts = np.bincount(np.asarray(st[alive], dtype=int), minlength=len(TBS))
        self.results.n_alive[self.ti] = alive.count()
        for s in TBS:
            self.results[s.name][self.ti] = counts[int(s)]
        return


class FluxCounts(ss.Analyzer):
    """
    Count state-to-state transitions of living agents over each timestep, plus the infections an ODE force of
    infection would predict from the same start-of-step state.

    Call ``snapshot()`` once after ``init_from_ode`` (before running) so the first step has a start state.
    After the run, ``self.flows`` is an (n_steps, n_tbs, n_tbs) array of counts from start state (row) to end
    state (column), ``self.stock`` the (n_steps, n_tbs) start-of-step stocks, and ``self.exp_inf_S`` /
    ``self.exp_inf_C`` the expected infections among SUSCEPTIBLE / CLEARED agents under the ODE.

    Args:
        beta_ode (float): ODE beta (per year), used for the expected infections
    """
    def __init__(self, beta_ode, **kwargs):
        super().__init__(**kwargs)
        self.beta_ode = beta_ode
        self.define_states(ss.FloatArr('prev_state', default=TBS.SUSCEPTIBLE), ss.BoolState('prev_alive', default=False))
        return

    def init_post(self):
        super().init_post()
        nt, ns = len(self.t), len(TBS)
        self.flows = np.zeros((nt, ns, ns), dtype=np.int64)
        self.stock = np.zeros((nt, ns), dtype=np.int64)
        self.exp_inf_S = np.zeros(nt)
        self.exp_inf_C = np.zeros(nt)
        self._pending = None
        return

    def snapshot(self):
        """Record the start-of-step state for the next step."""
        tb = self.sim.get_tb()
        alive = self.sim.people.alive
        self.prev_state[:] = tb.state[:]
        self.prev_alive[:] = alive[:]
        st = np.asarray(tb.state[alive], dtype=int)
        n = len(st)
        asy, sym = np.count_nonzero(st == TBS.ASYMPTOMATIC), np.count_nonzero(st == TBS.SYMPTOMATIC)
        foi = self.beta_ode * (tb.pars.trans_asymp * asy + sym) / n
        dt = self.sim.t.dt_year
        cleared = (tb.state == TBS.CLEARED) & alive
        rr = np.asarray(tb.rr_reinfection[cleared])
        self._pending = (np.count_nonzero(st == TBS.SUSCEPTIBLE) * (1 - np.exp(-foi * dt)), np.sum(1 - np.exp(-foi * rr * dt)))
        return

    def step(self):
        tb = self.sim.get_tb()
        ns = len(TBS)
        u = self.prev_alive.uids
        a = np.asarray(self.prev_state[u], dtype=int)
        b = np.asarray(tb.state[u], dtype=int)
        b = np.where(np.asarray(self.sim.people.alive[u]), b, int(TBS.DEAD))  # agents removed this step count as DEAD
        self.flows[self.ti] = np.bincount(a * ns + b, minlength=ns * ns).reshape(ns, ns)
        self.stock[self.ti] = np.bincount(a, minlength=ns)
        self.exp_inf_S[self.ti], self.exp_inf_C[self.ti] = self._pending
        self.snapshot()
        return


def make_abm(ode_pars, seed=0, n_agents=20_000, n_contacts=10, start=1900, stop=2100, dt=1/12, label=None, analyzers=None, tb_class=None, **tb_kwargs):
    """
    Build a tbsim.Sim matched to a TBODE parameter set.

    Matching choices: ss.RandomNet(dur=0) for mass action (beta converted by beta_edge_from_ode); constant
    ss.Births and ss.Deaths at rate mu (so TB deaths are *not* replaced: compare with TBODE(replace_tb_deaths=False));
    initial condition of init_sym in SYMPTOMATIC and the rest SUSCEPTIBLE (or use init_from_ode); reinfection of
    INFECTION and NON_INFECTIOUS turned off, as in the ODE (with constant hazards this changes only the random draws).

    Args:
        ode_pars (dict): TBODE parameters (e.g. ``TBODE().pars``)
        seed (int): random seed
        n_agents (int): number of agents
        n_contacts (int): RandomNet n_contacts
        start, stop (int): start and stop years
        dt (float): timestep in years
        label (str): sim label
        analyzers (list): additional analyzers (StateCounts is always included)
        tb_class (class): TB class to use (default tbsim.TB; e.g. tbode.TBExact)
        tb_kwargs: overrides passed to tbsim.TB
    """
    p = sc.objdict(ode_pars)
    rate_names = ['inf_cle', 'inf_non', 'inf_asy', 'non_rec', 'non_asy', 'asy_non', 'asy_sym', 'sym_asy', 'sym_dead']
    tb_pars = dict(
        beta=ss.peryear(beta_edge_from_ode(p.beta, n_contacts, dt)),
        init_prev=ss.bernoulli(0.0),
        init_prev_active=ss.bernoulli(p.init_sym),
        trans_asymp=p.trans_asymp,
        rr_reinfection_cleared=p.rr_reinfection_cleared,
        rr_reinfection_rec=p.rr_reinfection_rec,
        rr_reinfection_treat=p.rr_reinfection_treat,
        rr_reinfection_inf=0.0,
        rr_reinfection_non=0.0,
        **{k: ss.peryear(p[k]) for k in rate_names},
    )
    tb_pars.update(tb_kwargs)
    tb = (tb_class or tbsim.TB)(**tb_pars)

    sim = tbsim.Sim(
        sim_pars=dict(start=ss.date(f'{start}-01-01'), stop=ss.date(f'{stop}-01-01'), dt=ss.years(dt), rand_seed=seed),
        n_agents=n_agents,
        tb_model=tb,
        networks=ss.RandomNet(n_contacts=n_contacts, dur=0),
        demographics=[ss.Births(birth_rate=ss.peryear(1000 * p.mu)), ss.Deaths(death_rate=ss.peryear(1000 * p.mu))],
        analyzers=[StateCounts(name='statecounts')] + sc.tolist(analyzers),
        label=label,
        verbose=0,
    )
    return sim


def init_from_ode(sim, x, seed=0):
    """
    Assign initial TB states to an initialised (not yet run) sim from TBODE state proportions ``x``.

    Each ODE compartment maps to an ABM state; CLEARED/RECOVERED/TREATED map to CLEARED with the matching
    per-agent rr_reinfection. The share of NON_INFECTIOUS agents that has already been ASYMPTOMATIC in this
    episode (NON - NON0 in the ODE) is marked as such for incidence counting. Build the sim with
    ``init_sym=0`` in its ODE pars so no other seeding happens.
    """
    tb = sim.get_tb()
    p = tb.pars
    uids = sim.people.auids
    rng = np.random.default_rng(seed)
    props = np.clip(np.asarray(x[:n_states], dtype=float), 0, None)
    props /= props.sum()
    comp = rng.choice(n_states, size=len(uids), p=props)

    to_abm = dict(SUSCEPTIBLE=TBS.SUSCEPTIBLE, INFECTION=TBS.INFECTION, CLEARED=TBS.CLEARED, RECOVERED=TBS.CLEARED,
                  NON_INFECTIOUS=TBS.NON_INFECTIOUS, ASYMPTOMATIC=TBS.ASYMPTOMATIC, SYMPTOMATIC=TBS.SYMPTOMATIC,
                  TREATMENT=TBS.TREATMENT, TREATED=TBS.CLEARED)
    rr = dict(CLEARED=p.rr_reinfection_cleared, RECOVERED=p.rr_reinfection_rec, TREATED=p.rr_reinfection_treat)
    for i, name in enumerate(ode_states):
        u = uids[comp == i]
        tb.state[u] = to_abm[name]
        if name in rr:
            tb.rr_reinfection[u] = rr[name]
        if to_abm[name] not in (TBS.SUSCEPTIBLE,):
            tb.ever_infected[u] = True
        if to_abm[name] in (TBS.INFECTION, TBS.NON_INFECTIOUS, TBS.ASYMPTOMATIC, TBS.SYMPTOMATIC):
            tb.ti_infected[u] = 0

    sc_an = sim.analyzers.statecounts
    non = uids[comp == ode_states.index('NON_INFECTIOUS')]
    frac_been = 1 - x[n_states] / x[ode_states.index('NON_INFECTIOUS')] if x[ode_states.index('NON_INFECTIOUS')] > 0 else 0
    sc_an.been_active[non[rng.random(len(non)) < frac_been]] = True
    sc_an.been_active[tb.state.isin((TBS.ASYMPTOMATIC, TBS.SYMPTOMATIC))] = True
    tb.step_bookkeeping()  # derive susceptible/infected/rel_sus/rel_trans from the assigned states
    return sim


def abm_cohort_progression(ode_pars, dt=1/12, n_agents=20_000, seed=0, years=10, tb_class=None, **tb_kwargs):
    """
    ABM version of ode_cohort_progression: everyone starts latent at t=0, no transmission, no demographics,
    no TB death. Returns {'years', 'B'} with B the cumulative fraction ever ASYMPTOMATIC at each whole year.
    """
    p = sc.objdict(ode_pars)
    rate_names = ['inf_cle', 'inf_non', 'inf_asy', 'non_rec', 'non_asy', 'asy_non', 'asy_sym', 'sym_asy']
    tb = (tb_class or tbsim.TB)(beta=ss.peryear(0), init_prev=ss.bernoulli(1.0), sym_dead=ss.peryear(0),
                  **{k: ss.peryear(p[k]) for k in rate_names}, **tb_kwargs)

    class EverAsy(ss.Analyzer):
        def __init__(self):
            super().__init__(name='everasy')
            self.define_states(ss.BoolState('ever_asy', default=False))

        def init_results(self):
            super().init_results()
            self.define_results(ss.Result('frac', dtype=float, scale=False))

        def step(self):
            tb = self.sim.get_tb()
            self.ever_asy[tb.state.isin((TBS.ASYMPTOMATIC, TBS.SYMPTOMATIC))] = True
            self.results.frac[self.ti] = self.ever_asy.mean()

    sim = tbsim.Sim(sim_pars=dict(start=ss.years(0), stop=ss.years(years), dt=ss.years(dt), rand_seed=seed),
                    n_agents=n_agents, tb_model=tb, networks=[], demographics=[], analyzers=EverAsy(), verbose=0)
    sim.run()
    t = np.asarray(sim.t.yearvec) + dt  # result at index ti is recorded at the end of step ti
    frac = np.asarray(sim.results.everasy.frac)
    yrs = np.arange(1, years + 1)
    return sc.objdict(years=yrs, B=np.interp(yrs, t, frac))
