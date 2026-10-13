"""
tbsim.TB variant whose natural-history step is sampled exactly from the continuous-time Markov chain.
"""

import numpy as np
from scipy.linalg import expm
import starsim as ss
import tbsim
from tbsim import TBS
from .ode import latent_classes

__all__ = ['TBExact', 'TBExactMassAction']

# Within-step chain: transient states, then absorbing ones. CLEARED is split by pathway so the right
# rr_reinfection can be assigned (cleared from latent vs recovered from NON_INFECTIOUS)
_chain = [TBS.INFECTION, TBS.NON_INFECTIOUS, TBS.ASYMPTOMATIC, TBS.SYMPTOMATIC, 'CLE_inf', 'CLE_non', TBS.DEAD]
_end_state = {0: TBS.INFECTION, 1: TBS.NON_INFECTIOUS, 2: TBS.ASYMPTOMATIC, 3: TBS.SYMPTOMATIC, 4: TBS.CLEARED, 5: TBS.CLEARED, 6: TBS.DEAD}


class TBExact(tbsim.TB):
    """
    tbsim.TB whose natural-history transitions are drawn from the exact per-step transition matrix expm(Q*dt).

    Standard tbsim.TB applies competing-risk transitions immediately and in a fixed order (INFECTION,
    NON_INFECTIOUS, ASYMPTOMATIC, SYMPTOMATIC), so agents can cascade forward through several states within a
    step but not backward. Here each agent's end-of-step state is sampled from the matrix exponential of the
    generator instead, which is exact for constant rates within the step (including multiple jumps). Transmission
    and bookkeeping are unchanged.

    Progression by time since infection (k_asy = k_non = k > 0) uses the Poisson-mixture latent classes of
    tbode.latent_classes: each agent entering INFECTION is assigned a class m, with probability Poisson(m; r0/k),
    and then progresses at the constant rate m*k, so each class has its own constant-rate generator. This gives
    each agent the same distribution of progression times as the hazard r0*exp(-k*tau) in tbsim.TB. The class is
    stored in ``latent_class`` (-1 until drawn); ``ti_infected`` is still set but isn't used for progression.

    Incidence is also counted from the sampled transitions, not only from end-of-step snapshots: for each agent starting
    the step in INFECTION or NON_INFECTIOUS, whether its path passed through ASYMPTOMATIC or SYMPTOMATIC is drawn given its
    start and end states, from an augmented chain in which entering either state marks the path. A new episode (first
    ASYMPTOMATIC/SYMPTOMATIC entry since leaving SUSCEPTIBLE, INFECTION or CLEARED, as in tb_LAI_TPT's IncidenceTracker)
    is counted in ``results.new_episodes`` and flagged in ``new_episode``; ``been_active`` holds the per-agent episode flag.
    The visit draw uses its own random stream after the end state is drawn, so the dynamics are unchanged.

    Requires rr_activation/rr_clearance/rr_death = 1 and k_asy = k_non.

    Args:
        infection_timing (str): how agents infected during a step are stepped. 'start' (default): through the
            full step's expm(Q*dt), i.e. as if infected at the start of the step. 'within': each gets a uniform
            infection time u (fraction of the step already passed) and is stepped through expm(Q*(1-u)*dt).
        kwargs: passed to tbsim.TB
    """
    def __init__(self, pars=None, infection_timing='start', **kwargs):
        super().__init__(pars, **kwargs)
        if infection_timing not in ('start', 'within'):
            raise ValueError(f"infection_timing must be 'start' or 'within', not {infection_timing!r}")
        self.infection_timing = infection_timing
        self.define_states(
            ss.IntArr('latent_class', default=-1),       # Poisson-mixture class while in INFECTION (k > 0); -1 = not yet drawn
            ss.BoolState('been_active', default=False),  # Entered ASYMPTOMATIC/SYMPTOMATIC in the current episode
            ss.BoolState('new_episode', default=False),  # Started an incident episode during this step
        )
        self._rng_tinf = ss.random(name='tb_rng_tinf')    # Infection time within the step ('within' only)
        self._rng_class = ss.random(name='tb_rng_class')  # Latent class at infection (k > 0 only)
        self._rng_visit = ss.random(name='tb_rng_visit')  # Whether the path visited ASYMPTOMATIC/SYMPTOMATIC within the step (incidence only)
        self._new_inf = ss.uids()  # Agents entering INFECTION from an uninfected state during this step
        return

    def set_prognoses(self, uids, sources=None):
        uids = ss.uids(uids)
        new = uids[~self.infected[uids]]  # Reinfected latent agents keep their state, so they aren't new infections
        super().set_prognoses(uids, sources)
        self.latent_class[uids] = -1  # (Re)infection restarts the clock: the class is drawn afresh in step_transitions
        self._new_inf = self._new_inf.concatenate(new)
        return

    def step(self):
        self._new_inf = ss.uids()  # Transmission in super().step() fills this via set_prognoses
        return super().step()

    def init_results(self):
        super().init_results()
        self.define_results(ss.Result('new_episodes', dtype=int, scale=False, label='New incident episodes (from transitions)'))
        return

    def _rate_dt(self, r):
        """Rate per timestep, as in tbsim.TB.transition."""
        return float(r.value) * (self.sim.t.dt / r.unit)

    def _generators(self):
        """
        Per-step generators over _chain, one per latent class (a single one when k = 0), plus the class weights.
        The classes differ only in their INFECTION -> NON_INFECTIOUS / ASYMPTOMATIC rates.
        """
        p = self.pars
        if p.k_asy != p.k_non:
            raise NotImplementedError('TBExact requires k_asy = k_non')
        for name in ['rr_activation', 'rr_clearance', 'rr_death']:
            if np.any(np.asarray(getattr(self, name).values) != 1):
                raise NotImplementedError(f'TBExact requires {name} = 1')

        r = {k: self._rate_dt(p[k]) for k in ['inf_cle', 'inf_non', 'inf_asy', 'non_rec', 'non_asy', 'asy_non', 'asy_sym', 'sym_asy', 'sym_dead']}
        if p.k_asy > 0:
            weights, rate_non, rate_asy = latent_classes(r['inf_non'], r['inf_asy'], p.k_asy * self.sim.t.dt_year)  # k per step
        else:
            weights, rate_non, rate_asy = np.ones(1), np.array([r['inf_non']]), np.array([r['inf_asy']])
        Qs = np.zeros((len(weights), len(_chain), len(_chain)))
        Qs[:, 0, 4], Qs[:, 0, 1], Qs[:, 0, 2] = r['inf_cle'], rate_non, rate_asy
        Qs[:, 1, 5], Qs[:, 1, 2] = r['non_rec'], r['non_asy']
        Qs[:, 2, 1], Qs[:, 2, 3] = r['asy_non'], r['asy_sym']
        Qs[:, 3, 2], Qs[:, 3, 6] = r['sym_asy'], r['sym_dead']
        idx = np.arange(len(_chain))
        Qs[:, idx, idx] = -Qs.sum(axis=2)
        return Qs, weights

    def step_transitions(self):
        Qs, weights = self._generators()
        cumP = np.cumsum(expm(Qs), axis=2)  # (n_classes, n_chain, n_chain)

        rngs = [self._rng_inf, self._rng_non, self._rng_asy, self._rng_sym]
        start_masks = [self.latent, self.non_infectious, self.asymptomatic, self.symptomatic]
        moves = []
        inf_cls, inf_new = None, None
        for i, (mask, rng) in enumerate(zip(start_masks, rngs)):
            uids = mask.uids
            if len(uids):
                u = rng.rvs(uids)
                if i == 0:
                    cls = self._latent_classes(uids, weights)
                    cum = cumP[cls, 0, :]
                    inf_cls = cls
                    if self.infection_timing == 'within':
                        cum, inf_new = self._new_inf_cum(uids, cls, Qs, cum)
                else:
                    cum = np.broadcast_to(cumP[0, i], (len(uids), len(_chain)))
                dest = np.minimum((u[:, None] >= cum).sum(axis=1), len(_chain) - 1)
                moves.append((i, uids, dest))

        self._count_episodes(moves, Qs, inf_cls, inf_new)
        self._apply_moves(moves)
        return

    @staticmethod
    def _augment(Qs):
        """
        Augmented generators over [_chain unmarked, _chain marked]: entering ASYMPTOMATIC or SYMPTOMATIC moves the path
        to the marked copy, where it stays. Then expm(Q_aug)[start, n + end] = P(visited ASY/SYM and ended in end).
        """
        C, n, _ = Qs.shape
        Qa = np.zeros((C, 2 * n, 2 * n))
        off = Qs.copy()
        idx = np.arange(n)
        off[:, idx, idx] = 0
        active = np.isin(np.arange(n), [2, 3])  # ASYMPTOMATIC, SYMPTOMATIC
        Qa[:, :n, :n] = off * ~active[None, None, :]  # Unmarked -> unmarked (non-active destinations)
        Qa[:, :n, n:] = off * active[None, None, :]   # Unmarked -> marked (entering ASY/SYM)
        Qa[:, n:, n:] = off                           # Marked stays marked
        diag = Qs[:, idx, idx]
        Qa[:, idx, idx] = diag
        Qa[:, n + idx, n + idx] = diag
        return Qa

    def _count_episodes(self, moves, Qs, inf_cls, inf_new):
        """
        For agents starting in INFECTION or NON_INFECTIOUS, draw whether the path visited ASYMPTOMATIC/SYMPTOMATIC given
        start and end states; count new episodes; update the episode flag. Uses only _rng_visit, so dynamics are unchanged.
        """
        n = len(_chain)
        Qa = self._augment(Qs)
        Pa = expm(Qa)  # (n_classes, 2n, 2n), full step
        self.new_episode[:] = False
        cand_uids, cand_pv = [], []
        for i, uids, dest in moves:
            if i not in (0, 1):
                continue
            if i == 0:
                rows = Pa[inf_cls, 0, :]  # Per agent, by latent class
                if inf_new is not None and inf_new[0].any():
                    is_new, remaining = inf_new
                    rows = rows.copy()
                    rows[is_new] = expm(Qa[inf_cls[is_new]] * remaining[:, None, None])[:, 0, :]
            else:
                rows = np.broadcast_to(Pa[0, 1, :], (len(uids), 2 * n))
            ar = np.arange(len(uids))
            marked, unmarked = rows[ar, n + dest], rows[ar, dest]
            with np.errstate(invalid='ignore', divide='ignore'):
                pv = np.where(marked + unmarked > 0, marked / (marked + unmarked), 0)
            cand_uids.append(uids)
            cand_pv.append(pv)
        visited = {}
        if cand_uids:
            all_uids = ss.uids.concatenate(cand_uids)
            all_pv = np.concatenate(cand_pv)
            v = self._rng_visit.rvs(all_uids) < all_pv
            k = 0
            for uids in cand_uids:
                visited[id(uids)] = v[k:k + len(uids)]
                k += len(uids)

        n_new = 0
        for i, uids, dest in moves:
            vis = visited.get(id(uids), np.zeros(len(uids), bool))
            old = np.asarray(self.been_active[uids])
            new_ep = vis & ~old
            n_new += int(new_ep.sum())
            self.new_episode[uids[new_ep]] = True
            start_active = i in (2, 3)
            end_active = np.isin(dest, [2, 3])
            ended = np.isin(dest, [4, 5, 6])  # CLEARED (either pathway) or DEAD: the episode is over
            self.been_active[uids] = (old | vis | start_active | end_active) & ~ended
        self.results['new_episodes'][self.ti] = n_new
        return

    def _latent_classes(self, uids, weights):
        """Latent class of each INFECTION agent, drawing it (once per infection) for agents that don't have one yet."""
        if len(weights) == 1:
            return np.zeros(len(uids), dtype=int)
        cls = np.asarray(self.latent_class[uids])
        todo = cls < 0
        if todo.any():
            draw = np.searchsorted(np.cumsum(weights), self._rng_class.rvs(uids[todo]), side='right')
            cls[todo] = np.minimum(draw, len(weights) - 1)
            self.latent_class[uids[todo]] = cls[todo]
        return cls

    def _new_inf_cum(self, uids, cls, Qs, cum):
        """
        Cumulative end-state probabilities for INFECTION agents, with those infected this step stepped through
        only the rest of the step after a uniform infection time: row 0 of expm(Q*(1-u)) per agent, with Q its class's generator.
        Also returns (is_new, remaining) for the incidence count, or None if nobody was infected this step.
        """
        is_new = np.isin(uids, self._new_inf)
        if not is_new.any():
            return cum, None
        new = uids[is_new]
        remaining = 1 - self._rng_tinf.rvs(new)  # Fraction of the step left after infection
        rows = expm(Qs[cls[is_new]] * remaining[:, None, None])[:, 0, :]
        cum = cum.copy()
        cum[is_new] = np.cumsum(rows, axis=1)
        return cum, (is_new, remaining)

    def _apply_moves(self, moves):
        """Apply the sampled end states after sampling everyone, so each agent moves exactly once from its start-of-step state."""
        p = self.pars
        for _, uids, dest in moves:
            for j, s in _end_state.items():
                u = uids[dest == j]
                if not len(u):
                    continue
                self.state[u] = s
                if j == 4:
                    self._enter_cleared(u, p.rr_reinfection_cleared)
                elif j == 5:
                    self._enter_cleared(u, p.rr_reinfection_rec)
                elif s == TBS.ASYMPTOMATIC:
                    self.ti_asymp[u] = self.ti
        return


class TBExactMassAction(TBExact):
    """
    TBExact with mass-action transmission in place of network transmission (use with no networks).

    Each step, every susceptible agent j is infected with probability 1 - exp(-beta*w*rel_sus_j*dt) ('exp') or
    min(1, beta*w*rel_sus_j*dt) ('linear'), where beta is the ODE (mass-action) beta and w = (trans_asymp*ASY + SYM)/N
    over living agents. 'exp' is what network transmission gives per step; 'linear' gives the ODE's infection flow
    beta*w*S*dt for a susceptible pool fixed at the start of the step.

    Args:
        beta_ode (float): mass-action beta (per year), as in TBODE
        infection_prob (str): 'exp' or 'linear'
        kwargs: passed to TBExact
    """
    def __init__(self, pars=None, beta_ode=None, infection_prob='exp', **kwargs):
        super().__init__(pars, **kwargs)
        if beta_ode is None:
            raise ValueError('TBExactMassAction needs beta_ode')
        if infection_prob not in ('exp', 'linear'):
            raise ValueError(f"infection_prob must be 'exp' or 'linear', not {infection_prob!r}")
        self.beta_ode = beta_ode
        self.infection_prob = infection_prob
        self._rng_ma = ss.random(name='tb_rng_mass_action')
        return

    def infect(self):
        """New infections this step under mass action; returns (new_cases, sources, networks) as ss.Infection.infect does."""
        alive = self.sim.people.auids
        infectious = np.asarray(self.infectious[alive]) * np.asarray(self.rel_trans[alive])
        w = infectious.sum() / len(alive)
        sus = alive[np.asarray(self.susceptible[alive]) & (np.asarray(self.rel_sus[alive]) > 0)]
        if w == 0 or not len(sus):
            return ss.uids(), None, None
        x = self.beta_ode * w * np.asarray(self.rel_sus[sus]) * self.t.dt_year
        p = 1 - np.exp(-x) if self.infection_prob == 'exp' else np.minimum(x, 1)
        new = sus[self._rng_ma.rvs(sus) < p]
        return new, None, None
