"""
tbsim.TB variant whose natural-history step is sampled exactly from the continuous-time Markov chain.
"""

import numpy as np
from scipy.linalg import expm
import starsim as ss
import tbsim
from tbsim import TBS

__all__ = ['TBExact']

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

    Only for the constant-rate case: rr_activation/rr_clearance/rr_death must be 1 and k_asy = k_non = 0.
    """
    def _rate_dt(self, r):
        """Rate per timestep, as in tbsim.TB.transition."""
        return float(r.value) * (self.sim.t.dt / r.unit)

    def step_transitions(self):
        p = self.pars
        if p.k_asy or p.k_non:
            raise NotImplementedError('TBExact supports constant hazards only (k_asy = k_non = 0)')
        for name in ['rr_activation', 'rr_clearance', 'rr_death']:
            if np.any(np.asarray(getattr(self, name).values) != 1):
                raise NotImplementedError(f'TBExact requires {name} = 1')

        r = {k: self._rate_dt(p[k]) for k in ['inf_cle', 'inf_non', 'inf_asy', 'non_rec', 'non_asy', 'asy_non', 'asy_sym', 'sym_asy', 'sym_dead']}
        Q = np.zeros((len(_chain), len(_chain)))
        Q[0, 4], Q[0, 1], Q[0, 2] = r['inf_cle'], r['inf_non'], r['inf_asy']
        Q[1, 5], Q[1, 2] = r['non_rec'], r['non_asy']
        Q[2, 1], Q[2, 3] = r['asy_non'], r['asy_sym']
        Q[3, 2], Q[3, 6] = r['sym_asy'], r['sym_dead']
        Q[np.diag_indices_from(Q)] = -Q.sum(axis=1)
        cumP = np.cumsum(expm(Q), axis=1)

        rngs = [self._rng_inf, self._rng_non, self._rng_asy, self._rng_sym]
        start_masks = [self.latent, self.non_infectious, self.asymptomatic, self.symptomatic]
        moves = []
        for i, (mask, rng) in enumerate(zip(start_masks, rngs)):
            uids = mask.uids
            if len(uids):
                u = rng.rvs(uids)
                dest = np.minimum((u[:, None] >= cumP[i][None, :]).sum(axis=1), len(_chain) - 1)
                moves.append((uids, dest))

        # Apply after sampling everyone, so each agent moves exactly once from its start-of-step state
        for uids, dest in moves:
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
