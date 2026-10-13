"""
Experiment 08 diagnostic: deterministic expected-value version of the TBExact step with latent classes (start-of-step force of
infection, then exact transitions, 'within' infection timing), iterated to equilibrium for sets A and B at weekly and monthly dt.
This is what TBExact would give without sampling noise or finite-population effects. Output: step_map.log.
"""

import json
import os
import numpy as np
from scipy.linalg import expm
import tbode

outdir = os.path.dirname(os.path.abspath(__file__))


def generator(ode):
    """Per-year generator over [latent classes..., CLE (from latent), REC (from NON), NON, ASY, SYM, DEAD]."""
    p = ode.pars
    nL = ode.n_latent
    iC, iR, iN, iA, iS, iD = range(nL, nL + 6)
    Q = np.zeros((nL + 6, nL + 6))
    for m in range(nL):
        Q[m, iC], Q[m, iN], Q[m, iA] = p.inf_cle, ode.rate_non[m], ode.rate_asy[m]
    Q[iN, iR], Q[iN, iA] = p.non_rec, p.non_asy
    Q[iA, iN], Q[iA, iS] = p.asy_non, p.asy_sym
    Q[iS, iA], Q[iS, iD] = p.sym_asy, p.sym_dead
    Q[np.diag_indices_from(Q)] = -Q.sum(axis=1)
    return Q


def equilibrium(ode, dt, years=3000):
    """Iterate the expected step to equilibrium; return prevalence and INFECTION as ratios to the ODE equilibrium."""
    p = ode.pars
    nL = ode.n_latent
    Q = generator(ode)
    n = len(Q)
    P = expm(Q * dt)
    A = np.zeros((2 * n, 2 * n))
    A[:n, :n] = Q * dt
    A[:n, n:] = np.eye(n)
    Pavg = expm(A)[:n, n:]  # Averaged over a uniform infection time within the step
    xo = ode.equilibrium()
    # ODE layout: S, latent..., CLE, REC, NON, ASY, SYM, TX, TXD, NON0 -> here: S + chain
    x = np.concatenate([[xo[0]], xo[1:1 + nL], [xo[ode.idx('CLEARED')], xo[ode.idx('RECOVERED')], xo[ode.idx('NON_INFECTIOUS')],
                                                xo[ode.idx('ASYMPTOMATIC')], xo[ode.idx('SYMPTOMATIC')], 0]])
    iC, iR, iA, iS = 1 + nL, 2 + nL, 4 + nL, 5 + nL
    prev_ode, inf_ode = xo[ode.idx('ASYMPTOMATIC')] + xo[ode.idx('SYMPTOMATIC')], xo[1:1 + nL].sum()
    surv = np.exp(-p.mu * dt)
    for _ in range(int(years / dt)):
        N = x[:-1].sum()
        h = p.beta * (p.trans_asymp * x[iA] + x[iS]) / N * dt  # Expected infections per susceptible-equivalent (beta_edge_from_ode)
        newS, newC, newR = x[0] * h, x[iC] * h * p.rr_reinfection_cleared, x[iR] * h * p.rr_reinfection_rec
        chain = x[1:].copy()
        chain[iC - 1] -= newC
        chain[iR - 1] -= newR
        new = np.zeros(n)
        new[:nL] = (newS + newC + newR) * ode.latent_weights
        out = chain @ P + new @ Pavg
        S = (x[0] - newS) * surv + N * (1 - surv)  # Background deaths, replaced by births (TB deaths not replaced)
        out[:-1] *= surv
        out[-1] = 0
        x = np.concatenate([[S], out])
        x /= x[:-1].sum()
    return (x[iA] + x[iS]) / prev_ode, x[1:1 + nL].sum() / inf_ode


if __name__ == '__main__':
    calib = json.load(open(os.path.join(outdir, 'calibration.json')))
    for s, c in calib.items():
        ode = tbode.TBODE(beta=c['beta'], **c['pars'])
        for lab, dt in [('1 week', 7/365), ('1 month', 1/12)]:
            prev, inf = equilibrium(ode, dt)
            print(f'Set {s}, {lab:7s}: deterministic TBExact step ÷ ODE: prevalence {prev:.4f}, INFECTION {inf:.4f}')
