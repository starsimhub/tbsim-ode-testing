"""
Experiment 10: deterministic expected-value version of the TBExact step, iterated to equilibrium, for every arm.

The step: infections from the start-of-step state, with per-step probability x = beta*w*rel_sus*dt ('linear') or 1 - exp(-x)
('exp'); then exact transitions (latent classes when k > 0), with new infections stepped through the full step ('start')
or the rest of it after a uniform infection time ('within'); then background deaths replaced by births. This is what TBExact
would give without sampling noise or finite-population effects. Output: step_map.log.
"""

import os
import json
import numpy as np
from scipy.linalg import expm
import tbode

outdir = os.path.dirname(os.path.abspath(__file__))
dir03 = os.path.join(outdir, '..', '03_baseline_calibrated')
dir08 = os.path.join(outdir, '..', '08_progression_tsi')


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


def equilibrium(ode, dt, infection_prob='linear', timing='within', years=3000):
    """Iterate the expected step to equilibrium; return prevalence and INFECTION as ratios to the ODE equilibrium."""
    p = ode.pars
    nL = ode.n_latent
    Q = generator(ode)
    n = len(Q)
    P = expm(Q * dt)
    A = np.zeros((2 * n, 2 * n))
    A[:n, :n] = Q * dt
    A[:n, n:] = np.eye(n)
    P_new = expm(A)[:n, n:] if timing == 'within' else P  # 'within': averaged over a uniform infection time
    xo = ode.equilibrium()
    x = np.concatenate([[xo[0]], xo[1:1 + nL], [xo[ode.idx('CLEARED')], xo[ode.idx('RECOVERED')], xo[ode.idx('NON_INFECTIOUS')],
                                                xo[ode.idx('ASYMPTOMATIC')], xo[ode.idx('SYMPTOMATIC')], 0]])
    iC, iR, iA, iS = 1 + nL, 2 + nL, 4 + nL, 5 + nL
    prev_ode, inf_ode = xo[ode.idx('ASYMPTOMATIC')] + xo[ode.idx('SYMPTOMATIC')], xo[1:1 + nL].sum()
    prob = (lambda v: 1 - np.exp(-v)) if infection_prob == 'exp' else (lambda v: np.minimum(v, 1))
    surv = np.exp(-p.mu * dt)
    for _ in range(int(years / dt)):
        N = x[:-1].sum()
        xx = p.beta * (p.trans_asymp * x[iA] + x[iS]) / N * dt
        newS, newC, newR = x[0] * prob(xx), x[iC] * prob(xx * p.rr_reinfection_cleared), x[iR] * prob(xx * p.rr_reinfection_rec)
        chain = x[1:].copy()
        chain[iC - 1] -= newC
        chain[iR - 1] -= newR
        new = np.zeros(n)
        new[:nL] = (newS + newC + newR) * ode.latent_weights
        out = chain @ P + new @ P_new
        S = (x[0] - newS) * surv + N * (1 - surv)  # Background deaths, replaced by births (TB deaths not replaced)
        out[:-1] *= surv
        out[-1] = 0
        x = np.concatenate([[S], out])
        x /= x[:-1].sum()
    return (x[iA] + x[iS]) / prev_ode, x[1:1 + nL].sum() / inf_ode


def settings():
    """ODE parameter sets: 'k0' (tbsim defaults, baseline beta, exps 03-07) and sets A and B (exp 08)."""
    calib = json.load(open(os.path.join(dir08, 'calibration.json')))
    beta0 = json.load(open(os.path.join(dir03, 'calibration.json')))['beta']
    out = dict(k0=dict(beta=beta0, replace_tb_deaths=False, init_sym=0.0))
    out.update({s: dict(beta=c['beta'], **c['pars']) for s, c in calib.items()})
    return out


if __name__ == '__main__':
    for s, pars in settings().items():
        ode = tbode.TBODE(**pars)
        for dtl, dt in [('1 week', 7/365), ('1 month', 1/12)]:
            for timing in (['start', 'within'] if s == 'k0' else ['within']):
                vals = {form: equilibrium(ode, dt, form, timing)[0] for form in ['linear', 'exp']}
                print(f"{s:3s} {dtl:7s} {timing:6s}: prevalence ÷ ODE  linear {vals['linear']:.4f}  exp {vals['exp']:.4f}")
