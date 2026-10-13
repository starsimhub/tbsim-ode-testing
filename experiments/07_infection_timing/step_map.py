"""
Experiment 07 diagnostic: deterministic expected-value version of the TBExact step (start-of-step force of infection, then exact
transitions, with 'start' or 'within' infection timing), iterated to equilibrium at daily, weekly and monthly dt. Output: step_map.log.
"""
import numpy as np, tbode
from scipy.linalg import expm
p = tbode.TBODE().pars
beta = 9.917337993130424
# states: 0 S, 1 INF, 2 CLE(inf, rr 1.0), 3 REC(non, rr 0.21), 4 NON, 5 ASY, 6 SYM, 7 DEAD(TB)
Q = np.zeros((8, 8))
Q[1, 2], Q[1, 4], Q[1, 5] = p.inf_cle, p.inf_non, p.inf_asy
Q[4, 3], Q[4, 5] = p.non_rec, p.non_asy
Q[5, 4], Q[5, 6] = p.asy_non, p.asy_sym
Q[6, 5], Q[6, 7] = p.sym_asy, p.sym_dead
Q[np.diag_indices_from(Q)] = -Q.sum(1)
A = np.zeros((16, 16))

def step(x, dt, timing, P, Pavg):
    N = x[:7].sum()
    w = (p.trans_asymp * x[5] + x[6]) / N
    h = beta * w * dt  # expected infections per susceptible-equivalent (beta_edge_from_ode: K*q = beta*dt)
    newS, newC, newR = x[0] * h, x[2] * h * p.rr_reinfection_cleared, x[3] * h * p.rr_reinfection_rec
    new = newS + newC + newR
    y = x.copy(); y[0] -= newS; y[2] -= newC; y[3] -= newR
    old = y.copy()
    out = old @ P  # all existing agents through the full step
    e1 = np.zeros(8); e1[1] = new
    out += e1 @ (P if timing == 'start' else Pavg)
    # background deaths and births (mu), births replace background deaths only (TB deaths not replaced)
    alive = out[:7]
    deaths = alive * (1 - np.exp(-p.mu * dt))
    alive = alive - deaths
    alive[0] += x[:7].sum() * (1 - np.exp(-p.mu * dt))  # births at rate mu*N (start-of-step N)
    out[:7] = alive; out[7] = 0
    return out / out[:7].sum()

ode = tbode.TBODE(beta=beta, replace_tb_deaths=False, init_sym=0.0)
x0 = ode.equilibrium()
xo = np.array([x0[0], x0[1], x0[2], x0[3], x0[4], x0[5], x0[6], 0])
print('ODE prev', (xo[5] + xo[6]) * 1e5)
for dtl, dt in [('day', 1/365), ('week', 7/365), ('month', 1/12)]:
    P = expm(Q * dt)
    A[:8, :8] = Q * dt; A[:8, 8:] = np.eye(8)
    Pavg = expm(A)[:8, 8:]
    for timing in ['start', 'within']:
        x = xo.copy()
        for i in range(int(3000 / dt)):
            x = step(x, dt, timing, P, Pavg)
        # finish convergence check
        x2 = x.copy()
        for i in range(int(200 / dt)):
            x2 = step(x2, dt, timing, P, Pavg)
        print(f'{dtl:5s} {timing:6s}: prev ratio {(x2[5]+x2[6])/(xo[5]+xo[6]):.4f}  (drift over last 200y {((x2[5]+x2[6])/(x[5]+x[6])-1):.1e})  INF ratio {x2[1]/xo[1]:.4f}  CLE ratio {(x2[2]+x2[3])/(xo[2]+xo[3]):.4f}')
