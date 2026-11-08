"""
Experiment 22, diagnosis: why agents on RandomNet(10, dur=0) run 2-3% below the network standalone ODE's prediction in every arm.

Each step, an intervention (after the network update, before transmission) records the edges and each agent's state; an analyzer
(after the step) counts each transmission transition. For each class of target, it compares the infections observed with:

- the edge-level expectation: sum over targets of 1 - prod over infectious contacts of exp(-rr * rel * b * D), with rel = kappa' for
  asymptomatic and 1 for symptomatic contacts (what the conversion assumes agents do);
- the homogeneous ODE's flow over the step: rr * beta * D * (kappa * A + Y) / N per target (the conversion's aim).

Run from the repo root:

    PYTHONPATH=v4 v4/.venv/bin/python experiments/22_network_conversion/diagnose.py
"""
import sys
import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss

outdir = sc.thispath(__file__)
sys.path.insert(0, str(outdir))
import run as r # noqa: E402
import tbode4 # noqa: E402

targets = sc.objdict(transmission=('susceptible', None, 1.0), reinfect_cle=('cleared', 'rr_reinfection_cleared', None),
                     reinfect_rec=('recovered', 'rr_reinfection_rec', None), reinfect_trd=('treated', 'rr_reinfection_treat', None),
                     reinfect_remote=('infection_remote', 'rr_reinfection_inf', None))


class Snapshot(ss.Intervention):
    """ The edges and states transmission is about to use, with the expected infections by class """
    def __init__(self, netp, homp, D, **kwargs):
        super().__init__(**kwargs)
        self.netp, self.homp, self.D = netp, homp, D
        self.rows = []
        self.pending = None
        return

    def step(self):
        sim = self.sim
        tb = sim.diseases.tb
        n = int(np.asarray(sim.people.uid).max()) + 1
        alive = np.asarray(sim.people.alive.uids)
        A, Y = np.zeros(n, bool), np.zeros(n, bool)
        A[np.asarray(tb.asymptomatic.uids)] = True
        Y[np.asarray(tb.symptomatic.uids)] = True
        e = sim.networks[0].edges
        p1, p2 = np.asarray(e.p1), np.asarray(e.p2)
        nA = np.bincount(p1[A[p2]], minlength=n) + np.bincount(p2[A[p1]], minlength=n) # Asymptomatic contacts of each agent
        nY = np.bincount(p1[Y[p2]], minlength=n) + np.bincount(p2[Y[p1]], minlength=n)
        deg = np.bincount(p1, minlength=n) + np.bincount(p2, minlength=n)
        b, D, N = self.netp['beta'], self.D, len(alive)
        lam_h = self.homp['beta']*(self.homp['trans_asymp']*A.sum() + Y.sum())/N # The homogeneous force of infection, per year
        row = dict(ti=sim.ti, N=N, A=int(A.sum()), Y=int(Y.sum()), mean_deg=deg[alive].mean(), var_deg=deg[alive].var(), self_loops=int((p1 == p2).sum()))
        for tname, (state, rrkey, _) in targets.items():
            if not hasattr(tb, state):
                continue
            u = np.asarray(getattr(tb, state).uids)
            rr_n = self.netp[rrkey] if rrkey else 1.0
            rr_h = self.homp[rrkey] if rrkey else 1.0
            haz = rr_n*b*D*(self.netp['trans_asymp']*nA[u] + nY[u])
            row[f'{tname}_n'] = len(u)
            row[f'{tname}_exp_edge'] = -np.expm1(-haz).sum()
            qY = -np.expm1(-rr_n*b*D) # v4's chance per symptomatic edge; an asymptomatic one is trans_asymp times it (rel_trans scales the chance)
            row[f'{tname}_exp_lin'] = (1 - (1 - qY)**nY[u]*(1 - self.netp['trans_asymp']*qY)**nA[u]).sum()
            row[f'{tname}_exp_ode'] = rr_h*lam_h*D*len(u)
            row[f'{tname}_uids'] = u
        self.pending = row
        return


class Count(ss.Analyzer):
    """ The infections by class in the step, matched to its snapshot """
    def step(self):
        snap = self.sim.interventions[0]
        row = snap.pending
        if row is None:
            return
        tb = self.sim.diseases.tb
        for tname in targets:
            if f'{tname}_uids' not in row:
                continue
            u = np.asarray(getattr(tb, tname).uids)
            row[f'{tname}_obs'] = len(u)
            row[f'{tname}_obs_from_snap'] = int(np.isin(u, row.pop(f'{tname}_uids')).sum()) # Those infected who were in the class at the snapshot
        snap.rows.append(row)
        snap.pending = None
        return


def run_diag(i=3, arm='conveq_month', seed=1, stop=2000):
    sets = sc.load(r.exp21.fn.sets)
    conv = sc.load(r.fn.conv)
    s, netp = sets[i], conv[f's{i}_{arm}'].pars
    dkey = r.arms[arm][1]
    sim = r.make_sim(s, netp, dkey, seed)
    sim.pars.stop = stop
    D = r.dts[dkey].years
    sim.pars.interventions = Snapshot(netp, s.pars, D)
    sim.pars.analyzers = [sim.pars.analyzers, Count()]
    sim.run()
    df = pd.DataFrame(sim.interventions[0].rows)
    df['set'], df['arm'], df['seed'] = i, arm, seed
    return df


def report(df):
    d = df[df.ti > 0]
    print(f"Mean degree {d.mean_deg.mean():.4f} (variance {d.var_deg.mean():.3f}); self-loops per step {d.self_loops.mean():.1f}; N {d.N.mean():.0f}")
    rows = []
    for (i, arm), g in d.groupby(['set', 'arm']):
        for tname in targets:
            if f'{tname}_obs' not in g:
                continue
            obs, sn = g[f'{tname}_obs'].sum(), g[f'{tname}_obs_from_snap'].sum()
            ee, eo, el = g[f'{tname}_exp_edge'].sum(), g[f'{tname}_exp_ode'].sum(), g[f'{tname}_exp_lin'].sum()
            if not obs:
                continue
            rows.append(dict(set=i, arm=arm, transition=tname, observed=obs, obs_per_sat=obs/ee, obs_per_lin=obs/el, se=np.sqrt(obs)/ee, sat_per_ode=ee/eo, lin_per_ode=el/eo))
    out = pd.DataFrame(rows)
    with pd.option_context('display.width', 200):
        print(out.round(4).to_string(index=False))
    return out


if __name__ == '__main__':
    kw = [dict(i=i, arm=arm, seed=s) for i in [3, 7] for arm in ['conveq_month'] for s in [1, 2]]
    dfs = sc.parallelize(run_diag, iterkwargs=kw, ncpus=tbode4.ncpus)
    df = pd.concat(dfs)
    sc.save(outdir / 'diagnose.obj', df.drop(columns=[c for c in df.columns if c.endswith('_uids')], errors='ignore'))
    report(df)
