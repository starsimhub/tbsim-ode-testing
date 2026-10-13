"""Diagnose exp 07's weekly 'within' shortfall: 'u0' = within's RNG use but start's physics; 'np' = within physics with an independent numpy RNG."""
import sys, os, numpy as np, sciris as sc
from scipy.linalg import expm
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tbode
from run import beta, ode_kw, last50, start, stop, n_agents

class TBExactU0(tbode.TBExact):
    def _new_inf_cum(self, uids, cls, Qs, cum):
        is_new = np.isin(uids, self._new_inf)
        if not is_new.any():
            return cum
        new = uids[is_new]
        self._rng_tinf.rvs(new)  # Draw as 'within' does, then ignore: full step
        return cum

class TBExactNP(tbode.TBExact):
    def _new_inf_cum(self, uids, cls, Qs, cum):
        is_new = np.isin(uids, self._new_inf)
        if not is_new.any():
            return cum
        if not hasattr(self, '_nprng'):
            self._nprng = np.random.default_rng(1000 + self.sim.pars.rand_seed)
        remaining = 1 - self._nprng.random(is_new.sum())
        rows = expm(Qs[cls[is_new]] * remaining[:, None, None])[:, 0, :]
        cum = cum.copy(); cum[is_new] = np.cumsum(rows, axis=1)
        return cum

classes = dict(u0=TBExactU0, np=TBExactNP)

def run_abm(x, variant, seed):
    ode = tbode.TBODE(beta=beta, **ode_kw)
    sim = tbode.make_abm(ode.pars, seed=seed, n_agents=n_agents, start=start, stop=stop, dt=7/365, tb_class=classes[variant], infection_timing='within')
    sim.init(); tbode.init_from_ode(sim, x, seed=seed); sim.run()
    return tbode.annualize(tbode.abm_outputs(sim))

if __name__ == '__main__':
    x = tbode.TBODE(beta=beta, **ode_kw).equilibrium()
    seeds = range(80)
    outs = sc.parallelize(run_abm, iterkwargs=[dict(x=x, variant=v, seed=s) for v in classes for s in seeds], ncpus=tbode.ncpus)
    res = {v: outs[i*len(seeds):(i+1)*len(seeds)] for i, v in enumerate(classes)}
    sc.save(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'outputs_variants.obj'), res)
    for v, oo in res.items():
        p = last50(oo, 'prev_active') / 500
        print(f'{v}: prevalence / ODE = {p.mean():.4f} ± {p.std(ddof=1)/np.sqrt(len(p)):.4f} ({len(p)} reps)')
