import numpy as np, starsim as ss, tbsim, sciris as sc
from tbsim import TBS
class Rec(ss.Analyzer):
    def init_post(self):
        super().init_post(); self.gaps = []; self.durs = []
    def step(self):
        tx = self.sim.interventions[0]
        u = tx.pending_relapse.uids
        if len(u): self.gaps.extend((tx.ti_relapse[u] - tx.ti_treatment_end[u]).tolist())
        on = self.sim.diseases.tb.on_treatment.uids
        if len(on): self.durs.extend((tx.ti_treatment_end[on] - tx.ti_treatment_start[on]).tolist())
tb = tbsim.TB(init_prev=0.3)
prod = tbsim.Tx(efficacy=0.5, adherence=1.0, p_relapse=0.5, dur_treatment=ss.constant(v=ss.days(180).days))
elig = lambda sim: ((sim.diseases.tb.state == TBS.SYMPTOMATIC) & ~sim.diseases.tb.on_treatment).uids
tx = tbsim.TxDelivery(product=prod, eligibility=elig)
sim = tbsim.Sim(diseases=tb, interventions=tx, analyzers=Rec(), n_agents=5000, start=2000, stop=2004, dt=ss.month, verbose=0, rand_seed=1)
sim.run()
a = sim.analyzers[0]
tx = sim.interventions[0]; prod = tx.product; print('dt.days', tx.dt.days)
print('course length in steps (end - start):', np.unique(np.round(a.durs, 3)))
print('relapse gap in steps (ti_relapse - ti_treatment_end):', np.unique(np.round(a.gaps, 4))[:10])
print('dur_relapse dist:', prod.pars.dur_relapse, ' sample rvs:', prod.pars.dur_relapse.rvs(ss.uids([0,1])))
