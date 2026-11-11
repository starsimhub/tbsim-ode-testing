"""
Exp 23: (1) the pitfall checks' few "violations" -- which transitions did those agents make in that step? (2) agents ÷ ODE per prevalent case.
    NUMBA_NUM_THREADS=1 PYTHONPATH=v4 v4/.venv/bin/python experiments/23_treatment/diagnose.py
"""
import os
os.environ.setdefault('NUMBA_NUM_THREADS', '1')
import sys
import importlib.util
import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss

outdir = sc.path(os.environ.get('TB_EXP23_FOLDER', sc.thispath(__file__))) # TB_EXP23_FOLDER redirects the outputs (exp 26 reruns this at newer pins)
here = sc.thispath(__file__)
spec = importlib.util.spec_from_file_location('exp23', here / 'run.py')
r = importlib.util.module_from_spec(spec)
spec.loader.exec_module(r)
import tbode4 # noqa: E402

names = ['tx_start', 'tx_cure', 'tx_fail', 'relapse', 'tb_death', 'tb_death_tx', 'deaths', 'symptoms', 'remission']


class Trace(tbode4.TxChecks):
    """ TxChecks, also listing every transition made in the step by each agent flagged by a check """
    def init_post(self):
        super().init_post()
        self.cases = []

    def step(self):
        tb = self.sim.diseases.tb
        get = lambda name: tbode4.treat._uids(tb, name)
        starts, fails = get('tx_start'), get('tx_fail')
        self._grow(np.concatenate([starts, get('tb_death'), get('tb_death_tx'), [0]]).astype(int))
        flagged = {('start_while_on_tx', int(u)) for u in np.setdiff1d(starts[self.prev_tx[starts]], fails)}
        flagged |= {('tb_death_from_success', int(u)) for u in np.concatenate([get('tb_death'), get('tb_death_tx')]).astype(int) if self.prev_succ[u]}
        for check, u in flagged:
            self.cases.append(dict(check=check, uid=u, ti=self.ti, made=[n for n in names if u in set(get(n).tolist())]))
        super().step()


def trace(scen, course, dt, rep):
    sim, _ = r.make_sim(scen, course, dt, 'abm', rand_seed=r.seed(scen, course, dt, rep))
    sim.pars.analyzers = [sim.pars.analyzers[0], Trace()] if isinstance(sim.pars.analyzers, list) else [Trace()]
    sim.run()
    return [dict(c, scen=scen, course=course, dt=dt, rep=rep) for c in sim.analyzers[1].cases]


if __name__ == '__main__':
    # (1) Rerun the groups with flagged cases (same seeds, so the same runs) and list each case's transitions in its step
    groups = [('eq', 'month', 'rate'), ('eq', 'month', 'fixed'), ('eq', 'week', 'rate'), ('ramp', 'month', 'rate'), ('ramp', 'month', 'fixed')]
    jobs = [dict(scen=s, dt=d, course=c, rep=i) for s, d, c in groups for i in range(r.reps)]
    cases = sum(sc.parallelize(trace, iterkwargs=jobs, ncpus=tbode4.ncpus), [])
    df = pd.DataFrame(cases)
    df['made'] = df['made'].apply(lambda m: ' + '.join(m))
    df.to_csv(outdir / 'diagnose_cases.csv', index=False)
    print('Flagged cases and the transitions each agent made in that step:')
    print(df.to_string(index=False))

    # (2) Agents ÷ ODE per prevalent case: each output ÷ prevalence, by rep, then ÷ the ODE's
    runs, ode = sc.load(r.fn.runs), sc.load(r.fn.ode)
    rows = []
    for scen in r.scenarios:
        for dt in r.dts:
            for course in r.courses:
                o = r.window_mean(next(x for x in ode.runs if x.scen == scen and x.dt == dt and x.arm == course).yearly, scen)
                vals = []
                for x in runs:
                    if x.scen == scen and x.dt == dt and x.course == course:
                        a = r.window_mean(x.yearly, scen)
                        vals.append({k: (a[k]/a.prevalence)/(o[k]/o.prevalence) for k in ['incidence', 'notifications', 'on_treatment', 'mortality', 'mortality_tx', 'relapse']})
                v = pd.DataFrame(vals)
                row = dict(scenario=scen, dt=dt, course=course)
                for k in v:
                    row[k] = v[k].mean(); row[k + '_se'] = v[k].std()/np.sqrt(len(v))
                rows.append(row)
    per = pd.DataFrame(rows)
    per.to_csv(outdir / 'per_prevalent_case.csv', index=False)
    with pd.option_context('display.width', 250, 'display.precision', 4):
        print('\nAgents ÷ ODE, each output per prevalent case (mean ± SE over 40 reps):')
        print(per.to_string(index=False))
    # Pooled prevalence ratio by scenario
    avo = pd.read_csv(outdir / 'agents_vs_ode.csv')
    for scen, d in avo.groupby('scenario'):
        print(f'{scen}: prevalence agents ÷ ODE pooled over the 4 arms {d.prevalence.mean():.4f} ± {np.sqrt((d.prevalence_se**2).sum())/len(d):.4f}')
    o = next(x for x in ode.runs if x.scen == 'ramp' and x.dt == 'month' and x.arm == 'fixed').yearly
    print('\nODE, scale-up, fixed course, per 100,000:'); print(o.loc[[1980, 1994, 2000, 2010, 2015, 2020, 2024, 2035], ['prevalence', 'incidence', 'notifications', 'mortality']].round(1).to_string())
