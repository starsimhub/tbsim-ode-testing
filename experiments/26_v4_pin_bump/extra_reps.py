"""
Exp 26: 60 extra reps (reps 40–99, new seeds) of the equilibrium scenario with a fixed course at weekly steps, the arm whose prevalence
read low in exp 23 (0.962 ± 0.017) and in exp 26's first 40 reps (0.952 ± 0.016). Uses exp 23's run.py unchanged. Run from the repo root:

    PYTHONPATH=v4 v4/.venv/bin/python experiments/26_v4_pin_bump/extra_reps.py
"""
import os
os.environ.setdefault('NUMBA_NUM_THREADS', '1')
os.environ.setdefault('TB_EXP23_FOLDER', os.path.dirname(os.path.abspath(__file__)))
import importlib.util
import numpy as np
import pandas as pd
import sciris as sc

here = sc.thispath(__file__)
spec = importlib.util.spec_from_file_location('exp23', here.parent / '23_treatment' / 'run.py')
r = importlib.util.module_from_spec(spec)
spec.loader.exec_module(r)
import tbode4 # noqa: E402

scen, course, dt, new_reps = 'eq', 'fixed', 'week', range(40, 100)

if __name__ == '__main__':
    fn = here / 'extra_reps.obj'
    if not fn.exists():
        extra = sc.parallelize(r.run_agents, iterkwargs=[dict(scen=scen, course=course, dt=dt, rep=k) for k in new_reps], ncpus=tbode4.ncpus)
        sc.save(fn, extra)
    extra = sc.load(fn)
    first = [x for x in sc.load(r.fn.runs) if x.scen == scen and x.course == course and x.dt == dt]
    ode = r.window_mean(next(y for y in sc.load(r.fn.ode).runs if y.scen == scen and y.dt == dt and y.arm == course).yearly, scen)
    rows = []
    for label, runs in [('reps 0-39', first), ('reps 40-99', extra), ('all 100', first + extra)]:
        v = np.array([r.window_mean(x.yearly, scen).prevalence/ode.prevalence for x in runs])
        rows.append(dict(reps=label, n=len(v), prevalence=v.mean(), se=v.std(ddof=1)/np.sqrt(len(v)), p10=np.percentile(v, 10), p90=np.percentile(v, 90)))
    df = pd.DataFrame(rows)
    df.to_csv(here / 'extra_reps.csv', index=False)
    print(f'Prevalence agents ÷ ODE, {scen}, {course} course, {dt}ly steps, {r.scenarios[scen].window} (mean ± SE over reps):')
    print(df.to_string(index=False))
