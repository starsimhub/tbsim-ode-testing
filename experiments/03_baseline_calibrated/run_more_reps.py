"""
Experiment 03, extra reps: 30 more daily-dt reps (seeds 10-39), added to the 10 in outputs.obj. See README.md.
"""

import os
import json
import numpy as np
import pandas as pd
import sciris as sc
import matplotlib
matplotlib.use('Agg')
import tbode
from run import run_abm, ode_kw, outputs, start, stop

outdir = os.path.dirname(os.path.abspath(__file__))
new_seeds = range(10, 40)
lab, dt = '1 day', 1/365
scale = dict(prev_active=1e5, incidence_py=1e5, tb_deaths_py=1e5, frac_asymp=1)


if __name__ == '__main__':
    T = sc.timer()
    beta = json.load(open(os.path.join(outdir, 'calibration.json')))['beta']
    ode = tbode.TBODE(beta=beta, **ode_kw)
    x = ode.equilibrium()
    eq = ode.outputs(x)

    data = sc.load(os.path.join(outdir, 'outputs.obj'))
    iterkwargs = [dict(beta=beta, x=x, dt=dt, seed=s) for s in new_seeds]
    new = sc.parallelize(run_abm, iterkwargs=iterkwargs, ncpus=tbode.ncpus)
    data['abms'][lab] = data['abms'][lab][:10] + new
    sc.save(os.path.join(outdir, 'outputs_more_reps.obj'), data)
    T.toc('ABM runs')

    # Mean +/- SE and median / 10-90% across reps of the last-50-year mean, for every dt (daily now has 40 reps)
    rows = []
    for k in outputs:
        for dlab, oo in data['abms'].items():
            vals = np.array([np.nanmean(o[k][o.time >= stop - 50]) for o in oo]) * scale[k]
            ode_val = eq[k][0] * scale[k]
            se = vals.std(ddof=1) / np.sqrt(len(vals))
            rows.append(dict(output=k, dt=dlab, n_reps=len(vals), ode=ode_val, mean=vals.mean(), se=se, ratio_mean=vals.mean() / ode_val, ratio_se=se / ode_val,
                             z=(vals.mean() - ode_val) / se, median=np.median(vals), p10=np.quantile(vals, 0.1), p90=np.quantile(vals, 0.9)))
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(outdir, 'summary_more_reps.csv'), index=False)

    summ = {f'dt={dlab} ({len(oo)} reps)': tbode.summarize(oo) for dlab, oo in data['abms'].items()}
    keys = ['INFECTION', 'CLEARED', 'NON_INFECTIOUS', 'ASYMPTOMATIC', 'prev_active', 'incidence_py', 'tb_deaths_py', 'frac_asymp']
    fig = tbode.plot_comparison({'equilibrium': data['ode_run']}, summ, keys=keys, title=f'Baseline at equilibrium, ODE beta={beta:.3f} (500/100k), 40 daily reps',
                                scale=dict(prev_active=1e5, incidence_py=1e5, tb_deaths_py=1e5))
    fig.savefig(os.path.join(outdir, 'baseline_equilibrium_more_reps.png'), dpi=110)

    with pd.option_context('display.width', 250, 'display.float_format', '{:.4g}'.format):
        print(df.to_string(index=False))
    T.toc('Total')
