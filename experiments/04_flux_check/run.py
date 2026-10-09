"""
Experiment 04: per-step flux check of the ABM against the ODE rates. See README.md.
"""

import os
import json
import pandas as pd
import sciris as sc
import tbode

outdir = os.path.dirname(os.path.abspath(__file__))
beta = json.load(open(os.path.join(outdir, '..', '03_baseline_calibrated', 'calibration.json')))['beta']
n_seeds = 4
dts = {'1 day': 1/365, '1 week': 7/365, '1 month': 1/12}


if __name__ == '__main__':
    T = sc.timer()
    iterkwargs = [dict(beta=beta, dt=dt, seed=s) for dt in dts.values() for s in range(n_seeds)]
    outs = sc.parallelize(tbode.run_flux, iterkwargs=iterkwargs, ncpus=tbode.ncpus)
    T.toc('ABM runs')

    pars = tbode.TBODE(beta=beta, replace_tb_deaths=False).pars
    tabs = [tbode.flux_tables(outs[i * n_seeds:(i + 1) * n_seeds], pars, dt, lab) for i, (lab, dt) in enumerate(dts.items())]
    df = pd.concat([t[0] for t in tabs], ignore_index=True)
    tr = pd.DataFrame([t[1] for t in tabs])
    df.to_csv(os.path.join(outdir, 'natural_history_flows.csv'), index=False)
    tr.to_csv(os.path.join(outdir, 'transmission.csv'), index=False)
    with pd.option_context('display.width', 250, 'display.max_rows', 200, 'display.float_format', '{:.4g}'.format):
        print(df.pivot_table(index='flow', columns='dt', values='ratio', sort=False)[list(dts)])
        print(tr.to_string(index=False))
    T.toc('Total')
