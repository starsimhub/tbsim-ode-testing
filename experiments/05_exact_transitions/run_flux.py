"""
Experiment 05, flux check: per-step natural-history transitions of tbode.TBExact vs the exact matrix expm(Q*dt),
alongside standard tbsim.TB (same seeds). Same setup as experiment 04. See README.md.
"""

import os
import json
import numpy as np
import pandas as pd
import sciris as sc
import tbode

outdir = os.path.dirname(os.path.abspath(__file__))
beta = json.load(open(os.path.join(outdir, '..', '03_baseline_calibrated', 'calibration.json')))['beta']
n_seeds = 20
dts = {'1 week': 7/365, '1 month': 1/12}
classes = {'tbsim.TB': None, 'TBExact': tbode.TBExact}


if __name__ == '__main__':
    T = sc.timer()
    configs = [(cl, lab, dt) for cl in classes for lab, dt in dts.items()]
    iterkwargs = [dict(beta=beta, dt=dt, seed=s, tb_class=classes[cl]) for cl, _, dt in configs for s in range(n_seeds)]
    outs = sc.parallelize(tbode.run_flux, iterkwargs=iterkwargs, ncpus=tbode.ncpus)
    T.toc('ABM runs')

    pars = tbode.TBODE(beta=beta, replace_tb_deaths=False).pars
    nh, tr = [], []
    for i, (cl, lab, dt) in enumerate(configs):
        df, t = tbode.flux_tables(outs[i * n_seeds:(i + 1) * n_seeds], pars, dt, lab)
        df.insert(0, 'model', cl)
        t['model'] = cl
        nh.append(df)
        tr.append(t)
    nh = pd.concat(nh, ignore_index=True)
    nh['z'] = (nh.ratio - 1) / nh.se_ratio
    tr = pd.DataFrame(tr)
    nh.to_csv(os.path.join(outdir, 'flux_natural_history.csv'), index=False)
    tr.to_csv(os.path.join(outdir, 'flux_transmission.csv'), index=False)

    # Main flows (exclude 'stays', which are ~1 by construction, and the tiny background-death flows)
    main = nh[~nh.flow.str.endswith('stays') & ~nh.flow.str.endswith('->DEAD') | (nh.flow == 'SYMPTOMATIC->DEAD')]
    with pd.option_context('display.width', 250, 'display.max_rows', 200, 'display.float_format', '{:.3f}'.format):
        print(main.pivot_table(index='flow', columns=['model', 'dt'], values='ratio', sort=False))
        print(main.pivot_table(index='flow', columns=['model', 'dt'], values='z', sort=False))
        for (cl, lab), g in main.groupby(['model', 'dt'], sort=False):
            chi2 = np.nansum(g.z ** 2)
            print(f'{cl:9s} {lab:8s}: sum z^2 = {chi2:7.1f} over {g.z.notna().sum()} flows; max |z| = {np.nanmax(np.abs(g.z)):.1f}')
        print(tr[['model', 'dt', 'ratio_S', 'se_S', 'ratio_C', 'se_C']].to_string(index=False))
    T.toc('Total')
