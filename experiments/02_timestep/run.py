"""
Experiment 02: does the ABM-ODE baseline gap shrink with the ABM timestep? See README.md.
"""

import os
import numpy as np
import pandas as pd
import sciris as sc
import matplotlib
matplotlib.use('Agg')
import tbode

outdir = os.path.dirname(os.path.abspath(__file__))
start, stop = 1900, 1960
betas = [12, 20]
dts = {'1 month': (1/12, 10), '1 week': (7/365, 10), '1 day': (1/365, 5)}  # label: (dt in years, n_seeds)
n_agents = 20_000
report_years = [1901, 1910, 1950, 1959]


def run_abm(beta, dt, seed):
    ode = tbode.TBODE(beta=beta, replace_tb_deaths=False)
    sim = tbode.make_abm(ode.pars, seed=seed, n_agents=n_agents, start=start, stop=stop, dt=dt)
    sim.run()
    return tbode.abm_outputs(sim)


if __name__ == '__main__':
    T = sc.timer()
    configs = [(b, lab, dt, ns) for b in betas for lab, (dt, ns) in dts.items()]
    iterkwargs = [dict(beta=b, dt=dt, seed=s) for b, _, dt, ns in configs for s in range(ns)]
    outs = sc.parallelize(run_abm, iterkwargs=iterkwargs, ncpus=tbode.ncpus)
    T.toc('ABM runs')

    rows, i = [], 0
    summs = {b: {} for b in betas}
    for b, lab, dt, ns in configs:
        summ = tbode.summarize(outs[i:i + ns])
        i += ns
        summs[b][f'dt={lab}'] = summ
        o = tbode.ode_outputs(tbode.TBODE(beta=b, replace_tb_deaths=False).run(start=start, stop=stop, dt=1/52))
        for y in report_years:
            ia = np.argmin(np.abs(summ.time - y))
            io = np.argmin(np.abs(o.time - y))
            for k in ['prev_active', 'INFECTION', 'CLEARED', 'ASYMPTOMATIC', 'SYMPTOMATIC', 'tb_deaths_ppy']:
                rows.append(dict(beta=b, dt=lab, year=y, output=k, abm_median=summ[k][1, ia], ode=o[k][io],
                                 ratio_abm_ode=summ[k][1, ia] / o[k][io], n_seeds=ns))

    for b in betas:
        o = tbode.ode_outputs(tbode.TBODE(beta=b, replace_tb_deaths=False).run(start=start, stop=stop, dt=1/52))
        fig = tbode.plot_comparison({'births=muN': o}, summs[b], title=f'ABM timestep vs ODE, beta={b}')
        fig.savefig(os.path.join(outdir, f'timestep_beta{b}.png'), dpi=110)

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(outdir, 'summary.csv'), index=False)
    piv = df[df.year.isin([1901, 1959])].pivot_table(index=['beta', 'year', 'output'], columns='dt', values='ratio_abm_ode')
    with pd.option_context('display.width', 250, 'display.float_format', '{:.3f}'.format):
        print(piv[list(dts)])
    T.toc('Total')
