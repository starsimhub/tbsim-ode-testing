"""
Experiment 01: simplest matched comparison of tbsim.TB (ABM) and TBODE. See README.md.
"""

import os
import numpy as np
import pandas as pd
import sciris as sc
import matplotlib
matplotlib.use('Agg')
import tbode

outdir = os.path.dirname(os.path.abspath(__file__))
start, stop, dt = 1900, 2100, 1/12
betas = [9, 12, 20]
n_seeds = 10
n_agents = 20_000
n_seeds_big = 5
n_agents_big = 100_000
report_years = [1901, 1910, 1950, 2000, 2050, 2099]


def run_abm(beta, seed, n_agents):
    """Run one ABM replicate and return its comparison outputs."""
    ode = tbode.TBODE(beta=beta, replace_tb_deaths=False)
    sim = tbode.make_abm(ode.pars, seed=seed, n_agents=n_agents, start=start, stop=stop, dt=dt)
    sim.run()
    return tbode.abm_outputs(sim)


def run_ode(beta, replace_tb_deaths):
    return tbode.ode_outputs(tbode.TBODE(beta=beta, replace_tb_deaths=replace_tb_deaths).run(start=start, stop=stop, dt=dt))


if __name__ == '__main__':
    T = sc.timer()

    # Runs: (label, beta, n_agents, n_seeds)
    configs = [(f'beta={b}', b, n_agents, n_seeds) for b in betas] + [(f'beta=9, N=100k', 9, n_agents_big, n_seeds_big)]
    iterkwargs = [dict(beta=b, seed=s, n_agents=n) for _, b, n, ns in configs for s in range(ns)]
    outs = sc.parallelize(run_abm, iterkwargs=iterkwargs, ncpus=tbode.ncpus)
    T.toc('ABM runs')

    # Group replicates by config
    abms, i = {}, 0
    for lab, b, n, ns in configs:
        abms[lab] = dict(beta=b, n_agents=n, outs=outs[i:i + ns])
        i += ns

    odes = {b: {rep: run_ode(b, rep) for rep in [False, True]} for b in betas}

    # Figures: one per beta
    for b in betas:
        ab = {lab: tbode.summarize(a['outs']) for lab, a in abms.items() if a['beta'] == b}
        od = {'births=muN (matched)': odes[b][False], 'TB deaths replaced (TB_ODE)': odes[b][True]}
        fig = tbode.plot_comparison(od, ab, title=f'Baseline, ODE beta={b}')
        fig.savefig(os.path.join(outdir, f'baseline_beta{b}.png'), dpi=110)
        fig = tbode.plot_comparison(od, ab, xlim=(start, start + 20), title=f'Baseline, ODE beta={b}: first 20 years')
        fig.savefig(os.path.join(outdir, f'baseline_beta{b}_early.png'), dpi=110)

    # Summary table: ABM median (10-90%) vs ODE at selected years
    rows = []
    for lab, a in abms.items():
        summ = tbode.summarize(a['outs'])
        o = odes[a['beta']][False]
        o_rep = odes[a['beta']][True]
        n_extinct = sum(out.prev_active[-1] == 0 for out in a['outs'])
        for y in report_years:
            ia = np.argmin(np.abs(summ.time - y))
            io = np.argmin(np.abs(o.time - y))
            for k in ['prev_active', 'INFECTION', 'CLEARED', 'SUSCEPTIBLE', 'tb_deaths_ppy']:
                rows.append(dict(config=lab, year=y, output=k, abm_p10=summ[k][0, ia], abm_median=summ[k][1, ia], abm_p90=summ[k][2, ia],
                                 ode=o[k][io], ode_tb_deaths_replaced=o_rep[k][io], ratio_abm_ode=summ[k][1, ia] / o[k][io],
                                 n_extinct=n_extinct, n_seeds=len(a['outs'])))
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(outdir, 'summary.csv'), index=False)
    with pd.option_context('display.width', 250, 'display.max_rows', 500, 'display.float_format', '{:.4g}'.format):
        print(df[df.output.isin(['prev_active', 'INFECTION', 'CLEARED'])].to_string(index=False))
    T.toc('Total')
