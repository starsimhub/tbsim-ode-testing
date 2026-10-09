"""
Experiment 05: experiment 03 with tbode.TBExact (exact per-step natural-history sampling). See README.md.
"""

import os
import json
import numpy as np
import pandas as pd
import sciris as sc
import matplotlib
matplotlib.use('Agg')
import tbode

outdir = os.path.dirname(os.path.abspath(__file__))
dir03 = os.path.join(outdir, '..', '03_baseline_calibrated')
beta = json.load(open(os.path.join(dir03, 'calibration.json')))['beta']
start, stop = 1900, 2000
n_agents, n_seeds = 50_000, 40
dts = {'1 day': 1/365, '1 week': 7/365, '1 month': 1/12} # Daily first, so the longest runs start first
ode_kw = dict(replace_tb_deaths=False, init_sym=0.0)
outputs = ['prev_active', 'incidence_py', 'tb_deaths_py', 'frac_asymp']
scale = dict(prev_active=1e5, incidence_py=1e5, tb_deaths_py=1e5, frac_asymp=1)


def run_abm(x, dt, seed):
    ode = tbode.TBODE(beta=beta, **ode_kw)
    sim = tbode.make_abm(ode.pars, seed=seed, n_agents=n_agents, start=start, stop=stop, dt=dt, tb_class=tbode.TBExact)
    sim.init()
    tbode.init_from_ode(sim, x, seed=seed)
    sim.run()
    return tbode.annualize(tbode.abm_outputs(sim))


if __name__ == '__main__':
    T = sc.timer()
    ode = tbode.TBODE(beta=beta, **ode_kw)
    x = ode.equilibrium()
    eq = ode.outputs(x)
    iterkwargs = [dict(x=x, dt=dt, seed=s) for dt in dts.values() for s in range(n_seeds)]
    outs = sc.parallelize(run_abm, iterkwargs=iterkwargs, ncpus=tbode.ncpus)
    abms = {lab: outs[i * n_seeds:(i + 1) * n_seeds] for i, lab in enumerate(dts)}
    ode_run = tbode.ode_outputs(tbode.TBODE(beta=beta, **ode_kw).run(start=start, stop=stop, dt=1/52, y0=x * 1e5))
    sc.save(os.path.join(outdir, 'outputs.obj'), dict(abms=abms, ode_run=ode_run))
    T.toc('ABM runs')

    # Mean +/- SE (and median, 10-90%) across reps of each rep's last-50-year mean, ABM / ODE; TBExact here and tbsim.TB from exp 03
    tb03 = pd.read_csv(os.path.join(dir03, 'summary_more_reps.csv'))
    rows = []
    for k in outputs:
        ode_val = eq[k][0] * scale[k]
        for lab, oo in abms.items():
            vals = np.array([np.nanmean(o[k][o.time >= stop - 50]) for o in oo]) * scale[k]
            se = vals.std(ddof=1) / np.sqrt(len(vals))
            rows.append(dict(output=k, model='TBExact', dt=lab, n_reps=len(vals), ode=ode_val, mean=vals.mean(), se=se, ratio_mean=vals.mean() / ode_val,
                             ratio_se=se / ode_val, z=(vals.mean() - ode_val) / se, median=np.median(vals), p10=np.quantile(vals, 0.1), p90=np.quantile(vals, 0.9)))
        for _, r in tb03[tb03.output == k].iterrows():
            rows.append(dict(r, model='tbsim.TB (exp 03)'))
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(outdir, 'summary.csv'), index=False)

    summ = {f'TBExact dt={lab}': tbode.summarize(oo) for lab, oo in abms.items()}
    keys = ['INFECTION', 'CLEARED', 'NON_INFECTIOUS', 'ASYMPTOMATIC', 'prev_active', 'incidence_py', 'tb_deaths_py', 'frac_asymp']
    fig = tbode.plot_comparison({'equilibrium': ode_run}, summ, keys=keys, title=f'TBExact at equilibrium, ODE beta={beta:.3f} (500/100k), {n_seeds} reps per dt',
                                scale=dict(prev_active=1e5, incidence_py=1e5, tb_deaths_py=1e5))
    fig.savefig(os.path.join(outdir, 'exact_equilibrium.png'), dpi=110)

    with pd.option_context('display.width', 250, 'display.float_format', '{:.4g}'.format):
        print(df[['output', 'model', 'dt', 'n_reps', 'ode', 'mean', 'se', 'ratio_mean', 'ratio_se', 'z', 'median', 'p10', 'p90']].to_string(index=False))
    T.toc('Total')
