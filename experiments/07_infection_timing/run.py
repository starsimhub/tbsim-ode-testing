"""
Experiment 07: TBExact with infection times drawn within the step ('within'), vs exp 05's 'start' timing. See README.md.
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
dir05 = os.path.join(outdir, '..', '05_exact_transitions')
beta = json.load(open(os.path.join(dir03, 'calibration.json')))['beta']
start, stop = 1900, 2000
n_agents, n_seeds = 50_000, 40
dts = {'1 week': 7/365, '1 month': 1/12} # Weekly first, so the longest runs start first
ode_kw = dict(replace_tb_deaths=False, init_sym=0.0)
outputs = ['INFECTION', 'CLEARED', 'NON_INFECTIOUS', 'ASYMPTOMATIC', 'SYMPTOMATIC', 'prev_active', 'incidence_py', 'tb_deaths_py', 'frac_asymp']
scale = dict(prev_active=1e5, incidence_py=1e5, tb_deaths_py=1e5)


def run_abm(x, dt, seed, infection_timing):
    """As exp 05's run_abm, with the TBExact infection timing as an argument."""
    ode = tbode.TBODE(beta=beta, **ode_kw)
    sim = tbode.make_abm(ode.pars, seed=seed, n_agents=n_agents, start=start, stop=stop, dt=dt, tb_class=tbode.TBExact, infection_timing=infection_timing)
    sim.init()
    tbode.init_from_ode(sim, x, seed=seed)
    sim.run()
    return tbode.annualize(tbode.abm_outputs(sim))


def last50(oo, k):
    """Each rep's mean over the last 50 years, in display units."""
    return np.array([np.nanmean(o[k][o.time >= stop - 50]) for o in oo]) * scale.get(k, 1)


if __name__ == '__main__':
    T = sc.timer()
    ode = tbode.TBODE(beta=beta, **ode_kw)
    x = ode.equilibrium()
    eq_states = tbode.ode_outputs(ode.outputs(x, time=np.zeros(1)))

    # The default 'start' timing must reproduce exp 05 exactly
    out05 = sc.load(os.path.join(dir05, 'outputs.obj'))
    check = run_abm(x, dts['1 month'], 0, 'start')
    ref = out05['abms']['1 month'][0]
    for k in outputs:
        assert np.allclose(check[k], ref[k], equal_nan=True), f"'start' timing doesn't reproduce exp 05 for {k}"
    print("Check passed: 'start' timing reproduces exp 05 (seed 0, monthly)")

    iterkwargs = [dict(x=x, dt=dt, seed=s, infection_timing='within') for dt in dts.values() for s in range(n_seeds)]
    outs = sc.parallelize(run_abm, iterkwargs=iterkwargs, ncpus=tbode.ncpus)
    abms = {lab: outs[i * n_seeds:(i + 1) * n_seeds] for i, lab in enumerate(dts)}
    ode_run = tbode.ode_outputs(tbode.TBODE(beta=beta, **ode_kw).run(start=start, stop=stop, dt=1/52, y0=x * 1e5))
    sc.save(os.path.join(outdir, 'outputs.obj'), dict(abms=abms, ode_run=ode_run))
    T.toc('ABM runs')

    # Mean +/- SE across reps of each rep's last-50-year mean, ABM / ODE, for 'within' (here) and 'start' (exp 05); plus the paired difference by seed
    rows = []
    for k in outputs:
        ode_val = eq_states[k][0] * scale.get(k, 1)
        for lab in dts:
            vals = {'within': last50(abms[lab], k), 'start': last50(out05['abms'][lab], k)}
            for timing, v in vals.items():
                se = v.std(ddof=1) / np.sqrt(len(v))
                rows.append(dict(output=k, dt=lab, timing=timing, n_reps=len(v), ode=ode_val, mean=v.mean(), se=se, ratio_mean=v.mean() / ode_val,
                                 ratio_se=se / ode_val, z=(v.mean() - ode_val) / se, median=np.median(v), p10=np.quantile(v, 0.1), p90=np.quantile(v, 0.9)))
            d = (vals['within'] - vals['start']) / ode_val
            rows.append(dict(output=k, dt=lab, timing='within - start (paired)', n_reps=len(d), ode=ode_val, ratio_mean=d.mean(), ratio_se=d.std(ddof=1) / np.sqrt(len(d)),
                             z=d.mean() / (d.std(ddof=1) / np.sqrt(len(d)))))
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(outdir, 'summary.csv'), index=False)

    summ = {f'TBExact {lab}, within': tbode.summarize(oo) for lab, oo in abms.items()}
    summ.update({f'TBExact {lab}, start (exp 05)': tbode.summarize(out05['abms'][lab]) for lab in dts})
    keys = ['INFECTION', 'CLEARED', 'NON_INFECTIOUS', 'ASYMPTOMATIC', 'prev_active', 'incidence_py', 'tb_deaths_py', 'frac_asymp']
    fig = tbode.plot_comparison({'equilibrium': ode_run}, summ, keys=keys, title=f'TBExact infection timing at equilibrium, ODE beta={beta:.3f} (500/100k), {n_seeds} reps per arm',
                                scale=scale)
    fig.savefig(os.path.join(outdir, 'infection_timing.png'), dpi=110)

    with pd.option_context('display.width', 250, 'display.float_format', '{:.4g}'.format):
        print(df[['output', 'dt', 'timing', 'n_reps', 'ode', 'mean', 'se', 'ratio_mean', 'ratio_se', 'z', 'median', 'p10', 'p90']].to_string(index=False))
    T.toc('Total')
