"""
Experiment 03: baseline with beta calibrated to 500/100k prevalence; ABM started at the ODE equilibrium. See README.md.
"""

import os
import numpy as np
import pandas as pd
import sciris as sc
import matplotlib
matplotlib.use('Agg')
import tbode

outdir = os.path.dirname(os.path.abspath(__file__))
target, target_lo, target_hi = 500e-5, 200e-5, 800e-5
start, stop = 1900, 2000
n_agents, n_seeds = 50_000, 10
dts = {'1 day': 1/365, '1 week': 7/365, '1 month': 1/12}
ode_kw = dict(replace_tb_deaths=False, init_sym=0.0)
prog_targets = dict(ProgressionY1=0.0523, ProgressionRatio=0.663)
outputs = ['prev_active', 'incidence_py', 'tb_deaths_py', 'frac_asymp']


def run_abm(beta, x, dt, seed):
    ode = tbode.TBODE(beta=beta, **ode_kw)
    sim = tbode.make_abm(ode.pars, seed=seed, n_agents=n_agents, start=start, stop=stop, dt=dt)
    sim.init()
    tbode.init_from_ode(sim, x, seed=seed)
    sim.run()
    return tbode.annualize(tbode.abm_outputs(sim))


def run_cohort(pars, dt, seed):
    return tbode.abm_cohort_progression(pars, dt=dt, n_agents=50_000, seed=seed)


if __name__ == '__main__':
    T = sc.timer()

    # 1. Calibrate beta on the ODE equilibrium
    beta = tbode.calibrate_beta(target, replace_tb_deaths=False)
    beta_lo = tbode.calibrate_beta(target_lo, replace_tb_deaths=False)
    beta_hi = tbode.calibrate_beta(target_hi, replace_tb_deaths=False)
    ode = tbode.TBODE(beta=beta, **ode_kw)
    x = ode.equilibrium()
    eq = ode.outputs(x)
    print(f'ODE beta: {beta:.4f} (200/100k: {beta_lo:.4f}, 800/100k: {beta_hi:.4f})')

    # 2. ABM runs from the ODE equilibrium
    iterkwargs = [dict(beta=beta, x=x, dt=dt, seed=s) for dt in dts.values() for s in range(n_seeds)]
    outs = sc.parallelize(run_abm, iterkwargs=iterkwargs, ncpus=tbode.ncpus)
    abms = {lab: outs[i * n_seeds:(i + 1) * n_seeds] for i, lab in enumerate(dts)}
    T.toc('ABM runs')

    # 3. Cohort progression: ODE exact, ABM at each dt (2 seeds of 50k each)
    prog_ode = tbode.ode_cohort_progression(ode.pars)
    cohort_kw = [dict(pars=ode.pars, dt=dt, seed=s) for dt in dts.values() for s in range(2)]
    cohorts = sc.parallelize(run_cohort, iterkwargs=cohort_kw, ncpus=tbode.ncpus)
    prog_abm = {}
    for i, lab in enumerate(dts):
        B = np.mean([c.B for c in cohorts[i * 2:(i + 1) * 2]], axis=0)
        prog_abm[lab] = dict(ProgressionY1=B[0], ProgressionRatio=B[0] / B[9], B10=B[9])

    # 4. Steady-state table: ODE equilibrium vs ABM mean over the last 50 years (median and 10-90% across seeds)
    ode_run = tbode.ode_outputs(tbode.TBODE(beta=beta, **ode_kw).run(start=start, stop=stop, dt=1/52, y0=x * 1e5))
    rows = []
    scale = dict(prev_active=1e5, incidence_py=1e5, tb_deaths_py=1e5, frac_asymp=1)
    for k in outputs:
        row = dict(output=k, ode_equilibrium=eq[k][0] * scale[k])
        for lab, oo in abms.items():
            vals = np.array([np.nanmean(o[k][o.time >= stop - 50]) for o in oo]) * scale[k]
            first = np.array([np.nanmean(o[k][o.time < start + 10]) for o in oo]) * scale[k]
            row[f'abm {lab}: first 10y'] = np.median(first)
            row[f'abm {lab}: last 50y median'] = np.median(vals)
            row[f'abm {lab}: last 50y p10'] = np.quantile(vals, 0.1)
            row[f'abm {lab}: last 50y p90'] = np.quantile(vals, 0.9)
            row[f'abm {lab}: ratio to ODE'] = np.median(vals) / (eq[k][0] * scale[k])
        rows.append(row)
    for k, tgt in prog_targets.items():
        row = dict(output=k, ode_equilibrium=prog_ode[k], target=tgt)
        for lab in dts:
            row[f'abm {lab}: last 50y median'] = prog_abm[lab][k]
            row[f'abm {lab}: ratio to ODE'] = prog_abm[lab][k] / prog_ode[k]
        rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(outdir, 'summary.csv'), index=False)

    # 5. ODE beta that reproduces each ABM's long-run prevalence ("effective beta")
    eff = {}
    for lab, oo in abms.items():
        prev = np.median([np.nanmean(o.prev_active[o.time >= stop - 50]) for o in oo])
        eff[lab] = (prev, tbode.calibrate_beta(prev, replace_tb_deaths=False))
    sc.savejson(os.path.join(outdir, 'calibration.json'), dict(
        beta=beta, beta_200=beta_lo, beta_800=beta_hi,
        equilibrium={k: float(eq[k][0]) for k in outputs + ['INFECTION', 'NON_INFECTIOUS', 'ASYMPTOMATIC', 'SYMPTOMATIC']},
        prevalence_to_incidence=float(eq.prev_active[0] / eq.incidence_py[0]),
        mortality_to_incidence=float(eq.tb_deaths_py[0] / eq.incidence_py[0]),
        progression_ode=dict(ProgressionY1=prog_ode.ProgressionY1, ProgressionRatio=prog_ode.ProgressionRatio, B10=prog_ode.B[1]),
        progression_abm=prog_abm,
        effective_beta={lab: dict(prev=v[0], beta=v[1]) for lab, v in eff.items()},
    ))

    # 6. Figure, with its underlying data saved so it can be redrawn without rerunning
    sc.save(os.path.join(outdir, 'outputs.obj'), dict(abms=abms, ode_run=ode_run, cohorts=prog_abm))
    summ ={f'dt={lab}': tbode.summarize(oo) for lab, oo in abms.items()}
    keys = ['INFECTION', 'CLEARED', 'NON_INFECTIOUS', 'ASYMPTOMATIC', 'prev_active', 'incidence_py', 'tb_deaths_py', 'frac_asymp']
    fig = tbode.plot_comparison({'equilibrium': ode_run}, summ, keys=keys, title=f'Baseline at equilibrium, ODE beta={beta:.3f} (500/100k)',
                                scale=dict(prev_active=1e5, incidence_py=1e5, tb_deaths_py=1e5))
    fig.savefig(os.path.join(outdir, 'baseline_equilibrium.png'), dpi=110)

    with pd.option_context('display.width', 250, 'display.max_columns', 30, 'display.float_format', '{:.4g}'.format):
        print(df.T.to_string())
    print('Effective ODE beta of each ABM:', {k: (round(v[0] * 1e5, 1), round(v[1], 4)) for k, v in eff.items()})
    T.toc('Total')
