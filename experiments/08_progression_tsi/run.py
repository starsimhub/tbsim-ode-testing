"""
Experiment 08: progression by time since infection (k_asy = k_non = 5.846), ODE with latent classes vs TBExact and tbsim.TB,
for parameter sets A (tbsim defaults) and B (tb_LAI_TPT natural history). See README.md.
"""

import os
import json
import numpy as np
import pandas as pd
import sciris as sc
import matplotlib
matplotlib.use('Agg')
import tbsim
import tbode

outdir = os.path.dirname(os.path.abspath(__file__))
lai_csv = os.path.join(outdir, '..', '..', '..', 'tb_LAI_TPT', 'results', 'hm_real_20260728', 'trajectory_selection', 'collapsed_samples_eps1.2.csv')
start, stop = 1900, 2000
n_agents, n_seeds = 50_000, 40
target = 500e-5
dts = {'1 week': 7/365, '1 month': 1/12}  # Weekly first, so the longest runs start first
models = {'TBExact': (tbode.TBExact, dict(infection_timing='within')), 'tbsim.TB': (tbsim.TB, {})}
outputs = ['INFECTION', 'CLEARED', 'NON_INFECTIOUS', 'ASYMPTOMATIC', 'SYMPTOMATIC', 'prev_active', 'incidence_py', 'tb_deaths_py', 'frac_asymp']
scale = dict(prev_active=1e5, incidence_py=1e5, tb_deaths_py=1e5)


def par_sets():
    """ODE parameters (without beta) for sets A and B."""
    row = pd.read_csv(lai_csv).iloc[0]
    k = float(row.k_asy)
    assert np.isclose(k, row.k_non)
    base = dict(k_asy=k, k_non=k, replace_tb_deaths=False, init_sym=0.0)
    pi = float(row.rr_reinfection_rec)
    rates = ['inf_cle', 'inf_non', 'inf_asy', 'non_rec', 'non_asy', 'asy_non', 'asy_sym', 'sym_asy', 'trans_asymp']
    B = dict(base, **{r: float(row[r]) for r in rates}, sym_dead=float(row.sym_dead_1995), rr_reinfection_cleared=pi, rr_reinfection_rec=pi, rr_reinfection_treat=pi)
    return dict(A=base, B=B)


def calibrate(pars):
    """beta for 500/100k, the epidemic threshold, the prevalence elasticity and the ODE equilibrium outputs."""
    beta = tbode.calibrate_beta(target, bounds=(5, 100), **pars)
    ode = tbode.TBODE(beta=beta, **pars)
    x = ode.equilibrium()
    up = tbode.TBODE(beta=beta * 1.001, **pars)
    elasticity = (up.outputs(up.equilibrium()).prev_active[0] / target - 1) / 0.001
    threshold = tbode.calibrate_beta(1e-7, bounds=(1, beta), **pars)
    cohort = tbode.ode_cohort_progression(ode.pars, k_asy=pars['k_asy'], k_non=pars['k_non'])
    eq = tbode.ode_outputs(ode.outputs(x, time=np.zeros(1)))
    return dict(beta=beta, threshold=threshold, elasticity=elasticity, x=x, n_latent=ode.n_latent, latent_weights=ode.latent_weights.tolist(),
                equilibrium={k: float(eq[k][0]) for k in outputs}, ProgressionY1=cohort.ProgressionY1, ProgressionRatio=cohort.ProgressionRatio)


def run_abm(pars, beta, x, model, dt, seed):
    ode = tbode.TBODE(beta=beta, **pars)
    tb_class, kw = models[model]
    sim = tbode.make_abm(ode.pars, seed=seed, n_agents=n_agents, start=start, stop=stop, dt=dt, tb_class=tb_class, **kw)
    sim.init()
    tbode.init_from_ode(sim, x, seed=seed, ode=ode)
    sim.run()
    return tbode.annualize(tbode.abm_outputs(sim))


def last50(oo, k):
    """Each rep's mean over the last 50 years, in display units."""
    return np.array([np.nanmean(o[k][o.time >= stop - 50]) for o in oo]) * scale.get(k, 1)


if __name__ == '__main__':
    T = sc.timer()
    sets = par_sets()
    calib = {s: calibrate(p) for s, p in sets.items()}
    for s, c in calib.items():
        print(f"Set {s}: beta={c['beta']:.3f}, threshold={c['threshold']:.3f}, elasticity={c['elasticity']:.2f}, {c['n_latent']} latent classes, "
              f"ProgressionY1={c['ProgressionY1']:.4f}, ProgressionRatio={c['ProgressionRatio']:.3f}")
    json.dump({s: dict(pars=sets[s], **{k: v for k, v in c.items() if k != 'x'}) for s, c in calib.items()}, open(os.path.join(outdir, 'calibration.json'), 'w'), indent=2)
    T.toc('Calibration')

    arms = [(s, m, d) for s in sets for m in models for d in dts]
    iterkwargs = [dict(pars=sets[s], beta=calib[s]['beta'], x=calib[s]['x'], model=m, dt=dts[d], seed=seed) for (s, m, d) in arms for seed in range(n_seeds)]
    outs = sc.parallelize(run_abm, iterkwargs=iterkwargs, ncpus=tbode.ncpus)
    abms = {arm: outs[i * n_seeds:(i + 1) * n_seeds] for i, arm in enumerate(arms)}
    ode_runs = {s: tbode.ode_outputs(tbode.TBODE(beta=calib[s]['beta'], **sets[s]).run(start=start, stop=stop, dt=1/52, y0=calib[s]['x'] * 1e5)) for s in sets}
    sc.save(os.path.join(outdir, 'outputs.obj'), dict(abms=abms, ode_runs=ode_runs, calib=calib, sets=sets))
    T.toc('ABM runs')

    # Mean +/- SE across reps of each rep's last-50-year mean: each model / ODE, and tbsim.TB / TBExact (ratio of means, SE by the delta method)
    rows = []
    for k in outputs:
        for s in sets:
            ode_val = calib[s]['equilibrium'][k] * scale.get(k, 1)
            for d in dts:
                vals = {m: last50(abms[(s, m, d)], k) for m in models}
                for m, v in vals.items():
                    se = v.std(ddof=1) / np.sqrt(len(v))
                    rows.append(dict(output=k, set=s, dt=d, model=m, n_reps=len(v), ode=ode_val, mean=v.mean(), se=se, ratio_mean=v.mean() / ode_val,
                                     ratio_se=se / ode_val, z=(v.mean() - ode_val) / se, median=np.median(v), p10=np.quantile(v, 0.1), p90=np.quantile(v, 0.9)))
                a, b = vals['tbsim.TB'], vals['TBExact']
                r = a.mean() / b.mean()
                r_se = r * np.sqrt((a.std(ddof=1) / a.mean())**2 / len(a) + (b.std(ddof=1) / b.mean())**2 / len(b))
                rows.append(dict(output=k, set=s, dt=d, model='tbsim.TB / TBExact', n_reps=len(a), ratio_mean=r, ratio_se=r_se))
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(outdir, 'summary.csv'), index=False)

    keys = ['INFECTION', 'CLEARED', 'NON_INFECTIOUS', 'ASYMPTOMATIC', 'prev_active', 'incidence_py', 'tb_deaths_py', 'frac_asymp']
    for s in sets:
        summ = {f'{m} {d}': tbode.summarize(abms[(s, m, d)]) for m in models for d in dts}
        fig = tbode.plot_comparison({'equilibrium': ode_runs[s]}, summ, keys=keys, scale=scale,
                                    title=f"Set {s}, k = {sets[s]['k_asy']:.3f}: ODE beta = {calib[s]['beta']:.3f} (500/100k), {n_seeds} reps per arm")
        fig.savefig(os.path.join(outdir, f'progression_set{s}.png'), dpi=110)

    with pd.option_context('display.width', 250, 'display.float_format', '{:.4g}'.format):
        print(df[['output', 'set', 'dt', 'model', 'n_reps', 'ode', 'mean', 'se', 'ratio_mean', 'ratio_se', 'z', 'median', 'p10', 'p90']].to_string(index=False))
    T.toc('Total')
