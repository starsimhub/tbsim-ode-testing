"""
Experiment 09: exp 08's set A, monthly TBExact, at N = 25k, 50k, 100k and 200k agents. See README.md.
"""

import os
import json
import numpy as np
import pandas as pd
import sciris as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import tbode

outdir = os.path.dirname(os.path.abspath(__file__))
dir08 = os.path.join(outdir, '..', '08_progression_tsi')
start, stop = 1900, 2000
dt = 1/12
det_step = 0.9952  # Deterministic TBExact step ÷ ODE, set A monthly (exp 08 step_map.log)
reps = {200_000: 24, 100_000: 48, 50_000: 94, 25_000: 190}  # Largest first, so the longest runs start first
reuse50k = 40  # Seeds 0-39 at 50k come from exp 08
outputs = ['prev_active', 'INFECTION', 'incidence_py', 'tb_deaths_py']
scale = dict(prev_active=1e5, incidence_py=1e5, tb_deaths_py=1e5)


def setup():
    """Set A parameters and calibration from exp 08, and the ODE equilibrium."""
    c = json.load(open(os.path.join(dir08, 'calibration.json')))['A']
    ode = tbode.TBODE(beta=c['beta'], **c['pars'])
    return c, ode, ode.equilibrium()


def run_abm(n_agents, seed):
    c, ode, x = setup()
    sim = tbode.make_abm(ode.pars, seed=seed, n_agents=n_agents, start=start, stop=stop, dt=dt, tb_class=tbode.TBExact, infection_timing='within')
    sim.init()
    tbode.init_from_ode(sim, x, seed=seed, ode=ode)
    sim.run()
    return tbode.annualize(tbode.abm_outputs(sim))


def last50(oo, k):
    """Each rep's mean over the last 50 years, in display units."""
    return np.array([np.nanmean(o[k][o.time >= stop - 50]) for o in oo]) * scale.get(k, 1)


if __name__ == '__main__':
    T = sc.timer()
    c, ode, x = setup()
    eq = tbode.ode_outputs(ode.outputs(x, time=np.zeros(1)))

    # Exp 08's 50k set A monthly TBExact runs are reused: check the current code reproduces seed 0 exactly
    out08 = sc.load(os.path.join(dir08, 'outputs.obj'))
    ref = out08['abms'][('A', 'TBExact', '1 month')]
    check = run_abm(50_000, 0)
    for k in outputs:
        assert np.allclose(check[k], ref[0][k], equal_nan=True), f"Current TBExact doesn't reproduce exp 08 for {k}"
    print('Check passed: current TBExact reproduces exp 08 (set A, monthly, seed 0)')

    iterkwargs = [dict(n_agents=n, seed=s) for n, r in reps.items() for s in range(reuse50k if n == 50_000 else 0, r)]
    outs = sc.parallelize(run_abm, iterkwargs=iterkwargs, ncpus=tbode.ncpus)
    abms = {}
    i = 0
    for n, r in reps.items():
        k = r - (reuse50k if n == 50_000 else 0)
        abms[n] = (list(ref) if n == 50_000 else []) + outs[i:i + k]
        i += k
    sc.save(os.path.join(outdir, 'outputs.obj'), dict(abms=abms, det_step=det_step, calib=c))
    T.toc('ABM runs')

    rows = []
    for k in outputs:
        ode_val = eq[k][0] * scale.get(k, 1)
        for n in sorted(reps):
            oo = abms[n]
            faded = np.array([np.any(o.prev_active[o.time >= start + 1] == 0) for o in oo])
            for subset, keep in [('all reps', np.ones(len(oo), bool)), ('excluding fade-outs', ~faded)]:
                v = last50([o for o, kp in zip(oo, keep) if kp], k)
                se = v.std(ddof=1) / np.sqrt(len(v))
                rows.append(dict(output=k, n_agents=n, subset=subset, n_reps=len(v), n_faded=int(faded.sum()), ode=ode_val, mean=v.mean(), se=se,
                                 ratio_mean=v.mean() / ode_val, ratio_se=se / ode_val, sd_ratio=v.std(ddof=1) / ode_val,
                                 shortfall_vs_det=det_step - v.mean() / ode_val if k == 'prev_active' else np.nan))
    df = pd.DataFrame(rows)

    # Fits of the prevalence shortfall (vs the deterministic step) to c/N and c/sqrt(N), weighted by 1/SE^2
    d = df[(df.output == 'prev_active') & (df.subset == 'all reps')]
    fits = {}
    for lab, f in [('c/N', lambda n: 1e5 / n), ('c/sqrt(N)', lambda n: np.sqrt(1e5 / n))]:
        xx, yy, w = f(d.n_agents.values), d.shortfall_vs_det.values, 1 / d.ratio_se.values**2
        coef = np.sum(w * xx * yy) / np.sum(w * xx**2)
        chi2 = np.sum(w * (yy - coef * xx)**2)
        fits[lab] = dict(coef_at_100k=coef, chi2=chi2, dof=len(xx) - 1)
    yy, w = d.shortfall_vs_det.values, 1 / d.ratio_se.values**2
    const = np.sum(w * yy) / np.sum(w)
    fits['constant'] = dict(coef_at_100k=const, chi2=np.sum(w * (yy - const)**2), dof=len(yy) - 1)
    df.to_csv(os.path.join(outdir, 'summary.csv'), index=False)
    json.dump(fits, open(os.path.join(outdir, 'fits.json'), 'w'), indent=2)

    fig, axs = plt.subplots(1, 2, figsize=(13, 5))
    ax = axs[0]
    ax.errorbar(1e5 / d.n_agents, d.ratio_mean, yerr=d.ratio_se, fmt='o', color='k', label='TBExact (mean ± SE)')
    ax.axhline(1, color='gray', ls=':', label='ODE')
    ax.axhline(det_step, color='C0', ls='--', label='Deterministic TBExact step')
    xs = np.linspace(0, 4.2, 50)
    for j, lab in enumerate(['c/N', 'c/sqrt(N)']):
        fx = xs if lab == 'c/N' else np.sqrt(xs)
        ax.plot(xs, det_step - fits[lab]['coef_at_100k'] * fx, color=f'C{j + 1}', label=f'Fit: shortfall ∝ {lab}')
    ax.set_xlabel('100,000 / N')
    ax.set_ylabel('Prevalence, ABM ÷ ODE (last 50 years)')
    ax.set_title('Set A, monthly dt: prevalence vs population size')
    ax.legend(fontsize=8, frameon=False)
    ax = axs[1]
    for j, n in enumerate(sorted(reps)):
        s = tbode.summarize(abms[n])
        ax.plot(s.time, s.prev_active[1] * 1e5, color=f'C{j}', label=f'N = {n // 1000}k (median)')
    ax.axhline(eq.prev_active[0] * 1e5, color='k', ls='--', label='ODE')
    ax.set_ylim(bottom=0)
    ax.set_title('Prevalence per 100k over time')
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, 'population_size.png'), dpi=110)

    with pd.option_context('display.width', 250, 'display.float_format', '{:.4g}'.format):
        print(df.to_string(index=False))
    print(json.dumps(fits, indent=2))
    T.toc('Total')
