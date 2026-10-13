"""
Experiment 10: TBExact with network vs mass-action transmission ('exp' or 'linear' per-step infection probability). See README.md.
"""

import os
import numpy as np
import pandas as pd
import sciris as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import tbode
from step_map import settings, equilibrium

outdir = os.path.dirname(os.path.abspath(__file__))
dir05 = os.path.join(outdir, '..', '05_exact_transitions')
dir09 = os.path.join(outdir, '..', '09_population_size')
start, stop = 1900, 2000
dts = {'1 week': 7/365, '1 month': 1/12}

# (label, setting, dt, transmission, infection_prob, timing, n_agents, n_reps, reused): longest runs first
arms = [
    ('7: k=0 weekly within, mass action exp', 'k0', '1 week', 'mass_action', 'exp', 'within', 50_000, 80, None),
    ('6: k=0 monthly start, network', 'k0', '1 month', 'network', 'exp', 'start', 50_000, 160, 'exp05'),
    ('1: A, network', 'A', '1 month', 'network', 'exp', 'within', 25_000, 190, 'exp09'),
    ('2: A, mass action exp', 'A', '1 month', 'mass_action', 'exp', 'within', 25_000, 190, None),
    ('3: A, mass action linear', 'A', '1 month', 'mass_action', 'linear', 'within', 25_000, 190, None),
    ('4: B, network', 'B', '1 month', 'network', 'exp', 'within', 25_000, 190, None),
    ('5: B, mass action linear', 'B', '1 month', 'mass_action', 'linear', 'within', 25_000, 190, None),
]
outputs = ['prev_active', 'INFECTION', 'incidence_py', 'tb_deaths_py']
scale = dict(prev_active=1e5, incidence_py=1e5, tb_deaths_py=1e5)


def run_abm(setting, dt, transmission, infection_prob, timing, n_agents, seed):
    pars = settings()[setting]
    ode = tbode.TBODE(**pars)
    x = ode.equilibrium()
    if transmission == 'network':
        kw = dict(tb_class=tbode.TBExact)
    else:
        kw = dict(tb_class=tbode.TBExactMassAction, network=False, beta_ode=pars['beta'], infection_prob=infection_prob)
    sim = tbode.make_abm(ode.pars, seed=seed, n_agents=n_agents, start=start, stop=stop, dt=dts[dt], infection_timing=timing, **kw)
    sim.init()
    tbode.init_from_ode(sim, x, seed=seed, ode=ode)
    sim.run()
    return tbode.annualize(tbode.abm_outputs(sim))


def last50(oo, k):
    """Each rep's mean over the last 50 years, in display units."""
    return np.array([np.nanmean(o[k][o.time >= stop - 50]) for o in oo]) * scale.get(k, 1)


def reused(arm):
    """Saved runs reused for an arm (seeds 0..n-1), or []."""
    if arm[-1] == 'exp05':
        return list(sc.load(os.path.join(dir05, 'outputs.obj'))['abms']['1 month'])
    if arm[-1] == 'exp09':
        return list(sc.load(os.path.join(dir09, 'outputs.obj'))['abms'][25_000])
    return []


if __name__ == '__main__':
    T = sc.timer()
    keys = ['setting', 'dt', 'transmission', 'infection_prob', 'timing', 'n_agents']

    # Reused arms: check that the current code reproduces seed 0 exactly
    for arm in arms:
        old = reused(arm)
        if old:
            new = run_abm(**dict(zip(keys, arm[1:7])), seed=0)
            for k in outputs:
                assert np.allclose(new[k], old[0][k], equal_nan=True), f"Current code doesn't reproduce {arm[-1]} for {arm[0]}, {k}"
            print(f'Check passed: current code reproduces {arm[-1]} seed 0 ({arm[0]})')

    iterkwargs, slices = [], {}
    for arm in arms:
        n_old = len(reused(arm))
        i0 = len(iterkwargs)
        iterkwargs += [dict(**dict(zip(keys, arm[1:7])), seed=s) for s in range(n_old, arm[7])]
        slices[arm[0]] = (i0, len(iterkwargs))
    outs = sc.parallelize(run_abm, iterkwargs=iterkwargs, ncpus=tbode.ncpus)
    abms = {arm[0]: reused(arm) + outs[slices[arm[0]][0]:slices[arm[0]][1]] for arm in arms}
    sc.save(os.path.join(outdir, 'outputs.obj'), dict(abms=abms, arms=arms))
    T.toc('ABM runs')

    # Summary: ABM / ODE (mean +/- SE across reps) and the deterministic step-map prediction for the arm's transmission form
    sets = settings()
    rows = []
    for arm in arms:
        lab, setting, dt, transmission, prob, timing, n_agents, n_reps, _ = arm
        ode = tbode.TBODE(**sets[setting])
        eq = tbode.ode_outputs(ode.outputs(ode.equilibrium(), time=np.zeros(1)))
        pred = {form: equilibrium(ode, dts[dt], form, timing)[0] for form in ['linear', 'exp']}
        oo = abms[lab]
        faded = int(sum(np.any(o.prev_active[o.time >= start + 1] == 0) for o in oo))
        for k in outputs:
            v = last50(oo, k) / (eq[k][0] * scale.get(k, 1))
            se = v.std(ddof=1) / np.sqrt(len(v))
            row = dict(arm=lab, output=k, setting=setting, dt=dt, transmission=transmission, infection_prob=prob, timing=timing, n_agents=n_agents,
                       n_reps=len(v), n_faded=faded, ratio_mean=v.mean(), ratio_se=se, sd_ratio=v.std(ddof=1))
            if k == 'prev_active':
                row.update(map_linear=pred['linear'], map_exp=pred['exp'], z_vs_map=(v.mean() - pred[prob]) / se)
            rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(outdir, 'summary.csv'), index=False)

    d = df[df.output == 'prev_active'].reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(11, 5))
    y = np.arange(len(d))[::-1]
    ax.errorbar(d.ratio_mean, y, xerr=d.ratio_se, fmt='o', color='k', label='TBExact (mean ± SE)')
    ax.scatter(d.map_linear, y, marker='|', s=200, color='C0', label="Deterministic step, linear (x)")
    ax.scatter(d.map_exp, y, marker='|', s=200, color='C3', label="Deterministic step, 1 − exp(−x)")
    ax.axvline(1, color='gray', ls=':', label='ODE')
    ax.set_yticks(y)
    ax.set_yticklabels([f"{a} ({n} reps)" for a, n in zip(d.arm, d.n_reps)])
    ax.set_xlabel('Prevalence, ABM ÷ ODE (last 50 years)')
    ax.legend(fontsize=8, frameon=False, loc='lower left')
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, 'transmission_step.png'), dpi=110)

    with pd.option_context('display.width', 250, 'display.float_format', '{:.4g}'.format):
        print(d[['arm', 'n_agents', 'n_reps', 'n_faded', 'ratio_mean', 'ratio_se', 'sd_ratio', 'map_linear', 'map_exp', 'z_vs_map']].to_string(index=False))
        print(df[df.output != 'prev_active'][['arm', 'output', 'ratio_mean', 'ratio_se']].to_string(index=False))
    T.toc('Total')
