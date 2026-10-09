"""
Experiment 06: run time of tbode.TBExact vs tbsim.TB at daily, weekly and monthly dt. See README.md.

Runs serially in one process (no parallelism), so arms don't compete for cores.
"""

import os
import json
import time
import numpy as np
import pandas as pd
import sciris as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import tbsim
import tbode

outdir = os.path.dirname(os.path.abspath(__file__))
beta = json.load(open(os.path.join(outdir, '..', '03_baseline_calibrated', 'calibration.json')))['beta']
start, stop = 1900, 1910
n_agents, n_reps = 50_000, 3
dts = {'1 day': 1/365, '1 week': 7/365, '1 month': 1/12}
models = {'tbsim.TB': tbsim.TB, 'TBExact': tbode.TBExact}
ode_kw = dict(replace_tb_deaths=False, init_sym=0.0)


def time_run(model, dt, seed, x):
    """Build, initialise and run one sim; return wall times and the per-loop-entry CPU breakdown."""
    ode = tbode.TBODE(beta=beta, **ode_kw)
    t0 = time.perf_counter()
    sim = tbode.make_abm(ode.pars, seed=seed, n_agents=n_agents, start=start, stop=stop, dt=dt, tb_class=models[model])
    sim.init()
    tbode.init_from_ode(sim, x, seed=seed)
    t1 = time.perf_counter()
    sim.run(profile=True)
    t2 = time.perf_counter()
    sim.loop.to_df()
    n_steps = len(sim.t.tvec)
    row = dict(model=model, dt=[k for k, v in dts.items() if v == dt][0], seed=seed, n_steps=n_steps,
               init_s=t1 - t0, run_s=t2 - t1, ms_per_step=(t2 - t1) / n_steps * 1e3)
    cpu = sim.loop.cpu_df.reset_index()[['label', 'module', 'calls', 'cpu_time', 'percent']]
    for k in ['model', 'dt', 'seed']:
        cpu.insert(0, k, row[k])
    return row, cpu


if __name__ == '__main__':
    T = sc.timer()
    x = tbode.TBODE(beta=beta, **ode_kw).equilibrium()
    time_run('tbsim.TB', 1/12, 999, x)  # Warm-up (imports, caches); discarded

    rows, cpus = [], []
    for seed in range(n_reps):          # Interleave arms within each rep, so drift in machine load affects all arms
        for lab, dt in dts.items():
            for model in models:
                row, cpu = time_run(model, dt, seed, x)
                rows.append(row)
                cpus.append(cpu)
                print(f'{model:9s} {lab:8s} seed {seed}: {row["run_s"]:7.1f} s run, {row["ms_per_step"]:6.2f} ms/step, {row["init_s"]:5.1f} s init')
    df = pd.DataFrame(rows)
    cpu = pd.concat(cpus, ignore_index=True)
    df.to_csv(os.path.join(outdir, 'timings.csv'), index=False)
    cpu.to_csv(os.path.join(outdir, 'cpu_breakdown.csv'), index=False)

    # Summary: median ms/step, implied time per 100-year run, and TBExact / tbsim.TB ratio
    summ = df.groupby(['dt', 'model'], sort=False).agg(ms_per_step=('ms_per_step', 'median'), ms_min=('ms_per_step', 'min'), ms_max=('ms_per_step', 'max'),
                                                         init_s=('init_s', 'median')).reset_index()
    summ['s_per_100y'] = summ.ms_per_step / 1e3 * summ.dt.map({k: 100 / v for k, v in dts.items()})
    ref = summ[summ.model == 'tbsim.TB'].set_index('dt').ms_per_step
    summ['ratio_to_tbsim_TB'] = summ.ms_per_step / summ.dt.map(ref)
    summ.to_csv(os.path.join(outdir, 'summary.csv'), index=False)

    # Per-module share of run time (median over reps): top entries by arm
    top = cpu.groupby(['dt', 'model', 'label'], sort=False).cpu_time.median().reset_index()
    top['ms_per_step'] = top.cpu_time / top.dt.map(df.groupby('dt').n_steps.first()) * 1e3
    top = top.sort_values(['dt', 'model', 'ms_per_step'], ascending=[True, True, False]).groupby(['dt', 'model']).head(6)

    # Figure: ms per step by arm (each rep shown), from timings.csv
    fig, ax = plt.subplots(figsize=(7, 4))
    for i, model in enumerate(models):
        d = df[df.model == model]
        xpos = np.array([list(dts).index(v) for v in d.dt]) + (i - 0.5) * 0.3
        ax.scatter(xpos, d.ms_per_step, color=f'C{i}', label=model, zorder=3)
        med = summ[summ.model == model]
        ax.bar(np.arange(len(dts)) + (i - 0.5) * 0.3, med.ms_per_step, width=0.28, color=f'C{i}', alpha=0.35)
    ax.set_xticks(range(len(dts)), list(dts))
    ax.set_ylabel('Run time per step (ms)')
    ax.set_title(f'Run time per step, {n_agents:,} agents (bars: median of {n_reps} reps)')
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, 'timings.png'), dpi=110)

    with pd.option_context('display.width', 250, 'display.float_format', '{:.3g}'.format):
        print(summ.to_string(index=False))
        print(top.to_string(index=False))
    T.toc('Total')
