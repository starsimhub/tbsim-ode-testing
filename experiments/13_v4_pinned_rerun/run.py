"""
Experiment 13: Cliff's v4 redo (v4/common.py, calibrate.py, run_*.py) rerun at the pinned versions, with 40 reps per equilibrium arm and a daily
step for A0. Run in the v4 environment, from the repo root. See README.md.

Cliff's scripts are used unchanged, with their results redirected to this folder by TB_V4_FOLDER (inherited by worker processes). The
equilibrium runs use run_equilibrium.run_one with this experiment's design; the one-step, growth and timing scripts run as they are.
Rerunning resumes: equilibrium arms already saved with enough reps are skipped.
"""
import os
import sys
import subprocess
import numpy as np
import sciris as sc

outdir = sc.thispath(__file__)
v4dir = outdir.parent.parent / 'v4'
os.environ['TB_V4_FOLDER'] = str(outdir) # Before importing Cliff's code, so that it and its workers write to this folder
(outdir / 'results').mkdir(exist_ok=True)
sys.path.insert(0, str(v4dir))
import common as c
import run_equilibrium as req
import tbode4

reps = 40
more = 200 # For the arms whose first 40 reps were about 2 SE from the ODE: seeds 40-199 added (the first 40 seeds are Cliff's, so those reps repeat his runs)
design = [ # (arm, dt, network, reps): every arm of Cliff's run_equilibrium.py, plus A0 at a daily step
    ('A0', 'day', None, reps), ('A0', 'week', None, reps), ('A0', 'month', None, more), ('A0', 'month', 'net', reps), ('A0', 'month', 'net_matched', reps),
    ('B0', 'week', None, more), ('B0', 'month', None, reps),
    ('Ak', 'month', None, reps),
    ('Bk', 'week', None, more), ('Bk', 'month', None, reps),
    ('Bremote', 'week', None, reps), ('Bremote', 'month', None, reps),
]


def script(name, *args):
    """ Run one of Cliff's scripts as it is, in this environment, with its results going to this folder """
    T = sc.timer()
    subprocess.run([sys.executable, str(v4dir / name), *args], cwd=v4dir, env=os.environ, check=True)
    print(f'{name}: {T.toc(output=True):.0f} s', flush=True)


def run_equilibrium():
    """ As run_equilibrium.py's main block, with this experiment's design """
    fn = outdir / 'results' / 'equilibrium.obj'
    data = sc.load(fn) if fn.exists() else sc.objdict(ode=sc.objdict(), abm=sc.objdict(), versions=c.versions())
    for arm in sorted({d[0] for d in design}): # The ODE, at a monthly step, as the reference
        if arm not in data.ode:
            data.ode[arm] = req.run_one(arm, 'month', None, 0, solver='ode')
    for arm, dt, network, n in design:
        key = f'{arm}_{dt}' + (f'_{network}' if network else '')
        done = data.abm.get(key, [])
        if len(done) >= n:
            continue
        T = sc.timer()
        new = sc.parallelize(req.run_one, iterkwargs=[dict(arm=arm, dt=dt, network=network, seed=s) for s in range(len(done), n)], ncpus=c.ncpus) # Seeds after those already run
        data.abm[key] = list(done) + list(new)
        sc.save(fn, data)
        prev = np.array([o.prevalence[-50:].mean() for o in data.abm[key]])/data.ode[arm].prevalence[-50:].mean()
        print(f'{key:20s} {n} reps ({n - len(done)} new) in {T.toc(output=True):.0f} s; prevalence 1950–2000, agents ÷ ODE: {prev.mean():.3f} ± {prev.std(ddof=1)/np.sqrt(n):.3f}', flush=True)
    return data


def summarize(data):
    """ Agents ÷ ODE for every arm and output (Cliff's analyze.ratios: mean, SE, 10th and 90th percentiles over reps) """
    import analyze
    df = analyze.ratios(data)
    df.to_csv(outdir / 'summary.csv', index=False)
    return df


if __name__ == '__main__':
    T = sc.timer()
    tbode4.save_versions(outdir, 'v4')
    os.replace(outdir / 'versions.json', outdir / 'versions_v4.json')
    script('calibrate.py')
    data = run_equilibrium()
    df = summarize(data)
    print(df[df.output == 'prevalence'].to_string(index=False), flush=True)
    script('run_onestep.py')
    script('run_growth.py')
    script('run_timing.py') # Last, and alone, so that other runs don't slow it
    T.toc('Experiment 13')
