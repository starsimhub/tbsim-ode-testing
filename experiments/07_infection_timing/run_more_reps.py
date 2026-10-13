"""
Experiment 07, more reps: weekly dt, for both 'within' and 'start' timing, for a range of seeds. See README.md.

Usage: python run_more_reps.py FIRST_SEED STOP_SEED (e.g. 40 80). Saves outputs_more_reps_{first}-{last}.obj, then pools
every saved batch with seeds 0-39 (exp 07 'within', exp 05 'start') into summary_more_reps.csv.

The first 40 weekly 'within' reps came out ~5% low in every output (prevalence 0.95 +/- 0.02) while 'start' was 1.00;
the extra reps (seeds 40-79, then 80-159) test whether that is rep-to-rep drift. Before running, checks that the current TBExact (which has since
gained latent classes for k > 0) reproduces exp 07's 'within' and exp 05's 'start' runs exactly at k = 0.
"""

import os
import sys
import glob
import numpy as np
import pandas as pd
import sciris as sc
from run import run_abm, last50, beta, ode_kw, outdir, dir05, dts, n_seeds, outputs, scale
import tbode

seeds = range(int(sys.argv[1]), int(sys.argv[2])) if __name__ == '__main__' else None
dt_lab = '1 week'


if __name__ == '__main__':
    T = sc.timer()
    ode = tbode.TBODE(beta=beta, **ode_kw)
    x = ode.equilibrium()
    eq_states = tbode.ode_outputs(ode.outputs(x, time=np.zeros(1)))
    out07 = sc.load(os.path.join(outdir, 'outputs.obj'))
    out05 = sc.load(os.path.join(dir05, 'outputs.obj'))

    checks = sc.parallelize(run_abm, iterkwargs=[dict(x=x, dt=dts[dt_lab], seed=0, infection_timing='within'),
                                                 dict(x=x, dt=dts['1 month'], seed=0, infection_timing='start')], ncpus=2)
    for chk, ref, lab in [(checks[0], out07['abms'][dt_lab][0], "exp 07 'within' weekly"), (checks[1], out05['abms']['1 month'][0], "exp 05 'start' monthly")]:
        for k in outputs:
            assert np.allclose(chk[k], ref[k], equal_nan=True), f'Current TBExact does not reproduce {lab} for {k}'
    print("Check passed: current TBExact reproduces exp 07 ('within', weekly) and exp 05 ('start', monthly), seed 0")

    iterkwargs = [dict(x=x, dt=dts[dt_lab], seed=s, infection_timing=timing) for timing in ['within', 'start'] for s in seeds]
    outs = sc.parallelize(run_abm, iterkwargs=iterkwargs, ncpus=tbode.ncpus)
    more = {'within': outs[:len(seeds)], 'start': outs[len(seeds):]}
    sc.save(os.path.join(outdir, f'outputs_more_reps_{seeds[0]}-{seeds[-1]}.obj'), dict(seeds=list(seeds), dt=dt_lab, abms=more))
    T.toc('ABM runs')

    # All weekly reps per timing: seeds 0-39 (exp 07 'within', exp 05 'start') plus every batch saved here, in seed order
    batches = sorted([sc.load(f) for f in glob.glob(os.path.join(outdir, 'outputs_more_reps_*.obj'))], key=lambda b: b['seeds'][0])
    allreps = {'within': list(out07['abms'][dt_lab]), 'start': list(out05['abms'][dt_lab])}
    batch_labels = [('seeds 0-39', 0, n_seeds)]
    for b in batches:
        assert b['seeds'][0] == n_seeds + sum(len(bb['seeds']) for bb in batches if bb['seeds'][0] < b['seeds'][0]), 'seed batches must be contiguous'
        batch_labels.append((f"seeds {b['seeds'][0]}-{b['seeds'][-1]}", len(allreps['within']), len(allreps['within']) + len(b['seeds'])))
        for timing in allreps:
            allreps[timing] += b['abms'][timing]
    rows = []
    for k in outputs:
        ode_val = eq_states[k][0] * scale.get(k, 1)
        vals = {timing: last50(oo, k) for timing, oo in allreps.items()}
        for timing, v in vals.items():
            se = v.std(ddof=1) / np.sqrt(len(v))
            rows.append(dict(output=k, dt=dt_lab, timing=timing, n_reps=len(v), ode=ode_val, mean=v.mean(), se=se, ratio_mean=v.mean() / ode_val,
                             ratio_se=se / ode_val, z=(v.mean() - ode_val) / se, median=np.median(v), p10=np.quantile(v, 0.1), p90=np.quantile(v, 0.9)))
        d = (vals['within'] - vals['start']) / ode_val
        se_d = d.std(ddof=1) / np.sqrt(len(d))
        rows.append(dict(output=k, dt=dt_lab, timing='within - start (paired)', n_reps=len(d), ode=ode_val, ratio_mean=d.mean(), ratio_se=se_d, z=d.mean() / se_d))
        for half, i0, i1 in batch_labels:
            sl = slice(i0, i1)
            for timing, v in vals.items():
                rows.append(dict(output=k, dt=dt_lab, timing=f'{timing}, {half}', n_reps=len(v[sl]), ode=ode_val, ratio_mean=v[sl].mean() / ode_val,
                                 ratio_se=v[sl].std(ddof=1) / np.sqrt(len(v[sl])) / ode_val))
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(outdir, 'summary_more_reps.csv'), index=False)
    with pd.option_context('display.width', 250, 'display.float_format', '{:.4g}'.format):
        print(df[['output', 'dt', 'timing', 'n_reps', 'ode', 'mean', 'se', 'ratio_mean', 'ratio_se', 'z', 'median', 'p10', 'p90']].to_string(index=False))
    T.toc('Total')
