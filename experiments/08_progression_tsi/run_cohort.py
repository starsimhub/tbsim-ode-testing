"""
Experiment 08, cohort progression: ProgressionY1 / ProgressionRatio from ABM cohorts (TBExact, tbsim.TB) vs the exact ODE cohort. See README.md.
"""

import os
import numpy as np
import pandas as pd
import sciris as sc
import tbsim
import tbode
from run import par_sets, outdir

n_agents, n_seeds = 100_000, 4
dts = {'1 day': 1/365, '1 week': 7/365, '1 month': 1/12}
arms = [('TBExact', '1 week'), ('TBExact', '1 month'), ('tbsim.TB', '1 day'), ('tbsim.TB', '1 week'), ('tbsim.TB', '1 month')]
classes = {'TBExact': tbode.TBExact, 'tbsim.TB': tbsim.TB}


def run_cohort(pars, model, dt, seed):
    ode = tbode.TBODE(**pars)
    return tbode.abm_cohort_progression(ode.pars, dt=dt, n_agents=n_agents, seed=seed, tb_class=classes[model])


if __name__ == '__main__':
    T = sc.timer()
    sets = par_sets()
    iterkwargs = [dict(pars=sets[s], model=m, dt=dts[d], seed=seed) for s in sets for (m, d) in arms for seed in range(n_seeds)]
    outs = sc.parallelize(run_cohort, iterkwargs=iterkwargs, ncpus=tbode.ncpus)
    keys = [(s, m, d) for s in sets for (m, d) in arms]
    cohorts = {key: outs[i * n_seeds:(i + 1) * n_seeds] for i, key in enumerate(keys)}
    exact = {s: tbode.ode_cohort_progression(tbode.TBODE(**p).pars, k_asy=p['k_asy'], k_non=p['k_non']) for s, p in sets.items()}
    sc.save(os.path.join(outdir, 'cohort.obj'), dict(cohorts=cohorts, exact=exact))

    # Pool the seeds (equal sizes): B(1), B(10) and their ratio, each / exact, with binomial SEs (delta method for the ratio)
    rows = []
    for (s, m, d), cc in cohorts.items():
        n = n_agents * len(cc)
        B1, B10 = np.mean([c.B[0] for c in cc]), np.mean([c.B[9] for c in cc])
        se1, se10 = np.sqrt(B1 * (1 - B1) / n), np.sqrt(B10 * (1 - B10) / n)
        ex = exact[s]
        ratio = B1 / B10
        cov = B1 * (1 - B10) / n  # B1 events are a subset of B10's
        se_ratio = ratio * np.sqrt(max(se1**2 / B1**2 + se10**2 / B10**2 - 2 * cov / (B1 * B10), 0))
        rows.append(dict(set=s, model=m, dt=d, n_agents=n, ProgressionY1=B1, Y1_vs_exact=B1 / ex.ProgressionY1, Y1_vs_exact_se=se1 / ex.ProgressionY1,
                         B10=B10, B10_vs_exact=B10 / ex.B[1], B10_vs_exact_se=se10 / ex.B[1],
                         ProgressionRatio=ratio, ratio_vs_exact=ratio / ex.ProgressionRatio, ratio_vs_exact_se=se_ratio / ex.ProgressionRatio))
    for s, ex in exact.items():
        rows.append(dict(set=s, model='exact (ODE)', dt='', ProgressionY1=ex.ProgressionY1, B10=ex.B[1], ProgressionRatio=ex.ProgressionRatio))
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(outdir, 'cohort.csv'), index=False)
    with pd.option_context('display.width', 250, 'display.float_format', '{:.4g}'.format):
        print(df.to_string(index=False))
    T.toc('Total')
