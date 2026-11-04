"""
Experiment 12, v0.10 side: the TBODE trajectories and equilibrium that the v4 ODE is checked against. Run in the v0.10 environment. See README.md.
"""
import os
import json
import numpy as np
import sciris as sc
import tbode

outdir = os.path.dirname(os.path.abspath(__file__))
beta_A = 9.917337993130424 # Exp 03's calibration (tbsim defaults, 500/100k)
setB = json.load(open(os.path.join(outdir, '..', '08_progression_tsi', 'calibration.json')))['B']
setB_pars = {k: v for k, v in setB['pars'].items() if k not in ('k_asy', 'k_non', 'replace_tb_deaths', 'init_sym')} # Set B's natural history, at k = 0
ode_kw = dict(replace_tb_deaths=False)


def trajectory(beta, **pars):
    """ State shares from a 1% symptomatic seed, 1900-2100, monthly """
    ode = tbode.TBODE(beta=beta, init_sym=0.01, **ode_kw, **pars)
    res = ode.run(start=1900, stop=2100, dt=1/12)
    return dict(year=res.time.tolist(), **{s: (res[s]/res.N).tolist() for s in tbode.states})


if __name__ == '__main__':
    beta = tbode.calibrate_beta(500e-5, **ode_kw)
    ode = tbode.TBODE(beta=beta, **ode_kw)
    x = ode.equilibrium()
    eq = {s: float(x[ode.idx(s)]) for s in tbode.states}
    out = dict(
        beta_A = beta,
        beta_A_exp03 = beta_A,
        equilibrium_A = eq,
        prev_active_A = float(ode.outputs(x).prev_active[0]),
        pars_B = setB_pars,
        beta_B = setB['beta'],
        trajectory_A = trajectory(beta),
        trajectory_B = trajectory(setB['beta'], **setB_pars),
    )
    sc.savejson(os.path.join(outdir, 'reference_v010.json'), out)
    tbode.save_versions(outdir, 'v0.10')
    os.replace(os.path.join(outdir, 'versions.json'), os.path.join(outdir, 'versions_v010.json'))
    print(f'beta = {beta:.6f} (exp 03: {beta_A:.6f}); prevalence {out["prev_active_A"]:.6f}')
