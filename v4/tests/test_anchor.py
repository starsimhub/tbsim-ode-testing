"""
The anchor between the tracks: tbsim v4's ODE reproduces the reference ODE and our v0.10 ODE (TBODE). Rerun after every bump of the v4 pins.

The TBODE values come from experiments/12_v4_setup/reference_v010.json, made in the v0.10 environment (export_reference.py there).
"""
from pathlib import Path
import numpy as np
import sciris as sc
import starsim as ss
import tbsim.compartmental as tbc
import tbode4

ref = sc.loadjson(Path(__file__).resolve().parents[2] / 'experiments' / '12_v4_setup' / 'reference_v010.json')


def test_reference_ode():
    """ With births replacing every death, v4 TB as an ODE is tbsim's TB_ODE (520 years, every state within 1e-4) """
    tb_ode = tbc.TB_ODE().run(start_time=1500, end_time=2020)
    sim = tbode4.make_sim(beta=9.0, start=1900, stop=2420, dt=ss.years(1), init=dict(symptomatic=1e3/1e5), replace_tb_deaths=True).run()
    for s in tbode4.states:
        v4 = np.asarray(sim.results.tb[s].values, dtype=float)
        assert np.allclose(v4, tb_ode[s.upper()], rtol=1e-4, atol=1e-2), f'{s} differs from TB_ODE'


def test_tbode_trajectories():
    """ From a 1% symptomatic seed, the v4 ODE's state shares are TBODE's (baseline and set B at k = 0): within 1e-4 relative after the first year, 1e-5 absolute throughout """
    for label in ['A', 'B']:
        pars = tbode4.odepars_to_tb(ref['pars_B']) if label == 'B' else {}
        v4 = tbode4.shares(tbode4.make_sim(beta=ref[f'beta_{label}'], start=1900, stop=2100, init=dict(symptomatic=0.01), **pars).run())
        tr = ref[f'trajectory_{label}']
        later = v4.year >= 1901
        for s in tbode4.states:
            b = np.array(tr[s.upper()])
            assert np.allclose(v4[s], b, rtol=0, atol=1e-5), f'{label} {s}: absolute difference'
            assert np.allclose(v4[s][later], b[later], rtol=1e-4, atol=1e-8), f'{label} {s}: relative difference after 1901'


def test_tbode_equilibrium():
    """ Started at TBODE's baseline equilibrium (beta calibrated to 500/100k), the v4 ODE stays there """
    init = {s: max(ref['equilibrium_A'][s.upper()], 0.0) for s in tbode4.states}
    v4 = tbode4.shares(tbode4.make_sim(beta=ref['beta_A'], start=1900, stop=2000, init=init).run())
    assert np.allclose(v4.prev_active, 500e-5, rtol=1e-8)


def test_beta_edge_round_trip():
    """ The matched per-contact beta gives the ODE's per-step infections on a network of mean degree K """
    beta_ode, K, dt = 12, 10, 1/12
    b = tbode4.beta_edge_from_ode(beta_ode, K, dt)
    assert np.isclose(K*(1 - np.exp(-b*dt)), beta_ode*dt)
