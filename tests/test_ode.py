"""
Regression checks for tbode.TBODE.
"""

import numpy as np
import tbsim.compartmental as tbc
import tbode


def test_matches_tb_ode():
    """With all extensions off, TBODE reproduces tbsim's TB_ODE."""
    ref = tbc.TB_ODE().run(start_time=1500, end_time=2020)
    res = tbode.TBODE(init_sym=1e3 / 1e5).run(start=1500, stop=2020, dt=1)
    for s in tbode.states:
        assert np.allclose(res[s], ref[s], rtol=1e-4, atol=1e-2), f'{s} differs from TB_ODE'


def test_equilibrium_and_calibration():
    """calibrate_beta hits the target, and the equilibrium is stationary when integrated forward."""
    target = 500e-5
    beta = tbode.calibrate_beta(target, replace_tb_deaths=False)
    ode = tbode.TBODE(beta=beta, replace_tb_deaths=False)
    x = ode.equilibrium()
    assert np.isclose(ode.outputs(x).prev_active[0], target, rtol=1e-6)
    res = ode.run(start=0, stop=50, dt=1, y0=x * 1e5)
    assert np.allclose(res.prev_active, target, rtol=1e-6)


def test_non0_is_subset():
    """The auxiliary never-yet-asymptomatic NON_INFECTIOUS compartment stays within NON_INFECTIOUS."""
    res = tbode.TBODE(beta=15).run(start=0, stop=100, dt=1)
    assert np.all(res.NON0 <= res.NON_INFECTIOUS + 1e-9) and np.all(res.NON0 >= -1e-9)


def test_beta_edge_round_trip():
    """The per-edge beta reproduces the ODE's per-step infection probability on a degree-K network."""
    beta_ode, K, dt = 12, 10, 1 / 12
    b = tbode.beta_edge_from_ode(beta_ode, K, dt)
    assert np.isclose(K * (1 - np.exp(-b * dt)), beta_ode * dt)
