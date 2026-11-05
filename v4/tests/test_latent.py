"""
Checks of tbode4.latent: the exact Poisson-mixture classes against direct integration of the time-varying hazard r·exp(−k·τ), and the
two-state chain against exp 14 (set B). Run from the repo root: v4/.venv/bin/python -m pytest v4/tests/test_latent.py
"""
import numpy as np
from scipy.integrate import solve_ivp
import tbode4

# Set B: tb_LAI_TPT's best fit (row 0 of collapsed_samples_eps1.2.csv)
B = dict(inf_cle=0.5365965386150265, inf_non=0.4103330692741019, inf_asy=0.2741560588304815, non_rec=0.1706711906227213,
         non_asy=0.2426759489829796, pi=0.2415635040650483, k=5.845726750296858)


def tv_hazard_curve(p, years):
    """ Ever reached asymptomatic TB, by integrating the time-varying hazard directly: states L, NON, ASY """
    def rhs(t, y):
        L, N, _ = y
        h = np.exp(-p['k']*t)
        return [-(p['inf_cle'] + (p['inf_non'] + p['inf_asy'])*h)*L,
                p['inf_non']*h*L - (p['non_asy'] + p['non_rec'])*N,
                p['inf_asy']*h*L + p['non_asy']*N]
    sol = solve_ivp(rhs, (0, max(years)), [1, 0, 0], t_eval=years, rtol=1e-11, atol=1e-13, method='LSODA')
    return sol.y[2]


def test_exact_classes():
    """ The classes reproduce the time-varying hazard model exactly, for set B and for a slowly declining hazard """
    years = [0.25, 1, 2, 5, 10, 30]
    for k in [B['k'], 0.3]:
        p = dict(B, k=k)
        ref = tv_hazard_curve(p, years)
        got = tbode4.ever_asy(p, 'exact', years=years)
        assert np.allclose(got, ref, rtol=1e-7, atol=1e-10), f'k={k}: {got} vs {ref}'
    return


def test_latent_hazard():
    """ The hazard among people still latent equals r0·exp(−k·τ) for the classes """
    t, X, idx = tbode4.cohort_curve(B, 'exact', tmax=5)
    h = tbode4.latent_hazard(X, idx)
    ref = (B['inf_non'] + B['inf_asy'])*np.exp(-B['k']*t)
    assert np.allclose(h, ref, rtol=1e-6, atol=1e-12)
    return


def test_two_state_set_b():
    """ With no remote progression and k_latent = 5.5333 (exp 14), ProgressionY1 matches the exact model's to 1e-4 """
    exact = tbode4.progression_targets(B, 'exact')
    two = tbode4.progression_targets(dict(B, k_latent=5.5333), 'two-state')
    assert abs(two['ProgressionY1']/exact['ProgressionY1'] - 1) < 1e-4
    assert abs(exact['ProgressionY1'] - 0.0511) < 1e-4 and abs(exact['ProgressionRatio'] - 0.6745) < 1e-4
    return


def test_curves_and_exits():
    """ With the cascade, curve B equals ever_asy without it; A >= B >= C; and the first exits from latent infection add up to 1 """
    p = dict(B, asy_non=1.7385946555, asy_sym=0.7435451323, sym_asy=0.6630323475)
    two = dict(p, k_latent=3.0, inf_non_remote=0.01*p['inf_non'], inf_asy_remote=0.01*p['inf_asy'])
    for structure, q in [('exact', p), ('two-state', two)]:
        t, A = tbode4.first_passage(q, structure, 'A', tmax=20)
        _, Bc = tbode4.first_passage(q, structure, 'B', tmax=20)
        _, C = tbode4.first_passage(q, structure, 'C', tmax=20)
        assert np.allclose(Bc[[48, 480]], tbode4.ever_asy(q, structure, years=(1, 10)), rtol=1e-9)
        assert np.all(A >= Bc - 1e-12) and np.all(Bc >= C - 1e-12)
        _, ex = tbode4.latent_exits(q, structure, tmax=20)
        assert np.allclose(ex.cleared + ex.non + ex.asy + ex.still_latent, 1, atol=1e-10)
        assert np.allclose(ex.non + ex.asy, A, atol=1e-10) # First progression = first exit to NON or ASY
    return


if __name__ == '__main__':
    test_exact_classes()
    test_latent_hazard()
    test_two_state_set_b()
    test_curves_and_exits()
    print('All passed')
