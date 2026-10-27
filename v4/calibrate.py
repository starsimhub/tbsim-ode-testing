"""
Calibrate the ODE beta of each arm to 500/100k active prevalence at equilibrium, and record each arm's cohort progression
"""
import sciris as sc
import common as c

arms = sc.objdict(
    A0     = dict(parset='A', model='tb'),                     # tbsim defaults, constant progression (exps 01–07)
    B0     = dict(parset='B', model='tb'),                     # tb_LAI_TPT natural history, constant progression
    Ak     = dict(parset='A', model='classes'),                # A + time-since-infection progression (k = 5.846), exact latent classes
    Bk     = dict(parset='B', model='classes'),                # B + the same: the reference set of exps 08–11
    Bremote = dict(parset='B', model='tb', k_latent=c.k_tsi, inf_non_remote=0.0), # B with v4 tbsim.TB's own front-loading: a remote latent stage
)

if __name__ == '__main__':
    out = {}
    for key, arm in arms.items():
        cal = c.get_calib(key, **arm)
        prog = c.progression(**arm)
        out[key] = sc.mergedicts(dict(beta=cal.beta, prevalence=c.prevalence(cal.equilibrium)), prog)
        print(f'{key:8s} beta = {cal.beta:7.3f}, prevalence = {1e5*c.prevalence(cal.equilibrium):.2f}/100k, ProgressionY1 = {prog.ProgressionY1:.4f}, ProgressionRatio = {prog.ProgressionRatio:.3f}')
    sc.savejson(c.folder / 'results' / 'progression.json', out, indent=2)
