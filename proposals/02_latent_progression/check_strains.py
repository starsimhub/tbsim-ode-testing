"""
How tbsim v4's two latent states (k_latent) handle latent reinfection once there are strains (drugs=). Numbers for options.md.

Question: when someone with a remote latent infection is reinfected with a strain they don't carry, does v4 add the new strain
(superinfection, and a reset to recent), replace it, or ignore it? The declared transitions don't say: `reinfect_remote` is
remote → recent at rr_reinfection_inf·foi, and superinfection is declared only for the recent state (`superinfection_infection`).

Check: an ODE run that starts with remote latent infections, half DS and half RIF, and active TB of both strains, with superinfection
of recent latent infections turned off (rr_superinfection_latent = 0) and no de-novo resistance (q_prog = 0). If anyone in the recent
state carries DS+RIF, it came through reinfect_remote, which then adds the new strain. Rerun with rr_superinfection_latent = 1 for
comparison.

Run from the repo root: v4/.venv/bin/python proposals/02_latent_progression/check_strains.py
"""
import json
from pathlib import Path
import starsim as ss
import tbsim
import tbode4

outdir = Path(__file__).resolve().parent


def run(rr_superinfection_latent):
    """ Two years as an ODE at monthly dt; the people in each latent state, by set of strains, at the end """
    tb = tbsim.TB(drugs=['RIF'], k_latent=ss.peryear(5.533), inf_non_remote=ss.peryear(0), inf_asy_remote=ss.peryear(0),
                  init=dict(infection_remote=0.4, asymptomatic=0.02, symptomatic=0.02), strains=dict(DS=0.5, RIF=0.5),
                  rr_superinfection_latent=rr_superinfection_latent, rr_reinfection_inf=1.0, q_prog=0, beta=20)
    sim = tbsim.Sim(tb=tb, solver='ode', start=2000, stop=2002, dt=ss.months(1)).run()
    res = sim.results.tb
    return {f'{state}[{level}]': float(res[state].select(resistance=level).values[-1])
            for state in ['infection', 'infection_remote'] for level in ['DS', 'RIF', 'DS+RIF']}


if __name__ == '__main__':
    out = {f'rr_superinfection_latent={rr}': run(rr) for rr in [0.0, 1.0]}
    for arm, vals in out.items():
        print(arm)
        for key, val in vals.items():
            print(f'    {key:26s} {val:10.1f}')
    (outdir / 'check_strains.json').write_text(json.dumps(out, indent=2) + '\n')
    tbode4.save_versions(outdir, 'v4')
    print(f'Recent DS+RIF with superinfection of recent latents off: {out["rr_superinfection_latent=0.0"]["infection[DS+RIF]"]:.1f} '
          '(> 0: reinfect_remote adds the new strain)')
