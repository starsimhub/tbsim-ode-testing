"""
Experiment 17, question 2, tb_LAI_TPT as run: v0.10's tbsim.TB at monthly dt, in a cohort at a constant force of infection. See README.md.

The cohort is built as tb_LAI_TPT's cohort_progression._run_once builds it (beta = 0, everyone in INFECTION at the start, no
demographics or networks, sym_dead = 0, k_asy = k_non = k, each set's calibrated rates), plus each set's reinfection risk ratio pi for
everyone no longer infected and for latent people. `TBConstantFOI` infects at a constant force of infection lam, in place of
transmission (at the same point in the step), through tbsim.TB's own set_prognoses: so reinfecting a latent agent resets its clock,
exactly as transmission does. Non-infectious TB isn't reinfected (in v0.10 that changes nothing).

The per-person outputs are defined once, in cohort.py, for both environments.

v0.10 environment (tb_LAI_TPT as run: allowed for exp 17, user 2026-11-05). Run from the repo root:

    .venv/bin/python experiments/17_latent_reinfection/run_lai_tpt.py check   # lam = 0 against cohort_progression; timing
    .venv/bin/python experiments/17_latent_reinfection/run_lai_tpt.py         # all runs -> lai_tpt.obj, versions_v010.json

Rerunning resumes: runs already saved are kept.
"""
import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss
import tbsim
from tbsim.tb import TBS
import tbode

outdir = Path(__file__).resolve().parent
sys.path.insert(0, str(outdir))
import cohort as co # noqa: E402
lai_dir = outdir.parents[2] / 'tb_LAI_TPT'
sys.path.insert(0, str(lai_dir))
import cohort_progression as cpr # noqa: E402 (tb_LAI_TPT's own module)

samples_csv = lai_dir / 'results/hm_real_20260728/trajectory_selection/collapsed_samples_eps1.2.csv'
n_sets = 12                 # tb_LAI_TPT's face-validity parameter sets: the first 12 rows
n_agents = 250_000
seeds = [0, 1, 2, 3]
years = 30
lams = [0.0, 0.10, 0.60]    # Force of infection, per year
fn = outdir / 'lai_tpt.obj'


class TBConstantFOI(tbsim.TB):
    """
    tbsim.TB with infection at a constant force of infection `lam` (per year) times each agent's relative risk, in place of
    transmission; and every move in step_transitions recorded, so the per-person outputs see each one (agents can move more than once
    per step).

    Args:
        lam (float): force of infection, per year
        latent_reinf (bool): whether latent (INFECTION) agents are reinfected
        kwargs (dict): passed to tbsim.TB (with beta = 0)
    """
    def __init__(self, lam=0.0, latent_reinf=True, **kwargs):
        super().__init__(name='tb', **kwargs)
        self.lam = lam
        self.latent_reinf = latent_reinf
        self._rng_foi = ss.random(name='tb_rng_foi')
        self.events = None
        return

    def step(self):
        """ Constant-force infection where transmission would be, then tbsim.TB's step (whose transmission does nothing: beta = 0) """
        self.events = {k: [] for k in co.event_names}
        if self.lam > 0:
            uid = np.asarray(self.sim.people.auids)
            st = np.asarray(self.state[ss.uids(uid)])
            latent = st == TBS.INFECTION
            elig = (st == TBS.SUSCEPTIBLE) | (st == TBS.CLEARED) | (latent & self.latent_reinf)
            u = ss.uids(uid[elig])
            rr = np.ones(len(u))
            rr[latent[elig]] = self.pars.rr_reinfection_inf
            cle = st[elig] == TBS.CLEARED
            rr[cle] = self.rr_reinfection[u[cle]]
            p = 1 - np.exp(-rr*self.lam*self.sim.t.dt_year)
            hit = self._rng_foi.rvs(u) < p
            u_hit = u[hit]
            was_latent = latent[elig][hit]
            self.events['reset'].append(np.asarray(u_hit[was_latent]))
            self.events['new'].append(np.asarray(u_hit[~was_latent]))
            self.set_prognoses(u_hit)
        super().step()
        return

    def transition(self, uids, to, rng):
        """ tbsim.TB's transition, recording each agent that moves """
        if len(uids) == 0:
            return
        before = np.asarray(self.state[uids]).copy()
        super().transition(uids, to, rng)
        after = np.asarray(self.state[uids])
        moved = before != after
        u, b, a = np.asarray(uids)[moved], before[moved], after[moved]
        if self.events is not None:
            from_latent = b == TBS.INFECTION
            self.events['progress'].append(u[from_latent & ((a == TBS.NON_INFECTIOUS) | (a == TBS.ASYMPTOMATIC))])
            self.events['clear'].append(u[from_latent & (a == TBS.CLEARED)])
            self.events['asy'].append(u[a == TBS.ASYMPTOMATIC])
            self.events['sym'].append(u[a == TBS.SYMPTOMATIC])
        return


class CohortTracker(ss.Analyzer):
    """ Feeds each step's events to cohort.Tracker (index i is the end of month i + 1, as cohort_progression.FirstPassage) """

    def init_post(self):
        super().init_post()
        self.tr = co.Tracker(len(self.sim.people), self.sim.t.npts, self.sim.t.dt_year, offset=1, infect_at_start=True) # TBConstantFOI infects before the transitions
        return

    def step(self):
        ev = self.sim.diseases.tb.events
        self.tr.step(self.sim.ti, {k: np.concatenate(v) if v else [] for k, v in ev.items()})
        return


def set_pars(row):
    """ One parameter set's TB parameters, as cohort_progression builds them, plus the reinfection risk ratio pi for everyone """
    cascade = {k: float(row[k]) for k in cpr.CASCADE_KEYS}
    k = float(row['k_asy'])
    pi = float(row['rr_reinfection_rec'])
    return dict(k_asy=k, k_non=k, rr_reinfection_cleared=pi, rr_reinfection_rec=pi, **cpr._rates(cascade, float(row['inf_non']), float(row['inf_asy'])))


def run_one(row, lam, on, seed, n=n_agents, nyears=years):
    """ One cohort: per-person outputs (cohort.Tracker) and FirstPassage's curve B """
    tb = TBConstantFOI(lam=lam, latent_reinf=on, beta=ss.peryear(0.0), init_prev=ss.bernoulli(1.0), **set_pars(row))
    sim = tbsim.Sim(n_agents=n, diseases=tb, networks=[], demographics=[], analyzers=[cpr.FirstPassage(), CohortTracker()], dt=ss.months(1),
                    start=ss.date('2000-01-01'), stop=ss.date(f'{2000 + nyears - 1}-12-31'), rand_seed=int(seed), verbose=0)
    T = sc.timer()
    sim.run()
    out = sim.analyzers[1].tr.results()
    out.B_fp = sim.analyzers[0].B
    out.time = T.toc(output=True)
    return out


def run_job(i, row, lam, on, seed):
    out = run_one(row, lam, on, seed)
    out.update(param=i, lam=lam, on=on, seed=seed)
    return out


def check():
    """ At lam = 0, the cohort equals tb_LAI_TPT's own (cohort_progression._run_once) draw for draw; and a timing run """
    rows = pd.read_csv(samples_csv).head(n_sets)
    row = rows.iloc[0]
    cascade = {k: float(row[k]) for k in cpr.CASCADE_KEYS}
    ref = cpr._run_once(cascade, float(row['inf_non']), float(row['inf_asy']), float(row['k_asy']), 20_000, 0, 12)
    ours = run_one(row, 0.0, True, 0, n=20_000, nyears=13)
    n = len(ref.B)
    print(f'lam = 0, param00, 20k agents, seed 0: max |FirstPassage B in our sim - cohort_progression B| = {np.max(np.abs(ours.B_fp[:n] - ref.B)):.2e}')
    print(f'  ProgressionY1: FirstPassage (end-of-step snapshots) {np.interp(1, ours.t, ours.B_fp):.5f} (calibration 0.05390); '
          f'counted from transitions {np.interp(1, ours.t, ours.B):.5f}')
    for lam in [0.10, 0.60]:
        o = run_one(row, lam, True, 0)
        print(f'lam = {lam}, on, {n_agents:,} agents: {o.time:.0f} s; B(10) {np.interp(10, o.t, o.B):.4f}, B(30) {o.B[-1]:.4f}, '
              f'resets per person {np.arange(len(o.resets_pre)) @ o.resets_pre/o.n:.3f}, episodes {np.arange(len(o.episodes)) @ o.episodes/o.n:.3f}')
    return


def main():
    post = pd.read_csv(samples_csv)
    rows = post.head(n_sets)
    data = sc.load(fn) if fn.exists() else sc.objdict(runs=sc.objdict())
    data.settings = dict(samples_csv=str(samples_csv), n_sets=n_sets, n_agents=n_agents, seeds=seeds, years=years, lams=lams, dt='monthly')
    data.weights = rows.sample_weight.values/rows.sample_weight.sum()
    data.params = rows.reset_index(drop=True)
    arms = [(lam, on) for lam in lams for on in ([True] if lam == 0 else [True, False])]
    jobs = [dict(i=i, row=rows.iloc[i], lam=lam, on=on, seed=s) for lam, on in arms for i in range(n_sets) for s in seeds]
    key = lambda j: f"p{j['i']}_lam{j['lam']}_{'on' if j['on'] else 'off'}_s{j['seed']}"
    todo = [j for j in jobs if key(j) not in data.runs]
    for chunk in [todo[i:i + 24] for i in range(0, len(todo), 24)]: # Saved as it goes
        T = sc.timer()
        res = sc.parallelize(run_job, iterkwargs=chunk, ncpus=tbode.ncpus)
        for j, r in zip(chunk, res):
            data.runs[key(j)] = r
        sc.save(fn, data)
        print(f'{len(data.runs)}/{len(jobs)} runs ({T.toc(output=True):.0f} s for {len(chunk)})', flush=True)
    tbode.save_versions(outdir, 'v0.10')
    os.replace(outdir / 'versions.json', outdir / 'versions_v010.json')
    return data


def perperson():
    """
    Seed 0 of each set's reinfection-on runs again, for each person's total time latent and number of latent spells (added to
    cohort.py after the main runs). The dynamics are unchanged, so these are the same runs as in lai_tpt.obj: checked on curve B.
    """
    fn_pp = outdir / 'lai_tpt_perperson.obj'
    rows = pd.read_csv(samples_csv).head(n_sets)
    saved = sc.load(fn)
    jobs = [dict(i=i, row=rows.iloc[i], lam=lam, on=True, seed=0) for lam in lams for i in range(n_sets)]
    T = sc.timer()
    res = sc.parallelize(run_job, iterkwargs=jobs, ncpus=tbode.ncpus)
    out = sc.objdict()
    for j, r in zip(jobs, res):
        key = f"p{j['i']}_lam{j['lam']}_on_s0"
        assert np.array_equal(r.B, saved.runs[key].B), f'{key} does not reproduce the saved run'
        out[key] = sc.objdict(time_latent=r.hists.time_latent, time_latent_open=r.hists.time_latent_open, n_spells=r.n_spells, n=r.n, dt_year=r.dt_year)
    sc.save(fn_pp, out)
    print(f'{len(jobs)} runs in {T.toc(output=True):.0f} s; all reproduce the saved runs (curve B identical)')
    return out


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    if which == 'check':
        check()
    elif which == 'all':
        main()
    elif which == 'perperson':
        perperson()
    else:
        raise ValueError(f'Unknown part "{which}"')
