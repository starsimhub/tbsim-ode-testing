"""
Exp 16, tb_LAI_TPT as run: its no-reinfection progression cohort (cohort_progression.py) on v0.10's tbsim.TB at monthly dt, for the 12
face-validity parameter sets, with an analyzer that records the outputs compared in compare.py. See README.md.

The sim is built exactly as cohort_progression._run_once builds it (imported from ../tb_LAI_TPT): beta = 0, everyone in INFECTION at the
start, no demographics or networks, sym_dead = 0, k_asy = k_non = k. Classification follows tb_LAI_TPT's trackers: end-of-step state
snapshots, so an agent that passes INFECTION -> NON -> ASY within one step counts as leaving INFECTION to ASY.

v0.10 environment (explicit exception to the no-v0.10-arms rule, user 2026-11-04). Run from the repo root:
    .venv/bin/python experiments/16_two_state_vs_lai_tpt/run_lai_tpt.py
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
lai_dir = outdir.parents[2] / 'tb_LAI_TPT'
sys.path.insert(0, str(lai_dir))
import cohort_progression as cpr # noqa: E402 (tb_LAI_TPT's own module)

samples_csv = lai_dir / 'results/hm_real_20260728/trajectory_selection/collapsed_samples_eps1.2.csv'
n_sets = 12                 # tb_LAI_TPT's face-validity parameter sets: the first 12 rows
n_agents = 100_000
seeds = [0, 1, 2, 3]
years = 50                  # cohort_progression runs to the end of 2000 + years, i.e. years + 1 years of monthly steps
check_kw = dict(n_agents=20_000, seeds=(0, 1), years=12) # tb_LAI_TPT's own settings, for the check against its FirstPassage


class CohortOutcomes(ss.Analyzer):
    """
    Per-step cohort outputs, recorded after each step's transitions (index i is the end of month i + 1, as in cohort_progression.FirstPassage):
    - A, B, C: shares that have ever reached NON/ASY/SYM, ASY/SYM, SYM (tb_LAI_TPT's curves A, B, C, for a cohort)
    - latent: share in INFECTION; exit_cle, exit_non, exit_asy: cumulative shares whose first exit from INFECTION was to each
      (SYM in the snapshot counts as ASY, which it passed through)
    - prog_from_latent, latent_start: per step, agents leaving INFECTION to NON/ASY/SYM, and agents in INFECTION at the start of the step
    - first_B: each agent's step of first reaching ASY/SYM (-1: never); latent_1y, non_1y: agents in INFECTION, or in NON and not yet
      ASY, at the end of year 1
    """

    def init_post(self):
        super().init_post()
        n = len(self.sim.people)
        npts = self.sim.t.npts
        self.tvec_years = (np.arange(npts) + 1)*self.sim.t.dt_year
        self.i_1y = int(np.argmin(np.abs(self.tvec_years - 1.0)))
        self.everA = np.zeros(n, bool)
        self.everB = np.zeros(n, bool)
        self.everC = np.zeros(n, bool)
        self.exited = np.zeros(n, bool)
        self.prev = np.asarray(self.sim.diseases.tb.state.values).copy()
        self.first_B = np.full(n, -1)
        self.out = sc.objdict({k: np.zeros(npts) for k in ['A', 'B', 'C', 'latent', 'exit_cle', 'exit_non', 'exit_asy', 'prog_from_latent', 'latent_start']})
        self.exit_counts = dict(cle=0, non=0, asy=0)
        self.latent_1y = self.non_1y = None
        return

    def step(self):
        st = np.asarray(self.sim.diseases.tb.state.values)
        ti = self.sim.ti
        n = len(st)
        was_latent = self.prev == int(TBS.INFECTION)
        left = was_latent & (st != int(TBS.INFECTION)) & ~self.exited
        to_cle = left & (st == int(TBS.CLEARED))
        to_non = left & (st == int(TBS.NON_INFECTIOUS))
        to_asy = left & ((st == int(TBS.ASYMPTOMATIC)) | (st == int(TBS.SYMPTOMATIC)))
        self.exited |= left
        self.exit_counts['cle'] += to_cle.sum()
        self.exit_counts['non'] += to_non.sum()
        self.exit_counts['asy'] += to_asy.sum()
        newB = ~self.everB & ((st == int(TBS.ASYMPTOMATIC)) | (st == int(TBS.SYMPTOMATIC)))
        self.first_B[newB] = ti
        self.everA |= (st == int(TBS.NON_INFECTIOUS)) | (st == int(TBS.ASYMPTOMATIC)) | (st == int(TBS.SYMPTOMATIC))
        self.everB |= newB
        self.everC |= st == int(TBS.SYMPTOMATIC)
        o = self.out
        o.A[ti], o.B[ti], o.C[ti] = self.everA.mean(), self.everB.mean(), self.everC.mean()
        o.latent[ti] = (st == int(TBS.INFECTION)).mean()
        o.exit_cle[ti], o.exit_non[ti], o.exit_asy[ti] = (self.exit_counts[k]/n for k in ['cle', 'non', 'asy'])
        o.prog_from_latent[ti] = (to_non | to_asy).sum()
        o.latent_start[ti] = was_latent.sum()
        if ti == self.i_1y:
            self.latent_1y = st == int(TBS.INFECTION)
            self.non_1y = (st == int(TBS.NON_INFECTIOUS)) & ~self.everB
        self.prev = st.copy()
        return

    def later_split(self):
        """ Progression to ASY between years 1 and 10 (share of the cohort), and the share of it from agents latent at year 1 """
        i10 = int(np.argmin(np.abs(self.tvec_years - 10.0)))
        later = (self.first_B > self.i_1y) & (self.first_B <= i10)
        n = len(self.first_B)
        return later.sum()/n, (later & self.latent_1y).sum()/max(later.sum(), 1), (later & self.non_1y).sum()/max(later.sum(), 1)


def run_one(row, n_agents, seed, years):
    """ One cohort sim, built as cohort_progression._run_once builds it, with FirstPassage (for the check) and CohortOutcomes """
    cascade = {k: float(row[k]) for k in cpr.CASCADE_KEYS}
    k = float(row['k_asy'])
    tb = tbsim.TB(beta=ss.peryear(0.0), init_prev=ss.bernoulli(1.0), k_asy=k, k_non=k, **cpr._rates(cascade, float(row['inf_non']), float(row['inf_asy'])))
    fp, co = cpr.FirstPassage(), CohortOutcomes()
    sim = tbsim.Sim(n_agents=n_agents, diseases=tb, networks=[], demographics=[], analyzers=[fp, co], dt=ss.months(1), start=ss.date('2000-01-01'),
                    stop=ss.date(f'{2000 + int(years)}-12-31'), rand_seed=int(seed), verbose=0)
    sim.run()
    fp, co = sim.analyzers[0], sim.analyzers[1]
    later, from_latent, from_non = co.later_split()
    return sc.objdict(t=co.tvec_years, B_fp=fp.B, later=later, share_from_latent=from_latent, share_from_non=from_non, **co.out)


def run_set(i, row, n_agents, seed, years):
    out = run_one(row, n_agents, seed, years)
    out.update(param=i, seed=seed)
    return out


if __name__ == '__main__':
    post = pd.read_csv(samples_csv)
    rows = post.head(n_sets)
    w = rows.sample_weight.values/rows.sample_weight.sum()

    # Check: at tb_LAI_TPT's settings, our analyzer's B equals its FirstPassage's, and the targets equal the calibration's saved values
    print('Check against tb_LAI_TPT (param00, 20k agents, seeds 0 and 1):')
    chk = [run_one(rows.iloc[0], check_kw['n_agents'], s, check_kw['years']) for s in check_kw['seeds']]
    B_ours, B_fp = np.mean([c.B for c in chk], axis=0), np.mean([c.B_fp for c in chk], axis=0)
    t = chk[0].t
    b1, b10 = cpr.curves_at_years(t, B_ours, [1, 10])
    print(f'  max |our B - FirstPassage B| = {np.max(np.abs(B_ours - B_fp)):.2e}; ProgressionY1 {b1:.5f}, ProgressionRatio {b1/b10:.5f} (calibration: 0.05390, 0.63825)')
    ref = cpr.progression_targets(rows.iloc[0].to_dict())
    print(f'  cohort_progression.progression_targets: {ref}')

    jobs = [dict(i=i, row=rows.iloc[i], n_agents=n_agents, seed=s, years=years) for i in range(n_sets) for s in seeds]
    T = sc.timer()
    res = sc.parallelize(run_set, iterkwargs=jobs, ncpus=tbode.ncpus)
    T.toc(f'{len(jobs)} cohort sims ({n_agents:,} agents, {years + 1} years, monthly)')
    out = sc.objdict(settings=dict(samples_csv=str(samples_csv), n_sets=n_sets, n_agents=n_agents, seeds=seeds, years=years, dt='monthly'),
                     weights=w, params=rows.reset_index(drop=True), runs=res,
                     check=dict(max_abs_diff_B=float(np.max(np.abs(B_ours - B_fp))), Y1=float(b1), Ratio=float(b1/b10), progression_targets=ref))
    sc.save(outdir / 'lai_tpt.obj', out)
    tbode.save_versions(outdir, 'v0.10')
    os.replace(outdir / 'versions.json', outdir / 'versions_v010.json')

    print('Per parameter set (mean over seeds): ProgressionY1, ProgressionRatio, B(50), share of years 1-10 progression from latent at year 1')
    for i in range(n_sets):
        r = [x for x in res if x.param == i]
        B = np.mean([x.B for x in r], axis=0)
        b1, b10 = cpr.curves_at_years(r[0].t, B, [1, 10])
        print(f'  param{i:02d} (weight {w[i]:.3f}): Y1 {b1:.4f}, Ratio {b1/b10:.4f}, B(50) {B[-1]:.4f}, from latent {np.mean([x.share_from_latent for x in r]):.3f}')
