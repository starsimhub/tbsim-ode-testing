"""
Experiment 29, v4 environment: do agents' fixed-age bands in 5-year age shares (SUMMARY.md, Finding 5) come from ss.Pregnancy or from the agent
solver? Agents with births from ss.Births at a constant crude rate (B0), against the cohort-component projection with the same births. Run from the
repo root after run.py and cohort.py:

    PYTHONPATH=v4 v4/.venv/bin/python experiments/29_pregnancy_timing/births_test.py

- B0: 10 reps of 100,000 agents from WPP's 1950 population; births at 40 per 1,000 per year (ss.Births), 50% girls; exp 28's deaths by age and sex;
  monthly steps, 1950–2050.
- The reference: cohort.py's projection with the same crude birth rate and share of girls.
"""
import os
os.environ.setdefault('NUMBA_NUM_THREADS', '1')
import importlib.util
import numpy as np
import sciris as sc
import starsim as ss
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

here = sc.thispath(__file__)
def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod
r29 = load('run29', here / 'run.py')
cc = load('cohort29', here / 'cohort.py')
r28 = r29.r28
cbr = 0.040 # Crude birth rate, per person per year
const_mx = float(os.environ.get('EXP29_CONST_MX', 0)) # B1: one death rate at every age and sex (per year) instead of WPP's, in both models
tag = '_const' if const_mx else ''
reps, n_agents = r29.reps, r29.n_agents


def run_b0(rep):
    deaths = ss.death(rate=ss.peryear(const_mx)) if const_mx else r28.death_rate()
    sim = ss.Sim(diseases=ss.SIS(beta=0, init_prev=0), demographics=[ss.Births(birth_rate=ss.peryear(cbr)), deaths],
                 people=ss.population('India', year=r29.start, n_agents=n_agents), total_pop=n_agents, start=r29.start, stop=r29.stop,
                 dt=ss.months(1), solver='abm', rand_seed=29900 + rep, verbose=0)
    sim.stratify(age=r28.edges, sex=r28.sexes)
    sim.run()
    r = sim.results
    yv = np.asarray(sim.t.yearvec, dtype=float)
    years = np.arange(r29.start, int(np.floor(yv[-1])) + 1)
    idx = np.array([r28.year_index(sim, y) for y in years])
    pyr = {}
    for s in r28.sexes:
        v, lows = r28.by_age_sex(r.sis.S, s)
        pyr[s] = r28.regroup(v, lows, r28.pyr_edges)[idx]
    ppl = sim.people
    u = ppl.auids
    fem = float(np.asarray(ppl.female.raw[u])[np.asarray(ppl.uid.raw[u]) >= n_agents].mean()) # Girls among those born in the run, alive at the end
    return sc.objdict(years=years, pyr=pyr, rep=rep, fem_born=fem)


def project_cbr():
    """ cohort.py's projection with births at the crude rate cbr (× everyone), 50% girls """
    orig = cc.wpp.asfr_at
    pop0 = cc.start_population()
    out = sc.objdict(years=np.arange(r29.start, r29.stop), pyr={s: [] for s in r28.sexes})
    mgroup = np.clip(np.digitize(cc.age_lo, r28.edges) - 1, 0, len(r28.edges) - 1)
    pop = pop0
    t = float(r29.start)
    n_steps = int(round((r29.stop - r29.start)/cc.dt))
    for i in range(n_steps + 1):
        if abs(t - round(t)) < 1e-6 and int(round(t)) in out.years:
            for s in r28.sexes:
                out.pyr[s].append(r28.regroup(pop[s][None, :], np.minimum(cc.age_lo, 100), r28.pyr_edges)[0])
        if i == n_steps:
            break
        mx = np.full((2, len(r28.edges)), const_mx) if const_mx else cc.wpp.mx_at(t + cc.dt/2)
        births = cbr*sum(v.sum() for v in pop.values())*cc.dt
        for k, s in enumerate(r28.sexes):
            v = pop[s]*np.exp(-mx[k][mgroup]*cc.dt)
            new = np.zeros(cc.n_ages); new[1:] = v[:-1]; new[-1] += v[-1]
            new[0] = births*0.5
            pop[s] = new
        t = r29.start + (i + 1)*cc.dt
    out.pyr = {s: np.array(v) for s, v in out.pyr.items()}
    return out


def main():
    T = sc.timer()
    runs = sc.parallelize(run_b0, iterkwargs=[dict(rep=r) for r in range(reps)], ncpus=r29.ncpus)
    ref = project_cbr()
    yrs = list(ref.years)
    sh = lambda pyr, idx: (lambda b: b/b.sum(1, keepdims=True))((pyr['female'] + pyr['male'])[idx])
    ag = np.array([sh(o.pyr, [list(o.years).index(y) for y in yrs]) for o in runs])
    cs = sh(ref.pyr, slice(None))
    ratio = ag.mean(0)/cs
    se = (ag/cs).std(0, ddof=1)/np.sqrt(len(ag))
    lows = np.array(r28.pyr_edges, dtype=float)
    top = lows < 85
    sc.save(here / f'births_test{tag}.obj', sc.objdict(runs=runs, ref=ref, ratio=ratio, se=se, years=yrs))
    sel = [yrs.index(y) for y in yrs if y >= 1955]
    lines = [f'{"B1" if const_mx else "B0"} (agents, ss.Births at {1000*cbr:.0f} per 1,000, deaths {const_mx if const_mx else "WPP by age and sex"}, 10 reps) ÷ the cohort-component projection with the same births, share in each 5-year group, 1955–2050',
             f'Girls among those born in the run (alive at the end): {np.mean([o.fem_born for o in runs]):.4f}']
    for j, lo in enumerate(lows[top]):
        rr = ratio[sel, j]
        lines.append(f'  {lo:3.0f}–{lo + 4:3.0f}: {rr.min():.3f}–{rr.max():.3f} (mean {rr.mean():.4f}; SE about {se[sel, j].mean():.3f})')
    txt = '\n'.join(lines)
    print(txt)
    with open(here / f'births_test{tag}.log', 'w') as fh:
        fh.write(txt + '\n')
    fig, ax = plt.subplots(figsize=(9, 5.5))
    im = ax.pcolormesh(yrs, lows[top] + 2.5, ratio[:, top].T, cmap='RdBu_r', vmin=0.95, vmax=1.05, shading='nearest')
    fig.colorbar(im, ax=ax, label='Share in each 5-year group, agents ÷ the projection')
    ax.set_title(f'Agents with ss.Births ({1000*cbr:.0f} per 1,000; mean of {reps} reps) ÷ the cohort-component projection with the same births', fontsize=9)
    ax.set_xlabel('Year'); ax.set_ylabel('Age group (5-year)')
    fig.tight_layout(); fig.savefig(here / f'births_test{tag}.png', dpi=110); plt.close(fig)
    print(f'births_test.py finished in {T.toc(output=True):.0f} s')
    return


if __name__ == '__main__':
    main()
