"""
Experiment 29, v4 environment: an independent cohort-component projection, outside Starsim, with exact aging, as a third reference for the age shares
of agents (P2b) and the ODE (R4x, R8x). See SUMMARY.md, Finding 5. Run from the repo root after run.py:

    PYTHONPATH=v4 v4/.venv/bin/python experiments/29_pregnancy_timing/cohort.py

The projection: single months of age (0 to 100 years, the last open) by sex, monthly steps from 1950 to 2051.
- Start: WPP's 1950 population by 5-year group and sex, spread evenly over each group's months.
- Each step: deaths at WPP's rate for each age group and sex (interpolated by year, at the step's middle); births at WPP's rate for each woman's
  current age group (as the ODE's births), 48.63% girls (WPP's female share in 1950, as P2b's and the ODE's with WPP's start); then everyone is one
  month older and the newborns are aged 0.
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
spec = importlib.util.spec_from_file_location('run29', here / 'run.py')
r29 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(r29)
r28, wpp = r29.r28, r29.wpp
sexes = r28.sexes
dt = 1/12
n_ages = 1200 + 1 # Months 0 to 1199, and 100+ (open)
age_mid = (np.arange(n_ages) + 0.5)/12 # Age at the middle of each month of age (years); the open group is taken as 100+
age_lo = np.arange(n_ages)/12


def start_population():
    pyr = ss.data.population('India', year=r29.start)
    lows = np.array(ss.data.age_edges(pyr.index), dtype=float)
    pop = {}
    for s in sexes:
        v = np.zeros(n_ages)
        for j, lo in enumerate(lows):
            hi = lows[j + 1] if j + 1 < len(lows) else None
            if lo >= 100:
                v[-1] += pyr[s].values[j]
            else:
                m = (age_lo >= lo) & (age_lo < (hi if hi is not None else 100))
                v[m] += pyr[s].values[j]/m.sum()
        pop[s] = v
    total = sum(v.sum() for v in pop.values())
    return {s: v/total for s, v in pop.items()}


def project(p_female=None):
    p_female = r29.wpp_female_share()[2] if p_female is None else p_female
    pop = start_population()
    mgroup = np.clip(np.digitize(age_lo, r28.edges) - 1, 0, len(r28.edges) - 1) # WPP's mortality group of each month of age
    fa = np.array(wpp.fert_ages, dtype=float)
    fgroup = np.clip(np.digitize(age_lo, fa) - 1, 0, len(fa) - 1)
    fertile = (age_lo >= 10) & (age_lo < 55)
    years = np.arange(r29.start, r29.stop)
    out = sc.objdict(years=years, pyr={s: [] for s in sexes}, total=[])
    t = float(r29.start)
    n_steps = int(round((r29.stop - r29.start)/dt))
    for i in range(n_steps + 1):
        if abs(t - round(t)) < 1e-6 and int(round(t)) in years: # 1 January: record by WPP's 5-year groups
            for s in sexes:
                out.pyr[s].append(r28.regroup(pop[s][None, :], np.minimum(age_lo, 100), r28.pyr_edges)[0])
            out.total.append(sum(v.sum() for v in pop.values()))
        if i == n_steps:
            break
        tm = t + dt/2
        mx = wpp.mx_at(tm) # [sex, mortality group]
        asfr = wpp.asfr_at(tm)
        births = float((pop['female']*np.where(fertile, asfr[fgroup], 0.0)).sum()*dt) # From women at the step's start
        for k, s in enumerate(sexes):
            v = pop[s]*np.exp(-mx[k][mgroup]*dt) # Deaths over the step
            new = np.zeros(n_ages)
            new[1:] = v[:-1]
            new[-1] += v[-1] # The open group keeps its own
            new[0] = births*(p_female if s == 'female' else 1 - p_female)
            pop[s] = new
        t = r29.start + (i + 1)*dt
    out.pyr = {s: np.array(v) for s, v in out.pyr.items()}
    out.total = np.array(out.total)
    return out


def shares(pyr):
    both = pyr['female'] + pyr['male']
    return both/both.sum(axis=1, keepdims=True)


def main():
    T = sc.timer()
    cc = project()
    data = sc.load(r29.fn)
    yrs = list(cc.years)
    def model_shares(o):
        idx = [list(o.years).index(y) for y in yrs]
        return shares({s: o.pyr[s][idx] for s in sexes})
    cs = shares(cc.pyr)
    ag = np.mean([model_shares(o) for o in data.P2b], axis=0)
    r8, r4 = model_shares(data.R8x[0]), model_shares(data.R4x[0])
    ex = 2*r8 - r4
    pop = dict(cc=cc.total/cc.total[0], P2b=r29.by_year(data.P2b, 'pop', yrs).mean(0), R8x=r29.by_year(data.R8x, 'pop', yrs)[0])
    pop['exact'] = 2*pop['R8x'] - r29.by_year(data.R4x, 'pop', yrs)[0]
    sc.save(here / 'cohort.obj', sc.objdict(cc=cc, cs=cs, agents=ag, R8x=r8, R4x=r4, exact=ex, pop=pop, years=yrs))
    lows = np.array(r28.pyr_edges, dtype=float)
    top = lows < 85
    lines = ['Cohort-component projection (exact aging, monthly) as the reference. Population ÷ 1950, ratio to it:']
    for y in [1960, 1980, 2000, 2020, 2050]:
        i = yrs.index(y)
        lines.append(f'  {y}: agents (P2b) {pop["P2b"][i]/pop["cc"][i]:.4f}, ODE R8x {pop["R8x"][i]/pop["cc"][i]:.4f}, ODE exact aging {pop["exact"][i]/pop["cc"][i]:.4f}')
    lines.append('Share in each 5-year group (0–84), ratio to the projection: largest gap over 1955–2050, and its age group and year')
    sel = [yrs.index(y) for y in yrs if y >= 1955]
    for name, m in [('agents (P2b)', ag), ('ODE R8x', r8), ('ODE exact aging', ex)]:
        rr = m[sel][:, top]/cs[sel][:, top]
        i, j = np.unravel_index(np.argmax(np.abs(rr - 1)), rr.shape)
        band = {int(lo): f'{rr[:, k].min():.3f}–{rr[:, k].max():.3f}' for k, lo in enumerate(lows[top]) if lo in [0, 5, 15, 35, 80]}
        lines.append(f'  {name}: largest {rr[i, j]:.3f} (ages {lows[top][j]:.0f}–{lows[top][j] + 4:.0f}, {yrs[sel[i]]}); ranges at 0–4, 5–9, 15–19, 35–39, 80–84: {band}')
    txt = '\n'.join(lines)
    print(txt)
    with open(here / 'cohort.log', 'w') as fh:
        fh.write(txt + '\n')
    fig, axs = plt.subplots(1, 3, figsize=(20, 5.5), sharey=True)
    for ax, (name, m) in zip(axs, [('Agents (P2b, mean of 10 reps)', ag), ('ODE, eighth-year ages (R8x)', r8), ('ODE, exact aging (extrapolated)', ex)]):
        im = ax.pcolormesh(yrs, lows[top] + 2.5, (m[:, top]/cs[:, top]).T, cmap='RdBu_r', vmin=0.95, vmax=1.05, shading='nearest')
        ax.set_title(f'{name} ÷ cohort-component projection', fontsize=10); ax.set_xlabel('Year')
    axs[0].set_ylabel('Age group (5-year)')
    fig.colorbar(im, ax=axs, label='Share of population in each 5-year group, ÷ the projection')
    fig.suptitle('Age shares against an independent cohort-component projection with exact aging (WPP rates, 48.63% girls at birth)', fontsize=11)
    fig.savefig(here / 'cohort_check.png', dpi=110, bbox_inches='tight'); plt.close(fig)
    print(f'cohort.py finished in {T.toc(output=True):.0f} s')
    return


if __name__ == '__main__':
    main()
