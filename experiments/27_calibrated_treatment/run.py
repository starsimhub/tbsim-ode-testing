"""
Experiment 27, v4 environment: do v4's agents match its ODE with treatment, time-varying TB mortality and over-detection, at 12 sets
drawn from exp 24's posterior? See README.md. Run from the repo root, in order:

    PYTHONPATH=v4 v4/.venv/bin/python experiments/27_calibrated_treatment/run.py sets     # draw the 12 sets, against each seed's 95% intervals
    PYTHONPATH=v4 v4/.venv/bin/python experiments/27_calibrated_treatment/run.py check    # v4's ODE against exp 24's fastode, at each set
    PYTHONPATH=v4 v4/.venv/bin/python experiments/27_calibrated_treatment/run.py          # every run, then the summary
    PYTHONPATH=v4 v4/.venv/bin/python experiments/27_calibrated_treatment/run.py summary  # summaries and figures from the saved runs
"""
import os
os.environ.setdefault('NUMBA_NUM_THREADS', '1') # One thread per process (Starsim's kernels otherwise use every core); before importing starsim
import sys
import importlib.util
import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss
import matplotlib.pyplot as plt

outdir = sc.thispath(__file__)
root = outdir.parent.parent
sys.path.insert(0, str(root / 'v4'))
import tbode4 # noqa: E402
from tbode4 import fastode as fo # noqa: E402


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


exp24 = root / 'experiments' / '24_treatment_calibration'
sys.path.insert(0, str(exp24))
cb = load_module('calib24', exp24 / 'calib.py') # Exp 24's parameter mapping, schedules and targets
run24 = load_module('run24', exp24 / 'run.py') # Its target intervals

# Settings (README)
seeds24 = ['real', 'real_seed2'] # Exp 24's two IMIS seeds: 6 sets from each
sets_per_seed, draw_seed = 6, 27
reps, n_agents = 20, 200_000
start, stop = 1980, 2041
dts = sc.objdict(month=ss.months(1), week=ss.weeks(1))
windows = sc.objdict(early=(1995, 2010), late=(2010, 2025), const=(2025, 2041)) # 1995–2009, 2010–2024, 2025–2040
outputs = ['prevalence', 'prop_asymptomatic', 'incidence', 'notifications', 'on_treatment', 'mortality', 'mortality_tx', 'relapse',
           'recent', 'remote', 'cleared']
stocks = ['prevalence', 'prop_asymptomatic', 'on_treatment', 'recent', 'remote', 'cleared']
per_case = ['incidence', 'notifications', 'on_treatment', 'mortality', 'mortality_tx'] # Also compared per prevalent case
run_targets = ['prevalence', 'prop_asymptomatic', 'prev_to_inc', 'mort_to_inc', 'prev_to_notif', 'inc_ratio', 'mort_ratio', 'prev_2040_ratio']
fn = sc.objdict(sets=outdir / 'sets.obj', runs=outdir / 'runs.obj', ode=outdir / 'ode.obj')


# Sets and models

def draw_sets():
    """ 6 distinct draws from each seed's posterior (weighted resampling without replacement), with their untreated equilibria """
    rng = np.random.default_rng(draw_seed)
    sets, rows = [], []
    for seed24 in seeds24:
        res = sc.load(exp24 / f'{seed24}.obj')
        idx = rng.choice(len(res.weights), size=sets_per_seed, replace=False, p=res.weights)
        for i in idx:
            theta = np.asarray(res.theta[i], dtype=float)
            nat, trend = cb.to_pars(theta)
            sets.append(sc.objdict(seed24=seed24, post_index=int(i), weight=float(res.weights[i]), theta=dict(zip(cb.names, map(float, theta))),
                                   nat=nat, trend=trend, outs=cb.simulate(theta)))
            rows.append(dict(set=len(sets) - 1, seed24=seed24, **sets[-1].theta))
    sc.save(fn.sets, sets)
    df = pd.DataFrame(rows)
    with pd.option_context('display.width', 250, 'display.max_columns', 30):
        print('Drawn sets:'); print(df.round(3).to_string(index=False))
        for seed24 in seeds24:
            post = pd.read_csv(exp24 / f'{seed24}_posterior.csv').set_index('parameter')
            d = df[df.seed24 == seed24][cb.names]
            print(f'\n{seed24}: drawn range against the posterior 95% interval')
            print(pd.DataFrame({'min drawn': d.min(), 'max drawn': d.max(), 'post 2.5%': post.q025, 'post 97.5%': post.q975}).round(3).to_string())
    return sets


def starting_state(s):
    """ The untreated equilibrium (TBTreat's init=), the starting NON0 share (for incidence), and the share of NON already asymptomatic """
    p = fo.make_pars(**s.nat)
    d = dict(zip(fo.names(p), fo.equilibrium(p)))
    init = {k: d[k] for k in fo.state_names if k not in ['non0', 'non1']}
    init['non_infectious'] = d['non0'] + d['non1']
    init = {k: v for k, v in init.items() if v > 0}
    return init, d['non0'], d['non1']/(d['non0'] + d['non1'])


def make_sim(s, dt, solver, stages=None, rand_seed=0):
    """ TBTreat for one set: a fixed course, coverage, efficacy and sym_dead over time as ss.timeseries (exp 24's knots, linear) """
    init, non0, share_flagged = starting_state(s)
    t = s.trend
    tx = tbode4.lai_tpt_treatment(row=dict(tx_covg_2010=t.tx_covg_2010, tx_covg_2024=t.tx_covg_2024, efficacy_1995=t.efficacy_1995,
                                           efficacy_2024=t.efficacy_2024, p_relapse=t.p_relapse), step=dts[dt].years, end=stop)
    sym_dead = ss.peryear(ss.timeseries({1995.0: t.sym_dead_1995, 2024.0: t.sym_dead_2024}, interpolation='linear'))
    nat = {k: v for k, v in s.nat.items() if k != 'sym_dead'}
    tb = tbode4.TBTreat(init=init, course='fixed', stages=stages, tx_rate=tx.tx_rate, tx_success=tx.tx_success, tx_p_relapse=t.p_relapse,
                        dur_relapse_tx=1.5, sym_dead=sym_dead, **nat)
    if solver == 'abm':
        kw = dict(n_agents=n_agents, total_pop=n_agents, rand_seed=rand_seed, analyzers=[tbode4.Incidence(share_flagged=share_flagged), tbode4.TxChecks()])
    else:
        kw = dict(total_pop=1)
    sim = ss.Sim(tb, solver=solver, start=start, stop=stop, dt=dts[dt], verbose=0, **kw)
    return sim, non0


def per_step(sim, incidence):
    """ Each output per step, per person (per year for flows), with tbode4.tx_outputs; stocks at index i are at t_i, flows the step from t_i """
    out = tbode4.tx_outputs(sim, incidence)
    r = sim.results.tb
    alive = np.asarray(sim.results.alive.values, dtype=float)
    v = lambda k: np.asarray(r[k].values, dtype=float)
    out.prop_asymptomatic = v('asymptomatic')/(v('asymptomatic') + v('symptomatic'))
    out.recent, out.remote, out.cleared = v('infection')/alive, v('infection_remote')/alive, v('cleared')/alive
    out.tx_start = v('tx_start')/alive # Per person per step (for yearly totals)
    out.deaths = (v('tb_death') + v('tb_death_tx'))/alive
    out.first = np.asarray(incidence.first if hasattr(incidence, 'first') else incidence, dtype=float)/alive
    return out


def yearly(out, dt_year):
    """ Yearly means per 100,000 (shares for prop_asymptomatic): stocks at each step's end, flows over the step, dated by its middle """
    year = np.floor(out.year[1:] - dt_year/2 + 1e-9).astype(int)
    scale = lambda k: 1.0 if k == 'prop_asymptomatic' else 1e5
    df = pd.DataFrame({k: (out[k][1:] if k in stocks else out[k][:-1])*scale(k) for k in outputs})
    df['year'] = year
    return df.groupby('year').mean()


def target_values(out, over_detection):
    """ The run's values of exp 24's targets, as tb_LAI_TPT defines them: point values on 1 January, totals over each calendar year """
    year = out.year
    at = lambda y: int(np.argmin(np.abs(year - y))) # The step nearest 1 January (exact at monthly steps)
    total = lambda k, y: out[k][(year >= y - 1e-9) & (year < y + 1 - 1e-9)].sum() # Flows indexed by their step's start
    prev = lambda y: out.prevalence[at(y)]
    inc = {y: total('first', y) for y in [2015, 2020, 2024]}
    deaths = {y: total('deaths', y) for y in [2015, 2024]}
    notif = (1 + over_detection)*total('tx_start', 2020)
    return dict(prevalence=prev(2020), prop_asymptomatic=out.prop_asymptomatic[at(2020)], prev_to_inc=prev(2020)/inc[2020],
                mort_to_inc=deaths[2024]/inc[2024], prev_to_notif=prev(2020)/notif, inc_ratio=inc[2024]/inc[2015],
                mort_ratio=deaths[2024]/deaths[2015], prev_2040_ratio=prev(2040)/prev(2020))


def seed_of(i, dt, rep):
    """ New seeds for this experiment: 27000 + 100·(2·set + step) + rep """
    return 27000 + 100*(2*i + list(dts).index(dt)) + rep


def run_agents(i, dt, rep):
    s = sc.load(fn.sets)[i]
    sim, _ = make_sim(s, dt, 'abm', rand_seed=seed_of(i, dt, rep))
    sim.run()
    inc, chk = sim.analyzers[0], sim.analyzers[1]
    out = per_step(sim, inc)
    alive = np.asarray(sim.results.alive.values, dtype=float)
    return sc.objdict(set=i, dt=dt, rep=rep, seed=seed_of(i, dt, rep), yearly=yearly(out, sim.t.dt_year), targets=target_values(out, s.trend.over_detection),
                      checks=dict(chk.c), course_steps=dict(chk.course_steps), alive_range=(float(alive.min()), float(alive.max())), dt_year=sim.t.dt_year)


def run_ode(i, dt, stages=None):
    s = sc.load(fn.sets)[i]
    sim, non0 = make_sim(s, dt, 'ode', stages=stages)
    sim.run()
    out = per_step(sim, tbode4.ode_incidence(sim, non0=non0).first)
    return sc.objdict(set=i, dt=dt, stages=stages, yearly=yearly(out, sim.t.dt_year), targets=target_values(out, s.trend.over_detection))


def mapping_check():
    """ v4's ODE (12 Erlang stages, and the exact fixed course) against exp 24's fastode trajectory: prevalence on 1 January, 1995–2040 """
    sets = sc.load(fn.sets)
    years = np.arange(1995, 2041, 1.0)
    rows = []
    for i, s in enumerate(sets):
        r = cb.trajectories(np.array([s.theta[k] for k in cb.names]), t_eval=np.concatenate([[cb.t0], years]))
        fast = (r.X[1:, fo.A] + r.X[1:, fo.Y])
        for stages in [cb.tx_stages, None]:
            sim, _ = make_sim(s, 'month', 'ode', stages=stages)
            sim.run()
            yv = np.asarray(sim.t.yearvec, dtype=float)
            res = sim.results.tb
            prev = (np.asarray(res['asymptomatic'].values) + np.asarray(res['symptomatic'].values))/np.asarray(sim.results.alive.values)
            v4 = np.array([prev[int(np.argmin(np.abs(yv - y)))] for y in years])
            ratio = v4/fast
            rows.append(dict(set=i, course=f'{stages} stages' if stages else 'exact fixed', max_abs_dev=float(np.max(np.abs(ratio - 1))),
                             year_of_max=int(years[np.argmax(np.abs(ratio - 1))]), ratio_2020=float(ratio[years == 2020][0]), ratio_2040=float(ratio[-1])))
    df = pd.DataFrame(rows)
    df.to_csv(outdir / 'mapping.csv', index=False)
    with pd.option_context('display.width', 200, 'display.precision', 5):
        print("v4's ODE prevalence ÷ exp 24's fastode, monthly, 1995–2040 (1 January each year):")
        print(df.to_string(index=False))
    return df


def run_all():
    T = sc.timer()
    n_sets = len(sc.load(fn.sets))
    ode = [run_ode(i, dt) for i in range(n_sets) for dt in dts]
    sc.save(fn.ode, ode)
    print(f'ODE runs done ({T.toc(output=True):.0f} s)', flush=True)
    partdir = outdir / 'runs_by_group'
    partdir.mkdir(exist_ok=True)
    runs = []
    groups = [(i, dt) for dt in dts for i in range(n_sets)] # Monthly first
    for g, (i, dt) in enumerate(groups): # Saved group by group, so a stopped run resumes where it left off
        part = partdir / f'set{i:02d}_{dt}.obj'
        if part.exists():
            group = sc.load(part)
        else:
            group = sc.parallelize(run_agents, iterkwargs=[dict(i=i, dt=dt, rep=r) for r in range(reps)], ncpus=tbode4.ncpus)
            sc.save(part, group)
        runs += group
        print(f'Group {g + 1} of {len(groups)} (set {i}, {dt}): done ({T.toc(output=True):.0f} s)', flush=True)
    sc.save(fn.runs, runs)
    tbode4.save_versions(outdir, 'v4')
    print(f'{len(runs)} agent runs done ({T.toc(output=True):.0f} s)')
    return


# Summaries

def window_mean(df, w):
    lo, hi = windows[w]
    return df[(df.index >= lo) & (df.index < hi)].mean()


def by_set(runs, ode):
    """ Agents ÷ ODE for each set, step and window: each output, and the per-case outputs ÷ prevalence; mean ± SE across reps """
    rows = []
    for o in ode:
        rs = [r for r in runs if r.set == o.set and r.dt == o.dt]
        for w in windows:
            ow = window_mean(o.yearly, w)
            ratios = []
            for r in rs:
                aw = window_mean(r.yearly, w)
                d = {k: aw[k]/ow[k] for k in outputs}
                d.update({f'{k}_per_case': (aw[k]/aw.prevalence)/(ow[k]/ow.prevalence) for k in per_case})
                ratios.append(d)
            ratios = pd.DataFrame(ratios)
            for k in ratios:
                rows.append(dict(set=o.set, dt=o.dt, window=w, output=k, ratio=ratios[k].mean(), se=ratios[k].sem(), reps=len(ratios)))
    return pd.DataFrame(rows)


def pooled(bs):
    """ The mean over sets of each set's mean ratio, with the SE between sets, and the largest |z| of a set against 1 """
    g = bs.groupby(['dt', 'window', 'output'])
    out = g.agg(mean=('ratio', 'mean'), se_between=('ratio', 'sem'), lo=('ratio', 'min'), hi=('ratio', 'max'), sets=('ratio', 'size')).reset_index()
    z = bs.assign(z=(bs.ratio - 1)/bs.se).groupby(['dt', 'window', 'output']).z.agg(lambda x: x.abs().max()).reset_index(name='max_abs_z')
    return out.merge(z, on=['dt', 'window', 'output'])


def target_table(runs, ode):
    """ Each set's target values: the ODE, agents (mean ± SE over reps), agents ÷ ODE, and whether each is inside the target's interval """
    rows = []
    for o in ode:
        rs = pd.DataFrame([r.targets for r in runs if r.set == o.set and r.dt == o.dt])
        for k in run_targets:
            mean, lo, hi = run24.target_interval(cb.targets[k])
            a, ov = rs[k].mean(), o.targets[k]
            rows.append(dict(set=o.set, dt=o.dt, target=k, ode=ov, agents=a, agents_se=rs[k].sem(), ratio=a/ov, ratio_se=rs[k].sem()/ov,
                             target_mean=mean, target_lo=lo, target_hi=hi, ode_in=lo <= ov <= hi, agents_in=lo <= a <= hi))
    return pd.DataFrame(rows)


def pitfall_checks(runs):
    """ The pitfall checks, pooled over sets and reps for each step """
    rows = []
    for dt in dts:
        rs = [r for r in runs if r.dt == dt]
        c = sc.objdict({k: sum(r.checks[k] for r in rs) for k in rs[0].checks})
        steps = sc.ddict(int)
        for r in rs:
            for k, v in r.course_steps.items():
                steps[int(k)] += v
        dty = rs[0].dt_year
        mean_steps = c.course_steps_sum/c.courses_timed
        mode_steps = max(steps, key=steps.get)
        rows.append(dict(dt=dt, starts=c.starts, start_while_on_tx=c.start_while_on_tx, start_not_symptomatic=c.start_not_symptomatic,
                         relapse_not_from_success=c.relapse_not_from_success, tb_death_from_success=c.tb_death_from_success,
                         bg_death_rate_success=c.bg_death_from_success/(c.success_person_steps*dty), mu=cb.mu,
                         course_years_mean=mean_steps*dty, share_courses_at_mode=steps[mode_steps]/c.courses_timed, mode_years=mode_steps*dty,
                         alive_min_share=min(r.alive_range[0] for r in rs)/n_agents, alive_max_share=max(r.alive_range[1] for r in rs)/n_agents))
    return pd.DataFrame(rows)


main_outputs = ['prevalence', 'prop_asymptomatic', 'incidence_per_case', 'notifications_per_case', 'on_treatment_per_case', 'mortality_per_case']
short = dict(prevalence='Prevalence', prop_asymptomatic='Share asymptomatic', incidence_per_case='Incidence\nper case',
             notifications_per_case='Notifications\nper case', on_treatment_per_case='On treatment\nper case', mortality_per_case='Mortality\nper case')
window_colors = dict(early='#1f77b4', late='#ff7f0e', const='#2ca02c')


def plot_agents_vs_ode(bs, sets):
    """ Agents ÷ ODE by set, for each main output and window, at each step; each set labelled by its k_latent and inf_asy """
    fig, axes = plt.subplots(2, len(main_outputs), figsize=(22, 7), sharey='row')
    order = np.argsort([s.theta['k_latent'] for s in sets]) # Sets ordered along the ridge
    for row, dt in enumerate(dts):
        for col, k in enumerate(main_outputs):
            ax = axes[row, col]
            for j, w in enumerate(windows):
                d = bs[(bs.dt == dt) & (bs.output == k) & (bs.window == w)].set_index('set').loc[order]
                x = np.arange(len(order)) + (j - 1)*0.25
                ax.errorbar(x, d.ratio, yerr=d.se, fmt='o', ms=3, capsize=1.5, color=window_colors[w], label=f'{windows[w][0]}–{windows[w][1] - 1}')
            ax.axhline(1, color='k', lw=0.6)
            ax.set_xticks(range(len(order)), [f"{sets[i].theta['k_latent']:.1f}\n{sets[i].theta['inf_asy']:.2f}" for i in order], fontsize=6)
            if row == 0:
                ax.set_title(short[k], fontsize=10)
            if col == 0:
                ax.set_ylabel(f'{dt.capitalize()}ly steps\nagents ÷ ODE (mean ± SE, {reps} reps)')
    axes[0, 0].legend(fontsize=7, frameon=False)
    fig.supxlabel('Sets, ordered by k_latent (top line of each label; inf_asy below)', fontsize=9)
    fig.suptitle(f'Exp 27: agents ÷ ODE by set and window; {n_agents:,} agents', fontsize=11)
    fig.tight_layout()
    fig.savefig(outdir / 'agents_vs_ode.png', dpi=130)
    plt.close(fig)
    return


def plot_over_time(runs, ode, sets, dt='month'):
    """ Prevalence, incidence, notifications and mortality over time at three sets along the ridge (lowest, middle and highest k_latent) """
    order = np.argsort([s.theta['k_latent'] for s in sets])
    pick = [order[0], order[len(order)//2], order[-1]]
    keys = sc.objdict(prevalence='Prevalence (per 100,000)', incidence='Incidence (per 100,000 per year)',
                      notifications='Treatment starts (per 100,000 per year)', mortality='TB mortality (per 100,000 per year)')
    fig, axes = plt.subplots(len(pick), len(keys), figsize=(18, 10), sharex=True)
    for r_i, i in enumerate(pick):
        o = next(x for x in ode if x.set == i and x.dt == dt).yearly
        for c_i, (k, label) in enumerate(keys.items()):
            ax = axes[r_i, c_i]
            d = pd.concat([x.yearly[k] for x in runs if x.set == i and x.dt == dt], axis=1)
            ax.plot(d.index, d.values, color='#7fa7d6', lw=0.5, alpha=0.5)
            ax.plot(d.index, d.mean(axis=1), color='#1f4e8c', lw=1.8, label=f'Agents (mean of {d.shape[1]})')
            ax.plot(o.index, o[k], color='k', lw=1.2, ls='--', label='ODE')
            for y in [1995, 2024]:
                ax.axvline(y, color='k', lw=0.5, ls=':')
            ax.set_ylim(bottom=0)
            if r_i == 0:
                ax.set_title(label, fontsize=10)
            if c_i == 0:
                ax.set_ylabel(f"Set {i}\nk_latent {sets[i].theta['k_latent']:.1f}, inf_asy {sets[i].theta['inf_asy']:.2f}", fontsize=9)
    axes[0, 0].legend(fontsize=8, frameon=False)
    fig.suptitle(f'Exp 27: three sets along the ridge, {dt}ly steps, {n_agents:,} agents', fontsize=11)
    fig.tight_layout()
    fig.savefig(outdir / f'over_time_{dt}.png', dpi=130)
    plt.close(fig)
    return


def plot_targets(tt):
    """ Each set's target values, agents (mean ± SE) and ODE, against the target's interval; monthly steps """
    fig, axes = plt.subplots(2, 4, figsize=(18, 7))
    d = tt[tt.dt == 'month']
    for ax, k in zip(axes.flat, run_targets):
        t = d[d.target == k].sort_values('set')
        ax.axhspan(t.target_lo.iloc[0], t.target_hi.iloc[0], color='grey', alpha=0.2, lw=0, label='Target 95% interval')
        ax.axhline(t.target_mean.iloc[0], color='grey', lw=0.8)
        ax.plot(t.set - 0.12, t.ode, 'k_', ms=10, mew=2, label='ODE')
        ax.errorbar(t.set + 0.12, t.agents, yerr=t.agents_se, fmt='o', ms=3, color='#1f4e8c', capsize=1.5, label='Agents')
        ax.set_title(run24.labels.get(k, k), fontsize=10)
        ax.set_xticks(range(len(t)))
        ax.set_xlabel('Set', fontsize=8)
    axes.flat[0].legend(fontsize=7, frameon=False)
    fig.suptitle(f'Exp 27: target values by set, monthly steps ({reps} reps per set)', fontsize=11)
    fig.tight_layout()
    fig.savefig(outdir / 'targets.png', dpi=130)
    plt.close(fig)
    return


def summary():
    sets, runs, ode = sc.load(fn.sets), sc.load(fn.runs), sc.load(fn.ode)
    bs = by_set(runs, ode)
    pl = pooled(bs)
    tt = target_table(runs, ode)
    pc = pitfall_checks(runs)
    bs.to_csv(outdir / 'by_set.csv', index=False)
    pl.to_csv(outdir / 'pooled.csv', index=False)
    tt.to_csv(outdir / 'targets.csv', index=False)
    pc.to_csv(outdir / 'pitfall_checks.csv', index=False)
    plot_agents_vs_ode(bs, sets)
    for dt in dts:
        plot_over_time(runs, ode, sets, dt)
    plot_targets(tt)
    tp = tt.groupby(['dt', 'target']).agg(ratio=('ratio', 'mean'), se_between=('ratio', 'sem'), lo=('ratio', 'min'), hi=('ratio', 'max'),
                                         ode_in=('ode_in', 'sum'), agents_in=('agents_in', 'sum')).reset_index()
    tp.to_csv(outdir / 'targets_pooled.csv', index=False)
    with pd.option_context('display.width', 250, 'display.max_columns', 30, 'display.max_rows', 500, 'display.precision', 4):
        keep = main_outputs + ['mortality_tx_per_case', 'recent', 'remote', 'cleared']
        print('Agents ÷ ODE, mean over the 12 sets ± SE between sets (min–max of set means; largest |z| of a set):')
        print(pl[pl.output.isin(keep)].to_string(index=False))
        print('\nTargets, agents ÷ ODE over sets, and the number of sets inside the target interval:'); print(tp.to_string(index=False))
        print('\nPitfall checks (pooled over sets and reps):'); print(pc.to_string(index=False))
    return


if __name__ == '__main__':
    args = sys.argv[1:]
    if 'sets' in args:
        draw_sets()
    elif 'check' in args:
        mapping_check()
    else:
        if 'summary' not in args:
            run_all()
        summary()
