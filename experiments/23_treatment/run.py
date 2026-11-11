"""
Experiment 23, v4 environment: treatment (tbode4.TBTreat) in v4's ODE and agents, with a fixed 6-month course against a rate. See README.md.
Run from the repo root:

    PYTHONPATH=v4 v4/.venv/bin/python experiments/23_treatment/run.py          # every run, then the summary
    PYTHONPATH=v4 v4/.venv/bin/python experiments/23_treatment/run.py summary  # summaries and figures from the saved runs
"""
import os
os.environ.setdefault('NUMBA_NUM_THREADS', '1') # One thread per process (Starsim's kernels otherwise use every core); before importing starsim
import sys
import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss
import matplotlib.pyplot as plt

outdir = sc.path(os.environ.get('TB_EXP23_FOLDER', sc.thispath(__file__))) # TB_EXP23_FOLDER redirects the outputs (exp 26 reruns this at newer pins)
root = sc.thispath(__file__).parent.parent
sys.path.insert(0, str(root / 'v4'))
import tbode4 # noqa: E402
from tbode4 import fastode as fo # noqa: E402

# Exp 21's set 1 natural history (the starter script's parameters)
pi, phi, inf_non, inf_asy = 0.3518, 0.01194, 0.3153, 0.3735
pars = dict(beta=22.10, trans_asymp=0.7488, rr_reinfection_cleared=pi, rr_reinfection_rec=pi, rr_reinfection_treat=pi, rr_reinfection_inf=pi,
            inf_cle=0.6638, inf_non=inf_non, inf_asy=inf_asy, k_latent=7.517, inf_non_remote=phi*inf_non, inf_asy_remote=phi*inf_asy,
            non_rec=0.2207, non_asy=0.2776, asy_non=1.612, asy_sym=0.7560, sym_asy=0.7593, sym_dead=0.3959, mu=1/70)
lai = tbode4.lai_tpt_treatment() # tb_LAI_TPT's best-fit coverage and efficacy over time
const = dict(coverage=lai.coverage(2024), success=lai.efficacy(2024), p_relapse=lai.p_relapse) # The equilibrium scenario: 2024's values

scenarios = sc.objdict(
    eq   = sc.objdict(start=2000, stop=2050, window=(2005, 2050), label='Constant coverage, at equilibrium'),
    ramp = sc.objdict(start=1980, stop=2036, window=(1995, 2036), label="tb_LAI_TPT's scale-up"),
)
courses = ['rate', 'fixed']
dts = sc.objdict(month=ss.months(1), week=ss.weeks(1))
ode_stages = [2, 4, 12] # Extra ODE-only arms: a course of k Erlang stages (k = 1 is the rate)
reps, n_agents = 40, 200_000 # Treated prevalence is about 80/100k: 200,000 agents give about 160 people with TB
fixed_stages_init = 50 # fastode's stand-in for a fixed course, for the equilibrium scenario's starting state
outputs = ['prevalence', 'incidence', 'relapse', 'notifications', 'on_treatment', 'mortality', 'mortality_tx']
fn = sc.objdict(runs=outdir / 'runs.obj', ode=outdir / 'ode.obj')


def seed(scen, course, dt, rep):
    """ New seeds for this experiment: 23000 + 1000·scenario + 100·arm + rep """
    return 23000 + 1000*list(scenarios).index(scen) + 100*(2*courses.index(course) + list(dts).index(dt)) + rep


def tx_kwargs(scen):
    if scen == 'eq':
        return dict(tx_rate=float(tbode4.theta_from_coverage(const['coverage'])), tx_success=const['success'], tx_p_relapse=const['p_relapse'])
    t = tbode4.lai_tpt_treatment() # New timeseries objects for each sim
    return dict(tx_rate=t.tx_rate, tx_success=t.tx_success, tx_p_relapse=t.p_relapse)


def starting_state(scen, course):
    """
    The starting shares (TBTreat's init=), the starting NON0 share (for incidence), and the share of non-infectious people already
    asymptomatic in this episode. Scale-up: the untreated equilibrium. Equilibrium scenario: fastode's treated equilibrium with the
    course as a rate, or as 50 Erlang stages for the fixed course (everyone on treatment then starts at the beginning of a course)
    """
    if scen == 'ramp':
        p = fo.make_pars(**pars)
    else:
        k = 1 if course == 'rate' else fixed_stages_init
        p = fo.make_pars(**pars, tx_rate=float(tbode4.theta_from_coverage(const['coverage'])), tx_success=const['success'],
                         tx_p_relapse=const['p_relapse'], tx_stages=k)
    x = fo.equilibrium(p)
    d = dict(zip(fo.names(p), x))
    init = {k: d[k] for k in fo.state_names if k not in ['non0', 'non1']}
    init['non_infectious'] = d['non0'] + d['non1']
    init['treated'], init['relapsing'] = d['treated'], d['relapsing']
    k = p['tx_stages']
    init['treatment_success'] = sum(d[f'treatment_success_{i}'] for i in range(k))
    init['treatment_fail'] = sum(d[f'treatment_fail_{i}'] for i in range(k))
    init = {k: v for k, v in init.items() if v > 0}
    return init, d['non0'], d['non1']/(d['non0'] + d['non1'])


stocks = ['prevalence', 'on_treatment']
flows_at_start = True # Where Starsim puts a step's flows: index i is the step from t_i to t_(i+1) from starsim d83f4f87 on (exp 26); the step ending at t_i before (exp 23)


def yearly(out, dt_year):
    """ Yearly means per 100,000 of each output (per step: flows per year, stocks at the step's end), dated by the middle of each step """
    year = np.floor(out.year[1:] - dt_year/2 + 1e-9).astype(int)
    df = pd.DataFrame({k: (out[k][1:] if (k in stocks or not flows_at_start) else out[k][:-1])*1e5 for k in outputs})
    df['year'] = year
    return df.groupby('year').mean()


def make_sim(scen, course, dt, solver, stages=None, rand_seed=0):
    s = scenarios[scen]
    init, non0, share_flagged = starting_state(scen, course)
    tb = tbode4.TBTreat(init=init, course=course, stages=stages, **tx_kwargs(scen), **pars)
    kw = dict(n_agents=n_agents, total_pop=n_agents, rand_seed=rand_seed,
              analyzers=[tbode4.Incidence(share_flagged=share_flagged), tbode4.TxChecks()]) if solver == 'abm' else dict(total_pop=1)
    sim = ss.Sim(tb, solver=solver, start=s.start, stop=s.stop, dt=dts[dt], verbose=0, **kw)
    return sim, non0


def run_agents(scen, course, dt, rep):
    sim, _ = make_sim(scen, course, dt, 'abm', rand_seed=seed(scen, course, dt, rep))
    sim.run()
    inc, chk = sim.analyzers[0], sim.analyzers[1]
    out = tbode4.tx_outputs(sim, inc)
    out.relapse = inc.relapses/np.asarray(sim.results.alive.values, dtype=float)/sim.t.dt_year # Relapses counted by the analyzer
    alive = np.asarray(sim.results.alive.values, dtype=float)
    return sc.objdict(scen=scen, course=course, dt=dt, rep=rep, seed=seed(scen, course, dt, rep), yearly=yearly(out, sim.t.dt_year),
                      checks=dict(chk.c), course_steps=dict(chk.course_steps), by_year={k: v.tolist() for k, v in chk.by_year.items()},
                      alive_range=(float(alive.min()), float(alive.max())), dt_year=sim.t.dt_year)


def run_ode(scen, course, dt, stages=None):
    sim, non0 = make_sim(scen, course, dt, 'ode', stages=stages)
    sim.run()
    out = tbode4.tx_outputs(sim, tbode4.ode_incidence(sim, non0=non0))
    name = course if not stages else f'stages_{stages}'
    r = sim.results.tb
    alive = np.asarray(sim.results.alive.values, dtype=float)
    deaths_tx, fails = r['tb_death_tx'].values.sum(), r['tx_fail'].values.sum()
    return sc.objdict(scen=scen, arm=name, dt=dt, yearly=yearly(out, sim.t.dt_year), dt_year=sim.t.dt_year,
                      share_fail_courses_dying=float(deaths_tx/(deaths_tx + fails)), alive_range=(float(alive.min()), float(alive.max())))


def ode_equilibria():
    """ Fixed course against a rate, and k Erlang stages, at equilibrium: v4's ODE run for 400 years at constant 2024 treatment """
    rows = []
    for arm, course, stages in [('rate', 'rate', None), ('fixed', 'fixed', None)] + [(f'stages_{k}', 'fixed', k) for k in ode_stages]:
        for dt in dts:
            init, non0, _ = starting_state('eq', 'rate')
            tb = tbode4.TBTreat(init=init, course=course, stages=stages, **tx_kwargs('eq'), **pars)
            sim = ss.Sim(tb, solver='ode', start=0, stop=400, dt=dts[dt], total_pop=1, verbose=0).run()
            out = tbode4.tx_outputs(sim, tbode4.ode_incidence(sim, non0=non0))
            n = int(round(1/sim.t.dt_year)) # The last year
            rows.append(dict(arm=arm, dt=dt, **{k: float(np.mean(out[k][-n:]))*1e5 for k in outputs}))
    return pd.DataFrame(rows)


def run_all():
    T = sc.timer()
    ode = sc.objdict(runs=[], equilibria=None)
    for scen in scenarios:
        for dt in dts:
            for course in courses:
                ode.runs.append(run_ode(scen, course, dt))
            for k in ode_stages:
                ode.runs.append(run_ode(scen, 'fixed', dt, stages=k))
    ode.equilibria = ode_equilibria()
    sc.save(fn.ode, ode)
    print(f'ODE runs done ({T.toc(output=True):.0f} s)')
    groups = [(s, d, c) for s in scenarios for d in dts for c in courses]
    partdir = outdir / 'runs_by_group'
    partdir.mkdir(exist_ok=True)
    runs = []
    for i, (s, d, c) in enumerate(groups): # Saved group by group, so a stopped run resumes where it left off
        part = partdir / f'{s}_{d}_{c}.obj'
        if part.exists():
            group = sc.load(part)
        else:
            jobs = [dict(scen=s, course=c, dt=d, rep=r) for r in range(reps)]
            group = sc.parallelize(run_agents, iterkwargs=jobs, ncpus=tbode4.ncpus)
            sc.save(part, group)
        runs += group
        print(f'Group {i + 1} of {len(groups)} ({s}, {d}, {c}): done ({T.toc(output=True):.0f} s)', flush=True)
    sc.save(fn.runs, runs)
    tbode4.save_versions(outdir, 'v4')
    print(f'{len(runs)} agent runs done ({T.toc(output=True):.0f} s)')
    return


# Summaries and figures

def window_mean(df, scen):
    lo, hi = scenarios[scen].window
    return df[(df.index >= lo) & (df.index < hi)].mean()


def agents_vs_ode(runs, ode):
    """ Agents ÷ ODE for each output, averaged over the scenario's window: mean ± SE across reps, and the 10–90% range """
    rows = []
    for scen in scenarios:
        for dt in dts:
            for course in courses:
                o = window_mean(next(r for r in ode.runs if r.scen == scen and r.dt == dt and r.arm == course).yearly, scen)
                ratios = pd.DataFrame([window_mean(r.yearly, scen)/o for r in runs if r.scen == scen and r.dt == dt and r.course == course])
                row = dict(scenario=scen, dt=dt, course=course, reps=len(ratios))
                for k in outputs:
                    row[k] = ratios[k].mean(); row[k + '_se'] = ratios[k].std()/np.sqrt(len(ratios))
                    row[k + '_p10'], row[k + '_p90'] = ratios[k].quantile(0.1), ratios[k].quantile(0.9)
                rows.append(row)
    return pd.DataFrame(rows)


def pitfall_checks(runs, ode):
    """ The pitfall checks, pooled over reps for each scenario, course and dt """
    rows = []
    for scen in scenarios:
        for dt in dts:
            for course in courses:
                rs = [r for r in runs if r.scen == scen and r.dt == dt and r.course == course]
                c = sc.objdict({k: sum(r.checks[k] for r in rs) for k in rs[0].checks})
                steps = sc.ddict(int)
                for r in rs:
                    for k, v in r.course_steps.items():
                        steps[int(k)] += v
                dty = rs[0].dt_year
                mean_steps = c.course_steps_sum/c.courses_timed
                sd_steps = np.sqrt(c.course_steps_sq/c.courses_timed - mean_steps**2)
                mode_steps = max(steps, key=steps.get)
                theta = float(tbode4.theta_from_coverage(const['coverage']))
                row = dict(scenario=scen, dt=dt, course=course, starts=c.starts,
                           start_while_on_tx=c.start_while_on_tx, start_not_symptomatic=c.start_not_symptomatic,
                           relapse_not_from_success=c.relapse_not_from_success, tb_death_from_success=c.tb_death_from_success,
                           restart_same_step=c.restart_same_step/max(1, c.fails),
                           restart_expected=(1 - np.exp(-theta*dty/2)) if (scen == 'eq' and course == 'fixed') else np.nan, # With a rate, starting and failing in one step also counts
                           bg_death_rate_success=c.bg_death_from_success/(c.success_person_steps*dty), mu=pars['mu'],
                           course_years_mean=mean_steps*dty, course_years_sd=sd_steps*dty,
                           share_courses_at_mode=steps[mode_steps]/c.courses_timed, mode_years=mode_steps*dty,
                           alive_min=min(r.alive_range[0] for r in rs), alive_max=max(r.alive_range[1] for r in rs))
                rows.append(row)
    return pd.DataFrame(rows)


def success_by_year(runs):
    """ Scale-up: the share of courses set to cure, by year of starting (agents, pooled), against efficacy at the middle of the year """
    rows = []
    for dt in dts:
        for course in courses:
            tot = sc.ddict(lambda: np.zeros(2))
            for r in runs:
                if r.scen == 'ramp' and r.dt == dt and r.course == course:
                    for y, v in r.by_year.items():
                        tot[int(y)] += v
            for y in sorted(tot):
                if 1995 <= y < 2036 and tot[y][0]:
                    rows.append(dict(dt=dt, course=course, year=y, courses=int(tot[y][0]), share_cure=tot[y][1]/tot[y][0],
                                     efficacy_mid_year=float(lai.efficacy(y + 0.5))))
    return pd.DataFrame(rows)


def fixed_vs_rate(ode):
    """ ODE: each arm ÷ the rate, by year (scale-up) """
    rows = []
    for dt in dts:
        base = next(r for r in ode.runs if r.scen == 'ramp' and r.dt == dt and r.arm == 'rate').yearly
        for r in ode.runs:
            if r.scen == 'ramp' and r.dt == dt and r.arm != 'rate':
                ratio = r.yearly/base
                ratio['arm'], ratio['dt'] = r.arm, dt
                rows.append(ratio.reset_index())
    return pd.concat(rows, ignore_index=True)


colors = sc.objdict(rate='#1f77b4', fixed='#d62728')
panel_labels = sc.objdict(prevalence='Prevalence\n(per 100,000)', incidence='Incidence\n(per 100,000 per year)',
                          notifications='Notifications\n(per 100,000 per year)', on_treatment='On treatment\n(per 100,000)',
                          mortality='TB mortality\n(per 100,000 per year)')


def plot_over_time(runs, ode, scen, dt):
    fig, axes = plt.subplots(1, 5, figsize=(20, 4))
    for ax, (key, label) in zip(axes, panel_labels.items()):
        for course in courses:
            d = pd.concat([r.yearly[key] for r in runs if r.scen == scen and r.dt == dt and r.course == course], axis=1)
            ax.fill_between(d.index, d.quantile(0.1, axis=1), d.quantile(0.9, axis=1), color=colors[course], alpha=0.15, lw=0)
            ax.plot(d.index, d.mean(axis=1), color=colors[course], lw=1.2, alpha=0.8, label=f'Agents, {course} (mean, 10–90%)')
            o = next(r for r in ode.runs if r.scen == scen and r.dt == dt and r.arm == course).yearly
            ax.plot(o.index, o[key], color=colors[course], lw=1.5, ls='--', label=f'ODE, {course}')
            part = dict(incidence='relapse', mortality='mortality_tx').get(key)
            if part and course == 'fixed':
                ax.fill_between(o.index, 0, o[part], color='grey', alpha=0.3, lw=0, label='From relapse (ODE, fixed)' if part == 'relapse' else 'On treatment (ODE, fixed)')
        if scen == 'ramp':
            for y in [1995, 2024]:
                ax.axvline(y, color='k', lw=0.6, ls=':')
        ax.set_title(label, fontsize=10)
        ax.set_xlabel('Year')
        ax.set_ylim(bottom=0)
    axes[0].legend(fontsize=7, frameon=False, loc='lower left')
    axes[1].legend(fontsize=7, frameon=False, loc='lower left', handles=[h for h in axes[1].get_legend_handles_labels()[0] if 'relapse' in str(h.get_label())])
    axes[4].legend(fontsize=7, frameon=False, loc='lower left', handles=[h for h in axes[4].get_legend_handles_labels()[0] if 'On treatment' in str(h.get_label())])
    fig.suptitle(f'Exp 23: {scenarios[scen].label}, {dt}ly steps; {n_agents:,} agents, {reps} reps per course', fontsize=10)
    fig.tight_layout()
    fig.savefig(outdir / f'over_time_{scen}_{dt}.png', dpi=130)
    plt.close(fig)
    return


def plot_fixed_vs_rate(fvr):
    fig, axes = plt.subplots(1, 5, figsize=(20, 3.6))
    arm_colors = dict(fixed='#d62728', stages_12='#ff7f0e', stages_4='#2ca02c', stages_2='#9467bd')
    for ax, (key, label) in zip(axes, panel_labels.items()):
        for arm, color in arm_colors.items():
            d = fvr[(fvr.arm == arm) & (fvr.dt == 'month')]
            ax.plot(d.year, d[key], color=color, lw=1.3, label=arm.replace('_', ' = k ').replace('stages', 'Erlang'))
        ax.axhline(1, color='k', lw=0.6)
        for y in [1995, 2024]:
            ax.axvline(y, color='k', lw=0.6, ls=':')
        ax.set_title(label.split('\n')[0] + ' ÷ rate', fontsize=10)
        ax.set_xlabel('Year')
        ax.set_xlim(1985, 2036)
    axes[0].legend(fontsize=8, frameon=False)
    fig.suptitle("Exp 23, ODE: each course ÷ a course at a rate (mean 6 months), tb_LAI_TPT's scale-up, monthly steps", fontsize=10)
    fig.tight_layout()
    fig.savefig(outdir / 'fixed_vs_rate.png', dpi=130)
    plt.close(fig)
    return


def plot_agents_vs_ode(avo):
    keys = ['prevalence', 'incidence', 'notifications', 'on_treatment', 'mortality']
    fig, axes = plt.subplots(1, 2, figsize=(13, 4), sharey=True)
    for ax, scen in zip(axes, scenarios):
        d = avo[avo.scenario == scen].reset_index(drop=True)
        for j, row in d.iterrows():
            x = np.arange(len(keys)) + (j - (len(d) - 1)/2)*0.18
            y = [row[k] for k in keys]; e = [row[k + '_se'] for k in keys]
            ax.errorbar(x, y, yerr=e, fmt='o' if row['dt'] == 'month' else 's', color=colors[row['course']], ms=4, capsize=2, label=f"{row['course']}, {row['dt']}ly")
        ax.axhline(1, color='k', lw=0.6)
        ax.set_xticks(range(len(keys)), [panel_labels[k].split('\n')[0] for k in keys], fontsize=8)
        ax.set_title(f'{scenarios[scen].label} ({scenarios[scen].window[0]}–{scenarios[scen].window[1] - 1})', fontsize=10)
    axes[0].set_ylabel(f'Agents ÷ ODE (mean ± SE, {reps} reps)')
    axes[0].legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(outdir / 'agents_vs_ode.png', dpi=130)
    plt.close(fig)
    return


def summary():
    runs, ode = sc.load(fn.runs), sc.load(fn.ode)
    avo = agents_vs_ode(runs, ode)
    checks = pitfall_checks(runs, ode)
    sby = success_by_year(runs)
    fvr = fixed_vs_rate(ode)
    eqs = ode.equilibria
    avo.to_csv(outdir / 'agents_vs_ode.csv', index=False)
    checks.to_csv(outdir / 'pitfall_checks.csv', index=False)
    sby.to_csv(outdir / 'success_by_year.csv', index=False)
    fvr.to_csv(outdir / 'fixed_vs_rate.csv', index=False)
    eqs.to_csv(outdir / 'ode_equilibria.csv', index=False)
    shares = pd.DataFrame([dict(scen=r.scen, arm=r.arm, dt=r.dt, share_fail_courses_dying=r.share_fail_courses_dying) for r in ode.runs])
    shares.to_csv(outdir / 'ode_fail_death_share.csv', index=False)
    for scen in scenarios:
        for dt in dts:
            plot_over_time(runs, ode, scen, dt)
    plot_fixed_vs_rate(fvr)
    plot_agents_vs_ode(avo)
    with pd.option_context('display.width', 250, 'display.max_columns', 50, 'display.precision', 4):
        cols = ['scenario', 'dt', 'course'] + [f'{k}' for k in outputs]
        print('Agents ÷ ODE (mean over reps):'); print(avo[cols].to_string(index=False))
        print('\nSE:'); print(avo[['scenario', 'dt', 'course'] + [k + '_se' for k in outputs]].to_string(index=False))
        print('\nODE equilibria (per 100,000):'); print(eqs.to_string(index=False))
        print('\nPitfall checks:'); print(checks.to_string(index=False))
        print('\nShare of failing courses ending in TB death (ODE):'); print(shares.to_string(index=False))
        print('\nScale-up, share of courses set to cure by start year (agents, monthly fixed):')
        print(sby[(sby.dt == 'month') & (sby.course == 'fixed')].to_string(index=False))
        d = fvr[fvr.dt == 'month']
        print('\nODE fixed ÷ rate, scale-up, monthly, selected years:')
        print(d[d.year.isin([1996, 1998, 2000, 2005, 2010, 2015, 2020, 2024, 2030, 2035])][['arm', 'year'] + [k for k in outputs]].to_string(index=False))
    return


if __name__ == '__main__':
    if 'summary' not in sys.argv[1:]:
        run_all()
    summary()
