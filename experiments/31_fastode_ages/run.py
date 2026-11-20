"""
Experiment 31, v4 environment: a fast standalone ODE with ages and demography (tbode4.fastode_age), against v4's ODE (exp 30), its
exact-aging limit, and itself at smaller steps. See README.md. Run from the repo root:

    PYTHONPATH=v4 v4/.venv/bin/python experiments/31_fastode_ages/run.py copy     # arm C: option 1 against exp 30's quarter-year Starsim ODE
    PYTHONPATH=v4 v4/.venv/bin/python experiments/31_fastode_ages/run.py steps    # arms X and X0: option 2 at each step, both force-of-infection timings
    PYTHONPATH=v4 v4/.venv/bin/python experiments/31_fastode_ages/run.py start    # arm S: option 2's own starting tables (200, 100, 50 years)
    PYTHONPATH=v4 v4/.venv/bin/python experiments/31_fastode_ages/run.py prior    # arm P: the whole draw (start, then 1950–2041) at 20 prior draws
    PYTHONPATH=v4 v4/.venv/bin/python experiments/31_fastode_ages/run.py summary  # tables and figures from the saved runs (and the agents comparison)

Environment variable: EXP31_NCPUS (cores; default tbode4.ncpus).
"""
import os
os.environ.setdefault('NUMBA_NUM_THREADS', '1') # One thread per process; before importing starsim
import sys
import importlib.util
import numpy as np
import pandas as pd
import sciris as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

here = sc.thispath(__file__)
root = here.parent.parent
sys.path.insert(0, str(root / 'v4'))
import tbode4 # noqa: E402
from tbode4 import fastode as fo, fastode_age as fa # noqa: E402


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


r28 = load_module('run28', root / 'experiments' / '28_age_demography' / 'run.py') # Exp 28's windows and by_window()
cb = load_module('calib24', root / 'experiments' / '24_treatment_calibration' / 'calib.py') # Exp 24's priors and parameterisation

exp30 = root / 'experiments' / '30_tb_corrected_demography'
exp27 = root / 'experiments' / '27_calibrated_treatment'
ncpus = int(os.environ.get('EXP31_NCPUS', tbode4.ncpus))
start, stop = 1950, 2041
sets_used = [0, 6]
starsim_dt = sc.objdict(month=1/12, week=7/365.25) # Starsim's steps in exp 30 (years)
steps = sc.objdict(quarter=0.25, month=1/12, half=1/24, week=1/52) # Option 2's steps (and age widths); 1/52 keeps WPP's age edges and the 6-month course on the grid
fois = ['mid', 'start']
burns = [200, 100, 50]
n_prior, prior_seed = 20, 31
windows = r28.windows # 1950–1994, 1995–2024, 2025–2040
groups5 = list(fa.groups5)
burden = ['alive', 'prevalence', 'prop_asymptomatic', 'incidence', 'mortality', 'notifications', 'on_treatment', 'recent', 'remote', 'cleared'] + \
         [f'prev_{g}' for g in groups5] + [f'share_{g}' for g in groups5] + ['prev_15+']
target_keys = ['prevalence', 'prevalence_15plus', 'prop_asymptomatic', 'prev_to_inc', 'mort_to_inc', 'prev_to_notif', 'inc_ratio', 'mort_ratio', 'prev_2040_ratio']
fn = {k: here / f'{k}.obj' for k in ['copy', 'steps', 'start', 'prior']} # A plain dict: on an objdict, fn.copy would be its copy() method
sets = sc.load(exp27 / 'sets.obj')
starts30 = sc.load(exp30 / 'start.obj')


#%% The model's inputs

def pars(nat, trend):
    """ fastode_age's parameters for one set: its natural history with mu = 0, and exp 24's schedules """
    nat = {k: v for k, v in nat.items() if k not in ['mu', 'sym_dead']}
    return fo.make_pars(**nat, mu=0, **fa.schedules(trend))


def untreated_equilibrium(nat):
    """ Exp 27's start without ages (fastode's untreated equilibrium, births replacing deaths at 1/70): shares, and the share of NON already asymptomatic """
    p = fo.make_pars(**nat)
    x = fo.equilibrium(p)
    if x is None:
        return None, None
    d = dict(zip(fo.names(p), x))
    init = {k: d[k] for k in fo.state_names if k not in ['non0', 'non1']}
    init['non_infectious'] = d['non0'] + d['non1']
    return init, d['non1']/(d['non0'] + d['non1'])


def result(out, h, secs, **kw):
    """ What is kept of a run: yearly means (as exp 30 kept Starsim's and agents' runs), the calendar summary (to compare steps), run time """
    s = sc.objdict(kw)
    s.yearly = fa.yearly(out, h)
    s.cal = fa.calendar(out)
    s.secs = secs
    return s


#%% Arm C: option 1 against Starsim's quarter-year ODE

def run_copy_one(i, dt):
    s = sets[i]
    st = starts30[f'set{i}']
    T = sc.timer()
    out = fa.run_rate(pars(s.nat, s.trend), 0.25, start, stop, starsim_dt[dt], table=st.table, share_flagged=st.share_flagged, rtol=1e-6, atol=1e-12)
    return result(out, starsim_dt[dt], T.toc(output=True), set=i, dt=dt, nfev=out.nfev)


def run_copy():
    jobs = [dict(i=i, dt=dt) for i in sets_used for dt in starsim_dt]
    res = sc.parallelize(run_copy_one, iterkwargs=jobs, ncpus=ncpus)
    sc.save(fn['copy'], res)
    tbode4.save_versions(here, 'v4')
    print('copy: ' + ', '.join(f'set {r.set} {r.dt} {r.secs:.1f} s' for r in res), flush=True)
    return res


#%% Arms X and X0: option 2 at each step

def run_step_one(i, step, foi):
    s = sets[i]
    st = starts30[f'set{i}']
    T = sc.timer()
    out = fa.run_cohort(pars(s.nat, s.trend), steps[step], start, stop, table=st.table, share_flagged=st.share_flagged, foi=foi)
    r = result(out, steps[step], T.toc(output=True), set=i, step=step, foi=foi)
    r.targets = fa.targets(out, s.trend.over_detection)
    return r


def run_steps():
    jobs = [dict(i=i, step=k, foi=f) for i in sets_used for k in steps for f in fois]
    res = sc.parallelize(run_step_one, iterkwargs=jobs, ncpus=ncpus)
    sc.save(fn['steps'], res)
    tbode4.save_versions(here, 'v4')
    print('steps: ' + ', '.join(f'set {r.set} {r.step} {r.foi} {r.secs:.1f} s' for r in res), flush=True)
    return res


#%% Arm S: option 2's own starting tables

def run_start_one(i, years, step='month'):
    s = sets[i]
    p = pars(s.nat, s.trend)
    init, sf = untreated_equilibrium(s.nat)
    T = sc.timer()
    st = fa.settle(p, steps[step], init, sf, years=years)
    table, flagged = fa.start_table(st)
    secs_start = T.toc(output=True)
    out = fa.run_cohort(p, steps[step], start, stop, table=table, share_flagged=flagged)
    r = result(out, steps[step], T.toc(output=True), set=i, years=years, step=step, drift=st.drift, flagged=flagged, secs_start=secs_start)
    r.targets = fa.targets(out, s.trend.over_detection)
    r.table = table
    return r


def run_start():
    jobs = [dict(i=i, years=y) for i in sets_used for y in burns]
    res = sc.parallelize(run_start_one, iterkwargs=jobs, ncpus=ncpus)
    sc.save(fn['start'], res)
    tbode4.save_versions(here, 'v4')
    print('start: ' + ', '.join(f'set {r.set} {r.years} y: start {r.secs_start:.1f} s, drift {r.drift:.1e}' for r in res), flush=True)
    return res


#%% Arm P: the whole draw at prior draws

def prior_draws():
    """ n_prior draws from exp 24's priors with an untreated endemic equilibrium (exp 27's start without ages) """
    rng = np.random.default_rng(prior_seed)
    out, tried = [], 0
    while len(out) < n_prior:
        theta = cb.from_unit(rng.random((1, cb.d)))[0]
        tried += 1
        nat, trend = cb.to_pars(theta)
        init, sf = untreated_equilibrium(nat)
        if init is not None:
            out.append(sc.objdict(draw=len(out), theta=dict(zip(cb.names, map(float, theta))), nat=nat, trend=trend, init=init, flagged=sf))
    return out, tried


def run_prior_one(draw, step):
    p = pars(draw.nat, draw.trend)
    T = sc.timer()
    st = fa.settle(p, steps[step], draw.init, draw.flagged, years=burns[0])
    table, flagged = fa.start_table(st)
    X = st.X
    endemic = float((X[fa.A] + X[fa.Y]).sum()/X[:st.layout.ns].sum()) > 1e-9
    if not endemic:
        return sc.objdict(draw=draw.draw, step=step, endemic=False, secs=T.toc(output=True))
    out = fa.run_cohort(p, steps[step], start, stop, table=table, share_flagged=flagged)
    r = result(out, steps[step], T.toc(output=True), draw=draw.draw, step=step, endemic=True)
    r.targets = fa.targets(out, draw.trend.over_detection)
    return r


def run_prior():
    draws, tried = prior_draws()
    print(f'prior: {n_prior} draws with an untreated endemic equilibrium, of {tried} tried', flush=True)
    jobs = [dict(draw=d, step=k) for d in draws for k in ['quarter', 'month', 'week']]
    res = sc.parallelize(run_prior_one, iterkwargs=jobs, ncpus=ncpus)
    sc.save(fn['prior'], sc.objdict(draws=draws, tried=tried, runs=res))
    tbode4.save_versions(here, 'v4')
    print('prior: ' + ', '.join(f'{r.draw} {r.step} {r.secs:.0f} s' for r in res), flush=True)
    return res


#%% Summaries

def window_ratio(a, b, keys=burden):
    """ a ÷ b over each window (means of yearly values), and the largest |a ÷ b − 1| in any year from 1951 """
    rows = []
    for w in windows:
        aw, bw = r28.window_mean(a, w), r28.window_mean(b, w)
        rows.append({'window': w, **{k: aw[k]/bw[k] for k in keys if k in aw and k in bw}})
    yrs = a.index[(a.index >= 1951) & (a.index < 2041)].intersection(b.index)
    worst = {}
    for k in keys:
        if k in a and k in b:
            v = a.loc[yrs, k]/b.loc[yrs, k]
            if k in ['notifications', 'on_treatment']:
                v = v[v.index >= 1996]
            worst[k] = float(np.nanmax(np.abs(v - 1)))
    rows.append({'window': 'max |year − 1|', **worst})
    return pd.DataFrame(rows).set_index('window')


def exact_ref(i):
    """ Exp 30's exact-aging reference: 2 × eighth-year − quarter-year, Starsim's ODE at weekly steps """
    ode = sc.load(exp30 / 'ode.obj')
    return 2*ode[f'set{i}_week_eighth'].yearly - ode[f'set{i}_week_quarter'].yearly


def fmt(df, digits=4):
    return df.round(digits).to_string()


def summarize():
    ode = sc.load(exp30 / 'ode.obj')
    txt = []
    add = lambda s: txt.append(s)
    show = ['alive', 'prevalence', 'prop_asymptomatic', 'incidence', 'mortality', 'notifications', 'on_treatment', 'recent', 'remote', 'cleared'] + [f'prev_{g}' for g in groups5]

    # Arm C
    copy = sc.load(fn['copy'])
    rows = []
    add('ARM C. Option 1 (quarter-year groups, constant-rate aging, continuous time, 12 Erlang stages) ÷ exp 30\'s quarter-year Starsim ODE (deterministic)\n')
    for r in copy:
        ref = ode[f'set{r.set}_{r.dt}_quarter'].yearly
        add(f'set {r.set}, against Starsim at {r.dt}ly steps (option 1 output at the same grid; {r.secs:.1f} s, {r.nfev} evaluations)\n' + fmt(window_ratio(r.yearly, ref)[show]) + '\n')
        wr = window_ratio(r.yearly, ref)
        rows.append(dict(set=r.set, dt=r.dt, **{f'{k}_{w}': wr.loc[w, k] for w in wr.index for k in show}))
    for i in sets_used:
        add(f'set {i}: Starsim\'s own monthly ÷ weekly quarter-year runs\n' + fmt(window_ratio(ode[f'set{i}_month_quarter'].yearly, ode[f'set{i}_week_quarter'].yearly)[show]) + '\n')
    pd.DataFrame(rows).to_csv(here / 'copy.csv', index=False)

    # Arms X and X0
    st = sc.load(fn['steps'])
    get = lambda i, k, f: next(r for r in st if r.set == i and r.step == k and r.foi == f)
    rows, trows = [], []
    add('ARMS X, X0. Option 2 (exact aging) at each step ÷ exp 30\'s exact-aging reference (2 × eighth − quarter, Starsim weekly; yearly means), and ÷ option 2 at weekly steps (1/52), mid-step force of infection (calendar summary: stocks on 1 January, flows per calendar year; steps.csv)\n')
    for i in sets_used:
        ex = exact_ref(i)
        wk = get(i, 'week', 'mid')
        for k in steps:
            for f in fois:
                r = get(i, k, f)
                a, b = window_ratio(r.yearly, ex), window_ratio(r.cal, wk.cal) # Against Starsim: yearly means, as exp 30 kept them; between steps: the calendar summary
                add(f'set {i}, step {k}, foi {f} ({r.secs:.2f} s) ÷ exact-aging reference\n' + fmt(a[show]) + '\n')
                rows.append(dict(set=i, step=k, foi=f, against='exact_ref', secs=r.secs, **{f'{c}_{w}': a.loc[w, c] for w in a.index for c in show}))
                rows.append(dict(set=i, step=k, foi=f, against='week_mid', secs=r.secs, **{f'{c}_{w}': b.loc[w, c] for w in b.index for c in show}))
                trows.append(dict(set=i, step=k, foi=f, **{t: r.targets[t]/wk.targets[t] for t in target_keys}))
        # Option 2's own limit: Richardson from half-month and weekly steps, assuming second order (the mid-step force of infection)
        hm, w = get(i, 'half', 'mid').yearly, wk.yearly # Yearly means, to set against the exact-aging reference
        rr = (steps.half/steps.week)**2
        lim = w + (w - hm)/(rr - 1)
        add(f'set {i}: option 2\'s limit (weekly and half-month, extrapolated as second order) ÷ exact-aging reference\n' + fmt(window_ratio(lim, ex)[show]) + '\n')
        add(f'set {i}: option 2 at weekly steps ÷ its extrapolated limit\n' + fmt(window_ratio(w, lim)[show]) + '\n')
    pd.DataFrame(rows).to_csv(here / 'steps.csv', index=False)
    tt = pd.DataFrame(trows)
    tt.to_csv(here / 'steps_targets.csv', index=False)
    add('Calibration outputs (exp 24\'s targets, with prevalence at 15+) ÷ option 2 at weekly steps, mid-step force of infection\n' + tt.round(4).to_string(index=False) + '\n')

    # Arm S
    ss_ = sc.load(fn['start'])
    add('ARM S. Option 2\'s own starting tables (monthly steps; frozen 1950 demography for 200, 100, 50 years)\n')
    srows = []
    for r in ss_:
        t30 = starts30[f'set{r.set}'].table
        fromtab = get(r.set, 'month', 'mid')
        diff = {k: float(np.max(np.abs(r.table[k] - t30[k]))) for k in r.table if k in t30}
        diff_adult = {k: float(np.max(np.abs(r.table[k][:, 1:] - t30[k][:, 1:]))) for k in r.table if k in t30}
        wr = window_ratio(r.yearly, fromtab.yearly)
        add(f'set {r.set}, {r.years} years: drift {r.drift:.1e}; start {r.secs_start:.2f} s; largest |share − exp 30\'s| any age {max(diff.values()):.1e}, ages 1+ {max(diff_adult.values()):.1e}; '
            f'share of NON already asymptomatic {r.flagged:.4f} (exp 30: {starts30[f"set{r.set}"].share_flagged:.4f})\n'
            f'  1950–2041 from this table ÷ from exp 30\'s table (both option 2, monthly)\n' + fmt(wr[show]) + '\n'
            f'  targets from this table ÷ from exp 30\'s table: ' + ', '.join(f'{t} {r.targets[t]/fromtab.targets[t]:.4f}' for t in target_keys) + '\n')
        srows.append(dict(set=r.set, years=r.years, drift=r.drift, secs_start=r.secs_start, max_share_diff=max(diff.values()), max_share_diff_1plus=max(diff_adult.values()),
                          **{f'target_{t}': r.targets[t]/fromtab.targets[t] for t in target_keys}))
    pd.DataFrame(srows).to_csv(here / 'start.csv', index=False)
    for i in sets_used:
        r200 = next(r for r in ss_ if r.set == i and r.years == 200)
        for y in burns[1:]:
            r = next(r for r in ss_ if r.set == i and r.years == y)
            add(f'set {i}: targets with a {y}-year start ÷ a 200-year start: ' + ', '.join(f'{t} {r.targets[t]/r200.targets[t]:.5f}' for t in target_keys))
    add('')

    # Arm P
    pr = sc.load(fn['prior'])
    runs = pr.runs
    add(f'ARM P. The whole draw (start, 200 years frozen, then 1950–2041) at {n_prior} prior draws with an untreated endemic equilibrium (of {pr.tried} tried)\n')
    prow = []
    for d in pr.draws:
        wk = next(r for r in runs if r.draw == d.draw and r.step == 'week')
        for k in ['quarter', 'month']:
            r = next(r for r in runs if r.draw == d.draw and r.step == k)
            if not (r.endemic and wk.endemic):
                prow.append(dict(draw=d.draw, step=k, endemic=f'{r.endemic}/{wk.endemic}'))
                continue
            wr = window_ratio(r.cal, wk.cal) # The calendar summary: yearly means shift with the step when outputs change fast
            prow.append(dict(draw=d.draw, step=k, endemic='yes', secs=r.secs, prev_2020_per100k=wk.targets['prevalence']*1e5,
                             **{f't_{t}': r.targets[t]/wk.targets[t] for t in target_keys},
                             **{f'max_{c}': wr.loc['max |year − 1|', c] for c in ['alive', 'prevalence', 'incidence', 'mortality', 'notifications'] + [f'prev_{g}' for g in groups5]}))
    pdf = pd.DataFrame(prow)
    pdf.to_csv(here / 'prior.csv', index=False)
    with pd.option_context('display.width', 400, 'display.max_columns', 60):
        add('Each draw\'s outputs at quarter-year or monthly steps ÷ at weekly steps (t_: targets; max_: largest |ratio − 1| in the calendar summary, 1951–2040)\n' + pdf.round(4).to_string(index=False) + '\n')
    ok = pdf[pdf.endemic == 'yes']
    for k in ['quarter', 'month']:
        o = ok[ok.step == k]
        tcols = [c for c in o if c.startswith('t_')]
        mcols = [c for c in o if c.startswith('max_')]
        add(f'{k}: largest |target ratio − 1| over draws {np.nanmax(np.abs(o[tcols].values - 1)):.4f}; largest yearly deviation {np.nanmax(o[mcols].values):.4f}; median run time {o.secs.median():.1f} s')
    add('Median run times of the whole draw (s; the laptop slept during a few runs): ' + ', '.join(f'{k} {np.median([r.secs for r in runs if r.step == k]):.1f}' for k in ['quarter', 'month', 'week']) + '\n')

    # Against agents: exp 30's saved runs ÷ option 2
    agents = sc.load(exp30 / 'runs.obj')
    arows = []
    for label, (k, f) in {'month_mid': ('month', 'mid'), 'month_start': ('month', 'start'), 'week_mid': ('week', 'mid'), 'week_start': ('week', 'start')}.items():
        refs = [sc.objdict(set=i, dt=dt, yearly=get(i, k, f).yearly) for i in sets_used for dt in starsim_dt]
        bw = r28.by_window(agents, refs)
        bw['model'] = label
        arows.append(bw)
    refs = [sc.objdict(set=i, dt=dt, yearly=exact_ref(i)) for i in sets_used for dt in starsim_dt]
    bw = r28.by_window(agents, refs); bw['model'] = 'exp30_exact_ref'
    arows.append(bw)
    ag = pd.concat(arows)
    ag.to_csv(here / 'agents.csv', index=False)
    keys = ['prevalence', 'incidence_per_case', 'mortality_per_case', 'notifications_per_case', 'alive', 'recent', 'remote']
    ag['value'] = ag.apply(lambda x: f'{x.ratio:.3f} ± {x.se:.3f}', axis=1)
    with pd.option_context('display.width', 400, 'display.max_columns', 60):
        piv = ag[ag.output.isin(keys)].pivot_table(index=['output', 'model'], columns=['set', 'dt', 'window'], values='value', aggfunc='first')
        add('AGAINST AGENTS. Exp 30\'s agents (set 0: 20 reps, set 6: 60) ÷ each model, mean ± SE over reps of window means\n' + piv.to_string() + '\n')

    text = '\n'.join(txt)
    (here / 'summary.log').write_text(text)
    print(text)
    plot(st, ode)
    return


def plot(st, ode):
    """ Option 2 ÷ the exact-aging reference by year, at each step, for prevalence, 0–14 prevalence, incidence and population """
    items = [('prevalence', 'Prevalence'), ('prev_0-14', 'Prevalence, 0–14'), ('incidence', 'Incidence'), ('alive', 'Population')]
    colors = dict(quarter='tab:red', month='tab:blue', half='tab:green', week='k')
    fig, axs = plt.subplots(len(sets_used), len(items), figsize=(4.4*len(items), 3.4*len(sets_used)), squeeze=False)
    for a, i in enumerate(sets_used):
        ex = exact_ref(i)
        yrs = ex.index[(ex.index >= 1951) & (ex.index < 2041)]
        for b, (k, label) in enumerate(items):
            ax = axs[a, b]
            for step in steps:
                for f, ls in [('mid', '-'), ('start', ':')]:
                    r = next(r for r in st if r.set == i and r.step == step and r.foi == f)
                    ax.plot(yrs, r.yearly.loc[yrs, k]/ex.loc[yrs, k], ls, color=colors[step], lw=1.1, label=f'{step}, {f}' if (a == 0 and b == 0) else None)
            ax.plot(yrs, ode[f'set{i}_week_quarter'].yearly.loc[yrs, k]/ex.loc[yrs, k], color='0.6', lw=1.6, alpha=0.7, label='Starsim, quarter-year groups' if (a == 0 and b == 0) else None)
            ax.axhline(1, color='k', lw=0.5); ax.axhspan(0.998, 1.002, color='0.9', zorder=0)
            ax.set_ylim(0.985, 1.015)
            ax.set_title(f'Set {i}: {label}', fontsize=9)
    axs[0, 0].legend(fontsize=6.5, ncol=2)
    fig.suptitle('Option 2 (exact aging) ÷ exp 30\'s exact-aging reference, by year (solid: mid-step force of infection; dotted: start of step); grey band ±0.2%', fontsize=10)
    fig.tight_layout()
    fig.savefig(here / 'steps.png', dpi=120)
    plt.close(fig)
    print('Saved steps.png')
    return


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'summary'
    dict(copy=run_copy, steps=run_steps, start=run_start, prior=run_prior, summary=summarize)[cmd]()
