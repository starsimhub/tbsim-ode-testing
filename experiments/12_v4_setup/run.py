"""
Experiment 12, v4 side: does tbsim v4's ODE reproduce TBODE, and how should agents be matched to it? Run in the v4 environment. See README.md.
"""
import os
import io
import json
import time
import contextlib
import numpy as np
import sciris as sc
import starsim as ss
import tbsim
import tbsim.compartmental as tbc
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import tbode4

outdir = os.path.dirname(os.path.abspath(__file__))
ref = sc.loadjson(os.path.join(outdir, 'reference_v010.json'))
states = tbode4.states
K = 10 # Mean degree of the network (n_contacts)
n_agents = 50_000
tol = 1e-4 # Criterion for matching the ODEs (relative)
dts = {'1 month': ss.months(1), '1 week': ss.weeks(1), '1 day': ss.days(1)}


def max_rel(a, b, floor=1e-8):
    """ The largest relative difference between two arrays, where the reference is above a floor """
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    ok = np.abs(b) > floor
    return float(np.max(np.abs(a[ok] - b[ok])/np.abs(b[ok]))) if ok.any() else 0.0


def eq_init():
    """ The share in each state at TBODE's baseline equilibrium (the empty treatment states come out as ±1e-35 from the root solve: set to 0) """
    return {s: max(ref['equilibrium_A'][tbode4.state_map[s]], 0.0) for s in states}


# A1: the reference ODE
def check_reference():
    """ v4 TB as an ODE, with births replacing every death, against TB_ODE over 520 years (run from 1900, since the model has no calendar time) """
    tb_ode = tbc.TB_ODE().run(start_time=1500, end_time=2020)
    sim = tbode4.make_sim(beta=9.0, start=1900, stop=2420, dt=ss.years(1), init=dict(symptomatic=1e3/1e5), replace_tb_deaths=True).run()
    r = sim.results.tb
    errs = {s: max_rel(np.asarray(r[s].values, dtype=float), tb_ode[s.upper()], floor=1e-2) for s in states}
    return dict(errs=errs, max=max(errs.values()), years=tb_ode.time.tolist(), v4={s: np.asarray(r[s].values, dtype=float).tolist() for s in states},
                tb_ode={s: tb_ode[s.upper()].tolist() for s in states})


# A2: our baseline over time
def check_trajectory(label):
    """ v4 ODE from a 1% symptomatic seed against TBODE's saved trajectory (shares), monthly, 1900-2100 """
    beta = ref[f'beta_{label}']
    pars = tbode4.odepars_to_tb(ref['pars_B']) if label == 'B' else {}
    sim = tbode4.make_sim(beta=beta, start=1900, stop=2100, dt=ss.months(1), init=dict(symptomatic=0.01), **pars).run()
    v4 = tbode4.shares(sim)
    tr = ref[f'trajectory_{label}']
    assert len(v4.year) == len(tr['year']) and np.allclose(v4.year, tr['year']), 'Time points differ'
    later = v4.year >= 1901 # After the first year: in the first months, compartments filling from 0 carry the integrators' first-step error
    errs_all = {s: max_rel(v4[s], tr[tbode4.state_map[s]]) for s in states}
    errs = {s: max_rel(v4[s][later], np.array(tr[tbode4.state_map[s]])[later]) for s in states}
    abs_max = max(float(np.max(np.abs(v4[s] - np.array(tr[tbode4.state_map[s]])))) for s in states)
    prev_v4 = v4.asymptomatic + v4.symptomatic
    prev_ref = np.array(tr['ASYMPTOMATIC']) + np.array(tr['SYMPTOMATIC'])
    return dict(errs=errs, max=max(errs.values()), errs_all=errs_all, max_all=max(errs_all.values()), abs_max=abs_max, year=v4.year.tolist(),
                prev_v4=prev_v4.tolist(), prev_ref=prev_ref.tolist(), rel_prev=((prev_v4 - prev_ref)/prev_ref).tolist())


# A3: the calibration
def check_stationary():
    """ Started at TBODE's equilibrium, the v4 ODE's prevalence over 200 years """
    sim = tbode4.make_sim(beta=ref['beta_A'], start=1900, stop=2100, dt=ss.months(1), init=eq_init()).run()
    v4 = tbode4.shares(sim)
    errs = {s: max_rel(v4[s], ref['equilibrium_A'][tbode4.state_map[s]]) for s in states}
    return dict(errs=errs, max=max(errs.values()), year=v4.year.tolist(), prev=v4.prev_active.tolist())


def v4_eq_prev(beta, years=3000):
    """ The v4 ODE's prevalence after a long run from a 1% symptomatic seed (yearly output), and its change over the last 100 years """
    sim = tbode4.make_sim(beta=beta, start=1900, stop=1900 + years, dt=ss.years(1), init=dict(symptomatic=0.01)).run()
    prev = tbode4.shares(sim).prev_active
    return prev[-1], prev[-1] - prev[-101]


def calibrate_v4(target=500e-5):
    """ beta for 500/100k on the v4 ODE alone """
    beta = brentq(lambda b: v4_eq_prev(b)[0] - target, 9.5, 10.5, xtol=1e-6)
    prev, drift = v4_eq_prev(beta)
    return dict(beta=beta, beta_tbode=ref['beta_A'], rel_diff=beta/ref['beta_A'] - 1, prev=prev, drift_last_100y=drift)


# B1: network degree
class EdgeProbe(ss.Analyzer):
    """ Keeps the network's edges at each of the first few steps """
    def __init__(self, n_steps=3, **kwargs):
        super().__init__(**kwargs)
        self.n_steps = n_steps
        self.snaps = []

    def step(self):
        if len(self.snaps) < self.n_steps:
            net = self.sim.networks[0]
            self.snaps.append((np.array(net.edges.p1), np.array(net.edges.p2)))


def check_network():
    """ Mean degree, its variance, edges per agent, and the share of edges kept from one step to the next """
    sim = tbode4.make_sim(beta=ref['beta_A'], solver='abm', n_agents=n_agents, n_contacts=K, init=eq_init(), stop=1901, analyzers=EdgeProbe()).run()
    probe = sim.analyzers[0] # The sim's copy
    out = []
    for p1, p2 in probe.snaps:
        deg = np.bincount(np.concatenate([p1, p2]), minlength=n_agents)
        out.append(dict(edges_per_agent=len(p1)/n_agents, mean_degree=float(deg.mean()), var_degree=float(deg.var())))
    e0 = set(zip(*probe.snaps[0]))
    e1 = set(zip(*probe.snaps[1]))
    return dict(steps=out, kept=len(e0 & e1)/len(e0))


# B2: per-contact beta
def check_network_ode():
    """ The network sim as an ODE (homogeneous mixing at its mean degree), with beta/K per contact, against the plain ODE """
    kw = dict(beta=ref['beta_A'], start=1900, stop=2100, dt=ss.months(1), init=dict(symptomatic=0.01))
    plain = tbode4.shares(tbode4.make_sim(**kw).run())
    sim = tbode4.make_sim(n_contacts=K, edge_beta='rate', **kw)
    with contextlib.redirect_stdout(io.StringIO()):
        sim.run()
    net = tbode4.shares(sim)
    errs = {s: max_rel(net[s], plain[s]) for s in states}
    return dict(errs=errs, max=max(errs.values()), notes=[str(n) for n in sim.notes])


arms = {'mass action': dict(n_contacts=None), 'network, beta/K': dict(n_contacts=K, edge_beta='rate'), 'network, matched': dict(n_contacts=K, edge_beta='step')}


def run_transmission(arm, seed, dt=ss.months(1), years=50):
    """ Infections in each step as agents, against what the ODE's force of infection predicts from the state at the start of the step """
    sim = tbode4.make_sim(beta=ref['beta_A'], solver='abm', n_agents=n_agents, init=eq_init(), start=1900, stop=1900 + years, dt=dt, seed=seed, **arms[arm]).run()
    r = sim.results.tb
    get = lambda s: np.asarray(r[s].values, dtype=float)
    N = sum(get(s) for s in states)
    p = tbsim.tb_pars()
    w = (p.trans_asymp*get('asymptomatic') + get('symptomatic'))/N # Infectious weight at each time point
    x = ref['beta_A']*w*sim.t.dt.years # The ODE's infections per susceptible in a step, at the start-of-step state
    S = get('susceptible')
    obs = get('transmission')[1:] # The flow during the step ending at each time point
    return dict(arm=arm, seed=seed, obs=obs.sum(), linear=(S[:-1]*x[:-1]).sum(), exp=(S[:-1]*-np.expm1(-x[:-1])).sum(),
                scale=float(N[0]/n_agents), mean_x=float(x[:-1].mean()), prev=(get('asymptomatic') + get('symptomatic'))/N)


def check_transmission(n_seeds=4):
    """ Pooled over seeds: observed ÷ expected (linear ODE flow, and 1 − e^(−x)), with Poisson SEs """
    iterkwargs = [dict(arm=arm, seed=s) for arm in arms for s in range(n_seeds)]
    outs = sc.parallelize(run_transmission, iterkwargs=iterkwargs, ncpus=tbode4.ncpus)
    summary = {}
    for arm in arms:
        rows = [o for o in outs if o['arm'] == arm]
        obs = sum(o['obs'] for o in rows)
        n_events = obs/rows[0]['scale'] # Agents infected, for the Poisson SE
        se = 1/np.sqrt(n_events)
        lin = sum(o['linear'] for o in rows)
        ex = sum(o['exp'] for o in rows)
        summary[arm] = dict(obs_over_linear=obs/lin, obs_over_exp=obs/ex, rel_se=se, n_infected_agents=n_events,
                                                            mean_x=float(np.mean([o['mean_x'] for o in rows])))
    return dict(summary=summary, runs=outs)


# B3: starting at the equilibrium
def check_init():
    """ The share in each state at the start, as agents, against the shares given to init= """
    sim = tbode4.make_sim(beta=ref['beta_A'], solver='abm', n_agents=n_agents, init=eq_init(), stop=1901).run()
    v4 = tbode4.shares(sim)
    init = eq_init()
    diffs = {s: float(v4[s][0] - init[s]) for s in states}
    return dict(diffs=diffs, max_abs_agents=max(abs(d) for d in diffs.values())*n_agents)


# B4: notes and run time
def check_dt_notes():
    """ For each dt: the notes v4 gives as agents, and the time per step at 50k agents (homogeneous mixing and the network) """
    out = {}
    for lab, dt in dts.items():
        years = {'1 month': 20, '1 week': 5, '1 day': 2}[lab]
        row = {}
        for net in [False, True]:
            sim = tbode4.make_sim(beta=ref['beta_A'], solver='abm', n_agents=n_agents, n_contacts=K if net else None, init=eq_init(), stop=1900 + years, dt=dt)
            sim.init()
            T = time.perf_counter()
            sim.run()
            elapsed = time.perf_counter() - T
            row['network' if net else 'mass action'] = dict(ms_per_step=1e3*elapsed/(len(sim.t.yearvec) - 1), notes=[str(n) for n in sim.notes])
        out[lab] = row
    return out


def run_check_dt():
    """ sim.check_dt() at monthly dt (homogeneous mixing, 50 years): its printed table, kept as text """
    sim = tbode4.make_sim(beta=ref['beta_A'], solver='abm', n_agents=n_agents, init=eq_init(), stop=1950).run()
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        sim.check_dt()
    return buf.getvalue()


def plot(res):
    """ A2: relative difference in prevalence over time (v4 ODE ÷ TBODE − 1); B2: observed ÷ expected infections """
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.5))
    for label in ['A', 'B']:
        t = res[f'A2_{label}']
        axs[0].plot(t['year'], t['rel_prev'], label=f'{"Baseline (tbsim defaults)" if label == "A" else "Set B natural history, k = 0"}')
    axs[0].axhline(0, color='k', lw=0.5)
    axs[0].set_xlabel('Year')
    axs[0].set_ylabel('v4 ODE ÷ TBODE − 1 (prevalence)')
    axs[0].set_title('A2: v4 ODE vs TBODE, from a 1% symptomatic seed')
    axs[0].legend(frameon=False)
    s = res['B2_transmission']['summary']
    labs = list(s)
    for i, key in enumerate(['obs_over_linear', 'obs_over_exp']):
        vals = [s[l][key] for l in labs]
        errs = [s[l]['rel_se']*v for l, v in zip(labs, vals)]
        axs[1].errorbar(np.arange(len(labs)) + 0.15*(2*i - 1), vals, yerr=errs, fmt='o', capsize=4,
                        label='÷ ODE flow (x)' if i == 0 else '÷ 1 − e^(−x)')
    axs[1].axhline(1, color='k', lw=0.5)
    pred = -np.expm1(-ref['beta_A']/K/12)/(ref['beta_A']/K/12) # Per-edge saturation with beta/K: (1 − e^(−b·dt))/(b·dt)
    axs[1].plot([0.7, 1.3], [pred, pred], color='C3', lw=1, ls='--', label='predicted for beta/K: (1 − e^(−b·dt))/(b·dt)')
    axs[1].set_xticks(range(len(labs)), labs)
    axs[1].set_ylabel('Infections as agents ÷ expected (± SE)')
    axs[1].set_title(f'B2: monthly dt, {n_agents//1000}k agents, 4 seeds × 50 years')
    axs[1].legend(frameon=False)
    sc.figlayout()
    fig.savefig(os.path.join(outdir, 'v4_setup.png'), dpi=120)
    return fig


if __name__ == '__main__':
    T = sc.timer()
    res = sc.objdict()
    res.A1 = check_reference()
    print(f'A1 v4 ODE vs TB_ODE: max relative difference {res.A1["max"]:.2e}')
    for label in ['A', 'B']:
        res[f'A2_{label}'] = check_trajectory(label)
        a2 = res[f'A2_{label}']
        print(f'A2 {label} v4 ODE vs TBODE: max relative difference {a2["max"]:.2e} from 1901 ({a2["max_all"]:.2e} including the first year); max absolute {a2["abs_max"]:.1e}')
    res.A3_stationary = check_stationary()
    print(f'A3 stationary at TBODE equilibrium: max relative drift {res.A3_stationary["max"]:.2e}')
    res.A3_calibration = calibrate_v4()
    print(f'A3 v4 calibration: beta = {res.A3_calibration["beta"]:.6f} (TBODE {ref["beta_A"]:.6f}; {res.A3_calibration["rel_diff"]:.1e})')
    T.toc('A')
    res.B1 = check_network()
    print(f'B1 network: {res.B1["steps"][0]}, edges kept between steps: {res.B1["kept"]:.3f}')
    res.B2_ode = check_network_ode()
    print(f'B2 network as an ODE vs plain ODE: max relative difference {res.B2_ode["max"]:.2e}; notes: {res.B2_ode["notes"]}')
    res.B2_transmission = check_transmission()
    print(f'B2 transmission: {sc.jsonify(res.B2_transmission["summary"])}')
    res.B3 = check_init()
    print(f'B3 init: largest difference {res.B3["max_abs_agents"]:.2f} agents')
    res.B4 = check_dt_notes()
    for lab, row in res.B4.items():
        print(f'B4 {lab}: ' + ', '.join(f'{k} {v["ms_per_step"]:.1f} ms/step' for k, v in row.items()))
    res.B4_check_dt = run_check_dt()
    T.toc('B')

    sc.save(os.path.join(outdir, 'outputs.obj'), res)
    summary = {k: v for k, v in res.items() if k not in ('A1', 'A2_A', 'A2_B', 'A3_stationary', 'B2_transmission')}
    summary.update(A1_max=res.A1['max'], A2_A_max=res.A2_A['max'], A2_B_max=res.A2_B['max'], A3_stationary_max=res.A3_stationary['max'],
                   B2_transmission=res.B2_transmission['summary'])
    sc.savejson(os.path.join(outdir, 'summary.json'), sc.jsonify(summary))
    plot(res)
    tbode4.save_versions(outdir, 'v4')
    os.replace(os.path.join(outdir, 'versions.json'), os.path.join(outdir, 'versions_v4.json'))
