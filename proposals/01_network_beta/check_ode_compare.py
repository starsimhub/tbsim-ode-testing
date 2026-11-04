"""
The setting of tbsim v4's tests/devtests/ode_compare.py (tbsim.TB as an ODE and as agents on RandomNet(n_contacts=10)), rerun with
two per-contact betas: v4's beta/K, and the matched beta_edge (tbode4.beta_edge_from_ode). Agents with homogeneous mixing are the control.
Run in the v4 environment, from the repo root: v4/.venv/bin/python proposals/01_network_beta/check_ode_compare.py
"""
import os
import numpy as np
import sciris as sc
import starsim as ss
import tbode4

outdir = os.path.dirname(os.path.abspath(__file__))
# As ode_compare.py: tbsim's default beta (9 per year), births replacing every death, 1% symptomatic at the start, 10k agents, 1950-2100
kw = dict(beta=9.0, start=1950, stop=2100, init=dict(symptomatic=0.01), replace_tb_deaths=True, total_pop=10_000)
n_agents = 10_000
window = (2050, 2100) # Prevalence averaged over these years, near equilibrium
arms = {
    'monthly, mass action':    dict(dt=ss.months(1), n_contacts=None, n_runs=40),
    'monthly, network, beta/K': dict(dt=ss.months(1), n_contacts=10, edge_beta='rate', n_runs=40),
    'monthly, network, matched': dict(dt=ss.months(1), n_contacts=10, edge_beta='step', n_runs=40),
    'weekly, network, beta/K':  dict(dt=ss.weeks(1), n_contacts=10, edge_beta='rate', n_runs=20),
    'weekly, network, matched': dict(dt=ss.weeks(1), n_contacts=10, edge_beta='step', n_runs=20),
}


def mean_prev(o):
    """ Mean active (ASY + SYM) prevalence over the window """
    sel = (o.year >= window[0]) & (o.year <= window[1])
    return float(o.prev_active[sel].mean())


def run_arm(arm, seed):
    a = {k: v for k, v in arms[arm].items() if k != 'n_runs'}
    sim = tbode4.make_sim(solver='abm', n_agents=n_agents, seed=seed, **kw, **a).run()
    o = tbode4.shares(sim)
    return dict(arm=arm, seed=seed, prev=mean_prev(o), year=o.year, prev_active=o.prev_active)


if __name__ == '__main__':
    T = sc.timer()
    ode = tbode4.shares(tbode4.make_sim(dt=ss.months(1), **kw).run())
    ode_prev = mean_prev(ode)
    iterkwargs = [dict(arm=arm, seed=s) for arm, a in arms.items() for s in range(a['n_runs'])]
    outs = sc.parallelize(run_arm, iterkwargs=iterkwargs, ncpus=tbode4.ncpus)
    summary = sc.objdict(ode_prev=ode_prev)
    print(f'ODE prevalence {window[0]}-{window[1]}: {ode_prev*1e5:.1f}/100k')
    for arm in arms:
        p = np.array([o['prev'] for o in outs if o['arm'] == arm])/ode_prev
        summary[arm] = dict(mean=p.mean(), se=p.std(ddof=1)/np.sqrt(len(p)), n_runs=len(p), died_out=int((p == 0).sum()))
        print(f'{arm:28s} agents ÷ ODE {p.mean():.3f} ± {summary[arm]["se"]:.3f} ({len(p)} runs, {summary[arm]["died_out"]} died out)')
    T.toc()
    sc.save(os.path.join(outdir, 'check_ode_compare.obj'), dict(ode=ode, outs=outs, summary=summary))
    sc.savejson(os.path.join(outdir, 'check_ode_compare.json'), sc.jsonify(summary))
    tbode4.save_versions(outdir, 'v4')
