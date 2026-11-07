"""
Experiment 18: Starsim v4's ODE of a network (starsim 5ba36ebb) against agents on the same network, and against our convention
(mass-action ODE, agents with the matched per-contact beta). See README.md. Run from the repo root, in this experiment's environment:

    PYTHONPATH=v4 experiments/18_network_ode/.venv/bin/python experiments/18_network_ode/run.py q3   # calibration (needed by q1, q2)
    PYTHONPATH=v4 experiments/18_network_ode/.venv/bin/python experiments/18_network_ode/run.py q1   # one step, per transition
    PYTHONPATH=v4 experiments/18_network_ode/.venv/bin/python experiments/18_network_ode/run.py q2   # equilibrium, 40 reps per arm

Rerunning resumes: parts already saved are kept.
"""
import sys
import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss
from scipy.optimize import brentq

outdir = sc.thispath(__file__)
root = outdir.parent.parent
sys.path.insert(0, str(root / 'v4'))
sys.path.insert(0, str(root / 'experiments' / '17_latent_reinfection'))
import common as c # noqa: E402 (Cliff's v4 helpers)
import run_equilibrium as req # noqa: E402
import tbode4 # noqa: E402
import run_v4 as rv17 # noqa: E402 (exp 17: the two-state model's parameters)

K = 10 # Mean degree of ss.RandomNet (fresh contacts every step)
dts = sc.objdict(month=ss.months(1), week=ss.weeks(1))
sets = ['A', 'B2'] # A: tbsim's defaults; B2: exp 17's two-state set B at phi = 0.01, every reinfection risk at pi
target = 0.005 # Active (asymptomatic + symptomatic) prevalence
n_agents, reps, seed0 = 50_000, 40, 18000
start, stop = 1900, 2000
fn = sc.objdict(q1=outdir / 'q1.obj', q2=outdir / 'q2.obj', q3=outdir / 'q3.obj')
contact_flows = ['transmission', 'reinfect_cle', 'reinfect_rec', 'reinfect_trd', 'reinfect_remote']


def set_kw(name):
    """ The parameter set and extra tbsim.TB parameters of each named set """
    if name == 'A':
        return 'A', {}
    pars = {k: v for k, v in rv17.two_state_pars(0.01).items() if k != 'fit'}
    return 'B', pars


def rr_of(name, flow):
    """ Each transition by contact's risk relative to a susceptible's """
    parset, pars = set_kw(name)
    p = sc.mergedicts(c.parsets[parset], pars)
    import tbsim
    dflt = tbsim.tb_pars()
    get = lambda k: float(p[k]) if k in p else float(dflt[k].raw if hasattr(dflt[k], 'raw') else dflt[k])
    return dict(transmission=1.0, reinfect_cle=get('rr_reinfection_cleared'), reinfect_rec=get('rr_reinfection_rec'),
                reinfect_trd=get('rr_reinfection_treat'), reinfect_remote=get('rr_reinfection_inf'))[flow]


def make_sim(name, beta, solver, dt, network, init, n=None, seed=0, start=start, stop=stop):
    """ A sim of one set: `beta` is per year (mass action) or per contact (network) """
    parset, pars = set_kw(name)
    tb = c.make_tb(parset, beta, init=init, model='tb', **pars)
    kw = dict(start=start, stop=stop, dt=dt, total_pop=n or 1e5, rand_seed=seed, verbose=0)
    if network:
        kw['networks'] = ss.RandomNet(n_contacts=ss.constant(K), dur=0)
    if n:
        kw['n_agents'] = n
    return ss.Sim(tb, solver=solver, **kw)


def ode_equilibrium(name, beta, dt, network, init, years=600, tol=1e-8, maxiter=20):
    """ The ODE's equilibrium shares at this step (the network ODE's correction depends on it), run from `init` until they settle """
    shares = init
    for i in range(maxiter):
        sim = make_sim(name, beta, 'ode', dt, network, shares, start=0, stop=years)
        sim.run()
        new = c.end_shares(sim)
        change = max(abs(new[k] - shares.get(k, 0)) for k in new)
        shares = new
        if change < tol:
            break
    return shares


def q3():
    """ beta for 500/100k: mass-action ODE (any step), and network ODE (per contact, at each step); with the matched beta for comparison """
    data = sc.load(fn.q3) if fn.q3.exists() else sc.objdict()
    for name in sets:
        parset, pars = set_kw(name)
        if f'{name}_mass' not in data:
            T = sc.timer()
            cal = c.calibrate(parset, model='tb', **pars)
            data[f'{name}_mass'] = sc.objdict(beta=cal.beta, equilibrium=cal.equilibrium, prevalence=cal.prevalence)
            sc.save(fn.q3, data)
            print(f'{name} mass action: beta = {cal.beta:.5f} ({T.toc(output=True):.0f} s)', flush=True)
        ma = data[f'{name}_mass']
        for dkey, dt in dts.items():
            key = f'{name}_net_{dkey}'
            if key in data:
                continue
            T = sc.timer()
            matched = tbode4.beta_edge_from_ode(ma.beta, K, dt.years)
            eq = {}
            def gap(b):
                eq[b] = ode_equilibrium(name, b, dt, True, ma.equilibrium)
                return c.prevalence(eq[b]) - target
            b = brentq(gap, 0.97*matched, 1.03*matched, xtol=1e-7*matched)
            shares = eq[b] if b in eq else ode_equilibrium(name, b, dt, True, ma.equilibrium)
            data[key] = sc.objdict(beta_edge=b, beta_eff=tbode4.beta_ode_from_edge(b, K, dt.years), matched=matched, equilibrium=shares,
                                   prevalence=c.prevalence(shares), prev_at_matched=c.prevalence(ode_equilibrium(name, matched, dt, True, ma.equilibrium)))
            sc.save(fn.q3, data)
            d = data[key]
            print(f'{key}: per-contact b = {b:.5f} (matched {matched:.5f}, ratio {b/matched:.5f}); K·(1 − e^(−b·dt))/dt = {d.beta_eff:.5f} '
                  f'vs mass-action beta {ma.beta:.5f}; network ODE prevalence at the matched b: {1e5*d.prev_at_matched:.2f}/100k ({T.toc(output=True):.0f} s)', flush=True)
    tbode4.save_versions(outdir, 'v4 (exp 18 environment)')
    return data


def flows_over(sim, nsteps):
    """ Each transition by contact's total over the first nsteps steps """
    r = sim.results.tb
    pop = float(sim.results.alive.values[0]) # Per person, since the agents and the ODE stand for different population sizes
    return {f: float(np.sum(r[f].values[:nsteps]))/pop for f in contact_flows if f in r}


def q1_one(name, arm, dkey, seed, n=1_000_000, nsteps=24):
    """ One arm of question 1: agents over nsteps steps from the ODE's equilibrium, and the matching ODE """
    cal = sc.load(fn.q3)
    dt = dts[dkey]
    stop_ = start + nsteps*dt.years
    if arm == 'network':
        net = cal[f'{name}_net_{dkey}']
        b, eq, ode = net.beta_edge, net.equilibrium, (net.beta_edge, True)
    else: # Our convention: agents with the matched beta on the network, against the mass-action ODE
        ma = cal[f'{name}_mass']
        b, eq, ode = tbode4.beta_edge_from_ode(ma.beta, K, dt.years), ma.equilibrium, (ma.beta, False)
    sim = make_sim(name, b, 'abm', dt, True, eq, n=n, seed=seed, start=start, stop=stop_)
    sim.run()
    out = sc.objdict(agents=flows_over(sim, nsteps), b=b)
    if seed == seed0:
        osim = make_sim(name, ode[0], 'ode', dt, ode[1], eq, start=start, stop=stop_)
        osim.run()
        out.ode = flows_over(osim, nsteps)
    return out


def q1():
    data = sc.load(fn.q1) if fn.q1.exists() else sc.objdict()
    seeds = [seed0 + s for s in range(4)]
    for name in sets:
        for dkey in dts:
            for arm in ['network', 'convention']:
                key = f'{name}_{dkey}_{arm}'
                if key in data:
                    continue
                T = sc.timer()
                res = sc.parallelize(q1_one, iterkwargs=[dict(name=name, arm=arm, dkey=dkey, seed=s) for s in seeds], ncpus=2) # 1M agents each: 2 at a time
                data[key] = sc.objdict(runs=res, ode=res[0].ode, b=res[0].b)
                sc.save(fn.q1, data)
                print(f'{key} ({T.toc(output=True):.0f} s)', flush=True)
    summarise_q1(data)
    return data


def summarise_q1(data):
    """ Agents ÷ ODE for each transition by contact, with the prediction from each edge's per-step probability """
    cal = sc.load(fn.q3)
    rows = []
    for key, d in data.items():
        name, dkey, arm = key.split('_')
        D = dts[dkey].years
        for f, ode_val in d.ode.items():
            if ode_val <= 0:
                continue
            ag = np.array([r.agents[f] for r in d.runs])
            rr = rr_of(name, f)
            x = rr*d.b*D
            agent_factor = -np.expm1(-x)/x # Agents' per-contact share of mass action, for this transition's risk
            if arm == 'network': # The ODE applies the susceptibles' factor to every transition by contact
                ode_factor = -np.expm1(-d.b*D)/(d.b*D)
            else: # The mass-action ODE with beta = K·(1 − e^(−b·D))/D: the susceptibles' factor, as the ODE's own
                ode_factor = -np.expm1(-d.b*D)/(d.b*D)
            rows.append(dict(set=name, dt=dkey, arm=arm, flow=f, rr=rr, agents_over_ode=ag.mean()/ode_val, se=ag.std(ddof=1)/np.sqrt(len(ag))/ode_val,
                             predicted=agent_factor/ode_factor, events_per_run=ag.mean()*1e6))
    df = pd.DataFrame(rows)
    df.to_csv(outdir / 'q1_summary.csv', index=False)
    with pd.option_context('display.width', 200, 'display.max_rows', 200):
        print(df.round(4).to_string(index=False))
    return df


def outputs(sim):
    """ Yearly means of the shares compared: active prevalence, latent, cleared (all three no-longer-infected states), TB deaths """
    r = sim.results.tb
    alive = sim.results.alive.values
    get = lambda s: r[s].values if s in r else 0
    out = sc.objdict(year=np.array(sim.t.yearvec))
    out.prevalence = (get('asymptomatic') + get('symptomatic'))/alive
    out.latent = (get('infection') + get('infection_remote'))/alive
    out.cleared = (get('cleared') + get('recovered') + get('treated'))/alive
    out.tb_deaths = get('tb_death')/alive/sim.t.dt_year
    return req.yearly(out)


def q2_one(name, arm, dkey, seed, solver='abm'):
    cal = sc.load(fn.q3)
    dt = dts[dkey]
    if arm == 'network':
        net = cal[f'{name}_net_{dkey}']
        b, eq, ode = net.beta_edge, net.equilibrium, (net.beta_edge, True)
    else:
        ma = cal[f'{name}_mass']
        b, eq, ode = tbode4.beta_edge_from_ode(ma.beta, K, dt.years), ma.equilibrium, (ma.beta, False)
    if solver == 'ode':
        sim = make_sim(name, ode[0], 'ode', dt, ode[1], eq)
    else:
        sim = make_sim(name, b, 'abm', dt, True, eq, n=n_agents, seed=seed)
    sim.run()
    return outputs(sim)


def ratio(runs, ode, key, years=(1950, 2000)):
    sel = lambda o: (o.year >= years[0]) & (o.year < years[1])
    return np.array([o[key][sel(o)].mean() for o in runs])/ode[key][sel(ode)].mean()


def q2():
    data = sc.load(fn.q2) if fn.q2.exists() else sc.objdict(ode=sc.objdict(), abm=sc.objdict())
    for name in sets:
        for dkey in dts:
            for arm in ['network']: # Our convention gives the same agents and an identical ODE (q3: the calibrated per-contact betas are equal)
                key = f'{name}_{dkey}_{arm}'
                if key not in data.ode:
                    data.ode[key] = q2_one(name, arm, dkey, 0, solver='ode')
                if key in data.abm and len(data.abm[key]) >= reps:
                    continue
                T = sc.timer()
                data.abm[key] = sc.parallelize(q2_one, iterkwargs=[dict(name=name, arm=arm, dkey=dkey, seed=seed0 + s) for s in range(reps)], ncpus=tbode4.ncpus)
                sc.save(fn.q2, data)
                r = ratio(data.abm[key], data.ode[key], 'prevalence')
                print(f'{key:22s} {reps} reps in {T.toc(output=True):.0f} s; prevalence agents ÷ ODE {r.mean():.3f} ± {r.std(ddof=1)/np.sqrt(reps):.3f}', flush=True)
    summarise_q2(data)
    return data


def summarise_q2(data):
    rows = []
    for key, runs in data.abm.items():
        name, dkey, arm = key.split('_')
        for out in ['prevalence', 'latent', 'cleared', 'tb_deaths']:
            r = ratio(runs, data.ode[key], out)
            rows.append(dict(set=name, dt=dkey, arm=arm, output=out, reps=len(r), mean=r.mean(), se=r.std(ddof=1)/np.sqrt(len(r)),
                             p10=np.percentile(r, 10), p90=np.percentile(r, 90)))
    df = pd.DataFrame(rows)
    df.to_csv(outdir / 'q2_summary.csv', index=False)
    with pd.option_context('display.width', 200, 'display.max_rows', 200):
        print(df.round(4).to_string(index=False))
    return df


def explain():
    """
    The equilibrium gap predicted from the reinfection factor alone: the mass-action ODE at the calibrated beta, with each
    reinfection risk rr scaled by agents' per-contact factor for it, (1 − e^(−rr·b·D))/(rr·b·D), over the factor the ODE applies to
    every transition by contact (the susceptibles', rr = 1). Writes explain.json
    """
    cal = sc.load(fn.q3)
    out = {}
    keys = dict(rr_reinfection_cleared='reinfect_cle', rr_reinfection_rec='reinfect_rec', rr_reinfection_treat='reinfect_trd', rr_reinfection_inf='reinfect_remote')
    for name in sets:
        parset, pars = set_kw(name)
        for dkey, dt in dts.items():
            D = dt.years
            b = cal[f'{name}_net_{dkey}'].beta_edge
            ode_f = -np.expm1(-b*D)/(b*D)
            rr = {k: rr_of(name, f) for k, f in keys.items() if not (k == 'rr_reinfection_inf' and name == 'A')} # Set A has no remote state
            mult = {k: (-np.expm1(-v*b*D)/(v*b*D))/ode_f for k, v in rr.items()}
            eq = c.equilibrium(parset, cal[f'{name}_mass'].beta, model='tb', **sc.mergedicts(pars, {k: rr[k]*mult[k] for k in rr}))
            out[f'{name}_{dkey}'] = dict(per_contact_b=b, ode_factor=ode_f, rr=rr, agents_over_ode_factor=mult, predicted_prevalence_ratio=c.prevalence(eq)/target)
            print(f'{name} {dkey}: predicted agents ÷ ODE prevalence {c.prevalence(eq)/target:.3f} (reinfection multipliers {dict((k, round(v, 4)) for k, v in mult.items())})', flush=True)
    sc.savejson(outdir / 'explain.json', out)
    return out


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'q3'
    dict(q1=q1, q2=q2, q3=q3, explain=explain)[which]()
