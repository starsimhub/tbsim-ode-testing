"""
Exp 26: do agents and the ODE agree on TB deaths on treatment when treatment varies in time (ss.timeseries)?

In the full model, agents' TB deaths on treatment read 0.983 ± 0.006 of the ODE per prevalent case through the scale-up, against
0.994 ± 0.006 at constant treatment. This isolates treatment: people in Y start treatment at tb_LAI_TPT's time-varying rate, a
time-varying share start courses that will cure (Ts), the rest courses that will fail (Tf); failing courses face TB death; every
course ends after 6 months, back to Y. Arms: the scale-up's timeseries ('ramp'), and 2024's values held constant ('const');
each with a fixed course or a rate, at monthly steps. Agents ÷ ODE over 1995–2035 for starts, failing starts, person-time on
failing courses, and TB deaths on treatment. Run from the repo root:

    PYTHONPATH=v4 v4/.venv/bin/python experiments/26_v4_pin_bump/timeseries_check.py
"""
import os
os.environ.setdefault('NUMBA_NUM_THREADS', '1')
import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss
import tbode4

outdir = sc.thispath(__file__)
sym_dead, D = 0.3959, 0.5 # Exp 21 set 1's sym_dead; course length (years)
n_agents, seeds = 400_000, range(8)
start, stop, window = 1980, 2036, (1995, 2036)


def run(scen, fixed, solver, seed=0):
    lai = tbode4.lai_tpt_treatment()
    if scen == 'ramp':
        rate, success = lai.tx_rate, lai.tx_success
    else:
        rate, success = ss.peryear(float(tbode4.theta_from_coverage(lai.coverage(2024)))), float(lai.efficacy(2024))
    Y, Ts, Tf = ss.states(3)
    dur = ss.constant(ss.years(D)) if fixed else ss.years(D)
    dis = ss.Disease(name='tx', Y=Y, Ts=Ts, Tf=Tf, init=dict(Y=1.0),
                     start=ss.Transition(Y, {Ts: success, Tf: None}, rate=rate),
                     cure=ss.Transition(Ts, Y, dur=dur),
                     fail=ss.Transition(Tf, Y, dur=dur),
                     death_tx=ss.death(Tf, rate=ss.peryear(sym_dead)))
    kw = dict(n_agents=n_agents, rand_seed=seed) if solver == 'abm' else dict(total_pop=1)
    sim = ss.Sim(dis, solver=solver, start=start, stop=stop, dt=ss.months(1), verbose=0, **kw).run()
    r = sim.results.tx
    v = lambda k: np.asarray(r[k].values, dtype=float)
    year = np.asarray(sim.t.yearvec, dtype=float)
    alive = v('Y') + v('Ts') + v('Tf')
    flow_in = (year >= window[0]) & (year < window[1]) # Flows at index i are the step from t_i (exp 26, finding 1)
    stock_in = (year > window[0]) & (year <= window[1]) # Stocks at the end of each of those steps
    return dict(scen=scen, fixed=fixed, solver=solver, seed=seed,
                starts=(v('start')/alive)[flow_in].sum(),
                on_fail=(v('Tf')/alive)[stock_in].mean(),
                deaths_tx=(v('death_tx')/alive)[flow_in].sum(),
                fails=(v('fail')/alive)[flow_in].sum())


if __name__ == '__main__':
    jobs = []
    for scen in ['ramp', 'const']:
        for fixed in [True, False]:
            jobs += [dict(scen=scen, fixed=fixed, solver='abm', seed=s) for s in seeds]
            jobs.append(dict(scen=scen, fixed=fixed, solver='ode'))
    df = pd.DataFrame(sc.parallelize(run, iterkwargs=jobs, ncpus=tbode4.ncpus))
    df['deaths_per_on_fail'] = df.deaths_tx/df.on_fail
    df['share_dying'] = df.deaths_tx/(df.deaths_tx + df.fails)
    df.to_csv(outdir / 'timeseries_check_runs.csv', index=False)
    cols = ['starts', 'on_fail', 'deaths_tx', 'deaths_per_on_fail', 'share_dying']
    rows = []
    for (scen, fixed), g in df.groupby(['scen', 'fixed']):
        ode = g[g.solver == 'ode'].iloc[0]
        a = g[g.solver == 'abm']
        row = dict(scen=scen, course='fixed' if fixed else 'rate', reps=len(a))
        for c in cols:
            ratio = a[c]/ode[c]
            row[c] = ratio.mean(); row[c + '_se'] = ratio.sem()
        rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(outdir / 'timeseries_check.csv', index=False)
    with pd.option_context('display.width', 250, 'display.precision', 4):
        print(f'Agents ÷ ODE, {window[0]}–{window[1] - 1}, monthly steps, {n_agents:,} agents (mean ± SE over {len(seeds)} seeds):')
        print(out.to_string(index=False))
