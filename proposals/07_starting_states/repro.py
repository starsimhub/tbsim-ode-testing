"""
Reproductions for proposal 07: agents' starting states given by age group in Starsim v4 (starsim d83f4f87). Starsim only. Run from the repo root:

    PYTHONPATH=v4 v4/.venv/bin/python proposals/07_starting_states/repro.py
"""
import os
os.environ.setdefault('NUMBA_NUM_THREADS', '1')
import numpy as np
import pandas as pd
import starsim as ss

ages = list(range(0, 101)) # Single years of age, the last 100+


def shares(sim, states, bins):
    """ Each state's share of the agents in each age band, at the start """
    ppl, d = sim.people, sim.diseases[0]
    u = np.asarray(ppl.auids)
    age = np.asarray(ppl.age.raw[u], dtype=float)
    out = {}
    for k in states:
        m = np.zeros(len(ppl.age.raw), bool)
        m[np.asarray(getattr(d, k).uids)] = True
        m = m[u]
        out[k] = [float(m[(age >= a) & (age < b)].mean()) for a, b in zip(bins[:-1], bins[1:])]
    return out


def make(init, strat, seed=1, n=20_000, disease=ss.SIR):
    sim = ss.Sim(diseases=disease(beta=0, init=init), people=ss.population('India', year=1950, n_agents=n), total_pop=n, start=1950, stop=1951,
                 solver='abm', verbose=0, rand_seed=seed)
    sim.stratify(age=strat)
    sim.init()
    return sim


print('1. Counts by group are rounded to whole agents: a share of 0.3% infected at every age, agents stratified by single years of age')
bins = [0, 15, 50, 65, 80, 101]
for seed in [1, 2, 3]:
    sim = make(dict(I=ss.by(pd.DataFrame(dict(age=ages, value=0.003)))), ages, seed=seed, disease=ss.SIS)
    I = np.asarray(sim.diseases.sis.I.uids)
    age = np.asarray(sim.people.age.raw[I], dtype=float)
    print(f'   seed {seed}: infected aged 0–14, 15–49, 50–64, 65–79, 80+: {[int(((age >= a) & (age < b)).sum()) for a, b in zip(bins[:-1], bins[1:])]}')
sim = make(dict(I=ss.by(pd.DataFrame(dict(age=ages, value=0.003)))), ages, disease=ss.SIS)
age = np.asarray(sim.people.age.raw[sim.people.auids], dtype=float)
print(f'   expected (0.3% of each band): {[round(0.003*((age >= a) & (age < b)).sum(), 1) for a, b in zip(bins[:-1], bins[1:])]}')

print('\n2. A table finer than the agents\' age groups, with two states given: I = 0.1 and R = 0.5 at every age, by single year; agents stratified by 5-year groups (expected S 0.4, I 0.1, R 0.5)')
fine = lambda v: ss.by(pd.DataFrame(dict(age=ages, value=v)))
bins = [0, 5, 15, 30, 50, 70, 101]
sim = make(dict(I=fine(0.1), R=fine(0.5)), list(range(0, 101, 5)))
print('   shares by age band', bins, {k: np.round(v, 3).tolist() for k, v in shares(sim, ['S', 'I', 'R'], bins).items()})
sim = make(dict(I=fine(0.1), R=fine(0.5)), ages)
print('   the same table, agents stratified by single years:', {k: np.round(v, 3).tolist() for k, v in shares(sim, ['S', 'I', 'R'], bins).items()})

print('\n3. The same table with a share for every state, the first included (S = 0.4, I = 0.1, R = 0.5 at every age, by single year)')
for strat, label in [(ages, 'agents stratified by single years'), (list(range(0, 101, 5)), 'agents stratified by 5-year groups')]:
    try:
        make(dict(S=fine(0.4), I=fine(0.1), R=fine(0.5)), strat)
        print(f'   {label}: initialized')
    except Exception as e:
        print(f'   {label}: {type(e).__name__}: {str(e).splitlines()[0][:200]}')
