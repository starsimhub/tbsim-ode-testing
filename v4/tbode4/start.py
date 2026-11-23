"""
Starting TB states by age and sex, the same in v4's agents and its ODE (CONVENTIONS.md 10; exp 30).

- `use_drawn_start()`: agents' starting counts in each group drawn (multinomial), not rounded to whole agents, which loses rare states in small groups
- `settled_start()`: each state's share by age and sex from a run of the ODE (e.g. with the starting year's demography frozen until it settles),
  averaged to the agents' age groups, as `init=` tables for both solvers
- `shares_by_group()`, `coarsen()`, `init_table()`: its parts

Moved here from exp 30's run.py, which keeps its own copy so it can be rerun as it was. Starsim behaviours worked around: proposal 07.
"""
import contextlib
import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss
from .demog import sexes, age_edges

__all__ = ['tb_states', 'drawn_counts', 'use_drawn_start', 'shares_by_group', 'coarsen', 'init_table', 'settled_start']

# TBTreat's states (tbsim.TB's, and treatment's)
tb_states = ['susceptible', 'infection', 'infection_remote', 'cleared', 'recovered', 'non_infectious', 'asymptomatic', 'symptomatic', 'treated',
             'treatment_success', 'treatment_fail', 'relapsing']


#%% Agents' starting counts, drawn

@contextlib.contextmanager
def drawn_counts(rng):
    """ Within the block, Starsim's exact counts by group (Compartmental.round_preserving) are drawn from a multinomial with rng instead """
    Comp = ss.solvers.compartmental.Compartmental
    orig = Comp.__dict__['round_preserving']
    def drawn(y):
        y = np.asarray(y, dtype=float)
        n = int(round(y.sum()))
        return rng.multinomial(n, y/y.sum()).astype(float) if n > 0 else np.zeros_like(y)
    Comp.round_preserving = staticmethod(drawn)
    try:
        yield
    finally:
        Comp.round_preserving = orig


def use_drawn_start(on=True, salt=30):
    """
    Draw agents' starting states in each group instead of rounding them (for every agent sim initialized afterwards in this process).

    Starsim's agent solver gives each state given by group (`init=` with `ss.by`) an exact count in each age and sex group,
    round(share × agents in the group), keeping the group's total. A rare state's expected count is below 1 in a small group, so it rounds to
    nobody, in every seed: with 100,000 agents by single year of age and sex, 65+ started with 0.49–0.73 of the ODE's active TB (exp 30).
    Drawn counts are right on average; in distribution, the same as drawing each agent's state from its group's shares.

    Args:
        on (bool): True to draw (idempotent), False to restore Starsim's rounding
        salt (int): with the sim's rand_seed, seeds the draws (30 reproduces exp 30's runs)
    """
    Agents = ss.solvers.agents.Agents
    current = Agents.init_post
    base = getattr(current, '_tbode4_orig', current)
    if not on:
        Agents.init_post = base
        return
    def init_post(self):
        rng = np.random.default_rng([int(self.sim.pars.rand_seed), salt])
        with drawn_counts(rng):
            return base(self)
    init_post._tbode4_orig = base
    Agents.init_post = init_post
    return


#%% Starting shares by age and sex, from the ODE

def _by_age_sex(res, sex):
    """ A stratified result for one sex by age group: values [time, group] and the groups' lower edges """
    b = res.select(sex=sex).by('age')
    return np.asarray(b.values, dtype=float), np.array(ss.data.age_edges(b.columns), dtype=float)


def shares_by_group(sim, idx=-1, states=None):
    """
    Each TB state's share of each age and sex group of a sim stratified by age and sex, at time index idx.

    Returns:
        shares ({state: array [sex, age group]}), the groups' lower edges, and the people in each group [sex, age group]
    """
    r = sim.results.tb
    if states is None: # tbsim.TB's states and TBTreat's treatment states (generated with exp 35's options)
        info = getattr(sim.diseases.tb, 'tx_info', None)
        states = tb_states + ([k for k in info.states if k not in tb_states] if info is not None else [])
    vals = {k: np.stack([_by_age_sex(r[k], s)[0][idx] for s in sexes]) for k in states if k in r}
    pop = sum(vals.values())
    lows = _by_age_sex(r['susceptible'], sexes[0])[1]
    return {k: v/pop for k, v in vals.items()}, lows, pop


def coarsen(shares, lows, pop, edges):
    """ Shares in groups with lower edges `lows`, averaged over each coarser group (lower edges `edges`), weighted by its people `pop` """
    j = np.clip(np.digitize(lows + 1e-9, edges) - 1, 0, len(edges) - 1)
    out = {}
    for k, v in shares.items():
        num = np.zeros((2, len(edges)))
        den = np.zeros((2, len(edges)))
        for sx in range(2):
            np.add.at(num[sx], j, v[sx]*pop[sx])
            np.add.at(den[sx], j, pop[sx])
        out[k] = np.divide(num, den, out=np.zeros_like(num), where=den > 0)
    return out


def init_table(shares, lows):
    """
    `init=` for the TB module: each state but susceptible as an `ss.by(age, sex)` table; susceptible takes the rest (agents refuse a share by group
    for every state). States with no one in them are left out.
    """
    out = {}
    for k, v in shares.items():
        if k == 'susceptible' or not np.any(v > 0):
            continue
        out[k] = ss.by(pd.DataFrame(dict(age=np.tile(lows, 2), sex=np.repeat(sexes, len(lows)), value=np.concatenate([v[0], v[1]]))))
    return out


def settled_start(sim, edges=None, drift_years=20, states=None):
    """
    A starting state by age and sex from the end of a run of the ODE, for both solvers.

    Run the ODE, stratified by age and sex in groups at least as fine as `edges`, with the starting year's demography frozen and no intervention,
    until its shares by age settle (200 years in exp 30). Its shares at the end are averaged to the groups `edges` (weighted by its people).
    Give agents age groups no coarser than `edges` (e.g. single years), or they read the table per agent and start with wrong shares when several
    states are given; and call `use_drawn_start()`. The ODE takes the table whatever its own groups, as long as they're as fine as `edges`.

    Args:
        sim (ss.Sim): the ODE, run
        edges (list): lower edges of the table's age groups (default: single years to 100+)
        drift_years (float): the drift is the largest change in any group's share over this many years before the end

    Returns:
        objdict: lows, shares, people (the run's groups at the end), drift ({state: largest change}), edges, table ({state: [sex, group]}),
        and init (the `init=` tables)
    """
    edges = np.asarray(age_edges(1.0) if edges is None else edges, dtype=float)
    sh, lows, pop = shares_by_group(sim, states=states)
    yv = np.asarray(sim.t.yearvec, dtype=float)
    j = int(np.argmin(np.abs(yv - (yv[-1] - drift_years))))
    before, _, _ = shares_by_group(sim, j, states=states)
    drift = {k: float(np.max(np.abs(sh[k] - before[k]))) for k in sh}
    table = coarsen(sh, lows, pop, edges)
    return sc.objdict(lows=lows, shares=sh, people=pop, drift=drift, edges=edges, table=table, init=init_table(table, edges))
