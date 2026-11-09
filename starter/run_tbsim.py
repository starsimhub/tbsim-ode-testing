"""
A starting point for running TBsim v4 as agents, with the conventions worked out in this repo. See README.md.

Starts 50,000 agents at the model's endemic equilibrium in 1950, runs to 2000 at monthly steps, saves the results, and plots
prevalence, incidence and TB mortality over time. Run from the repo root:

    v4/.venv/bin/python starter/run_tbsim.py
"""
import sys
import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss
import tbsim
import matplotlib.pyplot as plt

outdir = sc.thispath(__file__)
sys.path.insert(0, str(outdir.parent / 'v4'))
from tbode4 import fastode # noqa: E402  # A standalone copy of TBsim's ODE, used to find the starting equilibrium

# Settings
n_agents = 50_000
n_reps = 4 # Independent runs, each with its own random seed
seed0 = 1
start, stop = 1950, 2000 # Runs take about 20-25 years to spread out from their shared starting state (README)
dt = ss.months(1)
ncpus = 4 # Runs in parallel, one per core

# Parameters: one draw from the posterior calibrated to India's prevalence survey and progression data (rates per year)
pi = 0.3518 # Risk of reinfection for anyone infected before, relative to someone never infected
phi = 0.01194 # Progression from remote latent infection, relative to recent
inf_non, inf_asy = 0.3153, 0.3735 # Progression from recent latent infection to non-infectious and asymptomatic TB
pars = dict(
    beta = 22.10, # Transmission rate per symptomatic case (homogeneous mixing)
    trans_asymp = 0.7488, # Infectiousness of asymptomatic TB, relative to symptomatic
    rr_reinfection_cleared = pi, rr_reinfection_rec = pi, rr_reinfection_treat = pi, rr_reinfection_inf = pi,
    inf_cle = 0.6638, # Clearance of latent infection
    inf_non = inf_non, inf_asy = inf_asy,
    k_latent = 7.517, # Recent latent -> remote latent (so most progression happens in the first months)
    inf_non_remote = phi*inf_non, inf_asy_remote = phi*inf_asy,
    non_rec = 0.2207, # Non-infectious TB -> recovered
    non_asy = 0.2776, # Non-infectious -> asymptomatic TB
    asy_non = 1.612, # Asymptomatic -> non-infectious TB
    asy_sym = 0.7560, # Asymptomatic -> symptomatic TB
    sym_asy = 0.7593, # Symptomatic -> asymptomatic TB
    sym_dead = 0.3959, # TB death, from symptomatic TB
    mu = 1/70, # Background death rate; births replace every death, so the population stays constant
)


def equilibrium(pars):
    """ The endemic equilibrium of the ODE version of the model: the share of people in each state """
    p = fastode.make_pars(**pars)
    x = fastode.equilibrium(p)
    init = dict(zip(fastode.state_names, x))
    non0, non1 = init.pop('non0'), init.pop('non1') # The ODE splits non-infectious TB by whether the episode has been asymptomatic yet
    init['non_infectious'] = non0 + non1
    return init, non1/(non0 + non1), fastode.summary(x, p)


class Incidence(ss.Analyzer):
    """ Incidence as the first entry to asymptomatic TB in each episode, counted from the transitions in each step """
    def __init__(self, share_flagged=0.0, **kwargs):
        super().__init__(**kwargs)
        self.share_flagged = share_flagged # Share of people starting in non-infectious TB whose episode has already been asymptomatic
        return

    def init_post(self):
        super().init_post()
        tb = self.sim.diseases.tb
        self.flag = np.zeros(2*len(self.sim.people), bool) # Whether each agent has been asymptomatic in their current episode
        for state in ['asymptomatic', 'symptomatic']:
            self._set(np.asarray(getattr(tb, state).uids), True)
        non = np.asarray(tb.non_infectious.uids)
        rng = np.random.default_rng(int(self.sim.pars.rand_seed) + 7)
        self._set(non[rng.random(len(non)) < self.share_flagged], True)
        self.first = np.zeros(self.sim.t.npts)
        return

    def _grow(self, uids):
        """ Make room in the flags for agents born during the run """
        if len(uids) and uids.max() >= len(self.flag):
            self.flag = np.concatenate([self.flag, np.zeros(max(len(self.flag), uids.max() + 1 - len(self.flag)), bool)])
        return

    def _set(self, uids, value):
        self._grow(uids)
        self.flag[uids] = value
        return

    def step(self):
        tb = self.sim.diseases.tb
        get = lambda name: np.asarray(getattr(tb, name).uids) # The agents who made this transition in this step
        from_latent = np.concatenate([get('latent_asy'), get('remote_asy')]) # Latent infection -> asymptomatic TB: always a first entry
        worsening = get('worsening') # Non-infectious -> asymptomatic TB: a first entry only if not yet asymptomatic this episode
        self._grow(worsening)
        self.first[self.sim.ti] = len(from_latent) + (~self.flag[worsening]).sum()
        self._set(np.concatenate([from_latent, worsening, get('remission')]), True)
        self._set(get('recovery'), False) # Recovery ends the episode
        return


def make_sim(seed, init, share_flagged):
    tb = tbsim.TB(init=init, **pars) # Agents start in each state in the equilibrium's proportions
    return ss.Sim(tb, solver='abm', n_agents=n_agents, total_pop=n_agents, start=start, stop=stop, dt=dt, rand_seed=seed, verbose=0,
                  analyzers=Incidence(share_flagged=share_flagged)) # No network: everyone mixes homogeneously


def run_one(seed, init, share_flagged):
    """ One run, as yearly averages per 100,000 people """
    sim = make_sim(seed, init, share_flagged)
    sim.run()
    r = sim.results.tb
    alive = sim.results.alive.values
    per100k = 1e5/alive
    df = pd.DataFrame(dict(
        year = np.floor(np.array(sim.t.yearvec)).astype(int),
        prevalence = (r['asymptomatic'].values + r['symptomatic'].values)*per100k, # Active TB (asymptomatic + symptomatic), per 100,000
        incidence = sim.analyzers[0].first*per100k/sim.t.dt_year, # New episodes per 100,000 per year
        mortality = r['tb_death'].values*per100k/sim.t.dt_year, # TB deaths per 100,000 per year
    ))
    df = df[df.year < stop].groupby('year').mean().reset_index() # Average the monthly values over each full year
    df['seed'] = seed
    return df


def plot(df, ref):
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    labels = dict(prevalence='Active TB prevalence\n(per 100,000)', incidence='Incidence\n(per 100,000 per year)', mortality='TB mortality\n(per 100,000 per year)')
    for ax, (key, label) in zip(axes, labels.items()):
        for seed, d in df.groupby('seed'):
            ax.plot(d.year, d[key], color='#7fa7d6', lw=0.8, alpha=0.7, label='Each run' if seed == seed0 else None)
        mean = df.groupby('year')[key].mean()
        ax.plot(mean.index, mean.values, color='#1f4e8c', lw=2, label=f'Mean of {n_reps} runs')
        ax.axhline(ref[key], color='k', ls='--', lw=1, label='ODE equilibrium')
        ax.set_title(label, fontsize=10)
        ax.set_xlabel('Year')
        ax.set_ylim(bottom=0)
    axes[0].legend(fontsize=8, frameon=False)
    fig.suptitle(f'TBsim v4, {n_agents:,} agents, homogeneous mixing, no treatment', fontsize=10)
    fig.tight_layout()
    fig.savefig(outdir / 'tbsim_burn_in.png', dpi=150)
    return fig


if __name__ == '__main__':
    init, share_flagged, ode = equilibrium(pars)
    ref = dict(prevalence=ode.prevalence*1e5, incidence=ode.incidence*1e5, mortality=ode.mortality*1e5)
    print(f'ODE equilibrium per 100,000: prevalence {ref["prevalence"]:.0f}, incidence {ref["incidence"]:.0f}, TB mortality {ref["mortality"]:.0f}')

    T = sc.timer()
    seeds = [seed0 + i for i in range(n_reps)]
    dfs = sc.parallelize(run_one, iterarg=seeds, kwargs=dict(init=init, share_flagged=share_flagged), ncpus=ncpus)
    df = pd.concat(dfs, ignore_index=True)
    print(f'{n_reps} runs in {T.toc(output=True):.0f} s')

    df.to_csv(outdir / 'tbsim_burn_in.csv', index=False) # Yearly outputs of every run
    sc.save(outdir / 'tbsim_burn_in.obj', sc.objdict(pars=pars, init=init, ode_reference=ref, results=df)) # Everything needed to redraw the figure
    plot(df, ref)
    print(f'Mean of the runs over {start}-{stop - 1}, per 100,000:')
    print(df[['prevalence', 'incidence', 'mortality']].mean().round(0).to_string())
    print('Saved tbsim_burn_in.csv, tbsim_burn_in.obj and tbsim_burn_in.png')
