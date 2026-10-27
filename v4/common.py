"""
Shared code for the v4 redo of the key ODE↔ABM findings: parameter sets, the TB models, the ODE equilibrium and
calibration, sim builders, and the outputs compared.

Under Starsim v4, `tbsim.TB` is one declared model that runs as an ODE (`solver='ode'`) or as agents
(`solver='abm'`), so there is no separate ODE to keep in step with the ABM. `TBClasses` is the same model with
progression that declines with time since infection, represented exactly by Poisson-mixture latent classes
(as `tbode.TBODE` and `tbode.TBExact` did under Starsim 3).
"""
import json
import numpy as np
import sciris as sc
import starsim as ss
import tbsim
from scipy.optimize import brentq
from scipy.stats import poisson

__all__ = ['folder', 'ncpus', 'parsets', 'k_tsi', 'TBClasses', 'latent_classes', 'make_tb', 'make_sim', 'equilibrium', 'calibrate',
           'get_calib', 'outputs', 'progression', 'versions']

folder = sc.thispath(__file__)
ncpus = 4 # Runs are on Cliff's laptop: at most 4 cores
k_tsi = 5.845726750296858 # tb_LAI_TPT's best-fit k_asy = k_non (per year)

# Natural-history parameter sets, as in exps 01–11 (rates per year). A: tbsim's defaults. B: tb_LAI_TPT's best fit (row 0 of
# collapsed_samples_eps1.2.csv, pre-1995 sym_dead), with reinfection risk ratio pi for everyone no longer infected
parsets = sc.objdict(
    A = dict(),
    B = dict(
        inf_cle     = 0.5365965386150265,
        inf_non     = 0.4103330692741019,
        inf_asy     = 0.2741560588304815,
        non_rec     = 0.1706711906227213,
        non_asy     = 0.2426759489829796,
        asy_non     = 1.738594655537684,
        asy_sym     = 0.7435451323482238,
        sym_asy     = 0.6630323475395543,
        trans_asymp = 0.9407760328175468,
        sym_dead    = 0.2626559974664169,
        rr_reinfection_cleared = 0.2415635040650483,
        rr_reinfection_rec     = 0.2415635040650483,
        rr_reinfection_treat   = 0.2415635040650483,
    ),
)

# The states compared, as shares of everyone alive; 'cleared' sums the three states of people no longer infected
groups = sc.objdict(
    infection      = ['infection', 'infection_remote'], # Latent infection (TBClasses sums its classes into its result 'infection')
    non_infectious = ['non_infectious'],
    asymptomatic   = ['asymptomatic'],
    symptomatic    = ['symptomatic'],
    cleared        = ['cleared', 'recovered', 'treated'],
    prevalence     = ['asymptomatic', 'symptomatic'],
)


def latent_classes(inf_non, inf_asy, k, tol=1e-10):
    """
    The Poisson-mixture latent classes that give progression hazards `inf_non·exp(−k·τ)` and `inf_asy·exp(−k·τ)` exactly

    A new infection joins class m with probability Poisson(m; r0/k), r0 = inf_non + inf_asy, and progresses at the constant
    rate m·k (class 0 never progresses), to non-infectious or asymptomatic TB in the ratio inf_non : inf_asy. The weights are
    truncated where the remaining weight is below `tol` (added to the last class). As `tbode.latent_classes`.

    Returns:
        weights (array): the share of new infections in each class
        rates (array): each class's rate of progression, per year
    """
    lam = (inf_non + inf_asy)/k
    M = int(poisson.isf(tol, lam)) + 1
    m = np.arange(M + 1)
    weights = poisson.pmf(m, lam)
    weights[-1] += 1 - weights.sum()
    return weights, m*k


class TBClasses(ss.Disease):
    """
    `tbsim.TB`'s natural history (no treatment) with progression that falls with time since infection, as tbsim's `k_asy = k_non = k`

    Latent infection is split into the Poisson-mixture classes of `latent_classes()`: `infection_0` (never progresses),
    `infection_1`, ... The result `infection` is their sum. Reinfection of latent people is off, as in exps 08–11.

    Args:
        k (float): the rate per year at which the hazard of progression falls
        kwargs (dict): parameters, as for `tbsim.TB` (rates per year)
    """
    def __init__(self, k=k_tsi, init=None, name='tb', **kwargs):
        p = ss.pars(tbsim.tb_pars(), **kwargs)
        weights, rates = latent_classes(p.inf_non.raw.value, p.inf_asy.raw.value, k)
        frac_non = p.inf_non.raw.value/(p.inf_non.raw.value + p.inf_asy.raw.value)
        super().__init__(init=init, name=name)
        classes = [ss.State(f'infection_{m}', label=f'Latent infection, class {m}') for m in range(len(weights))]
        self.define_states(
            ss.State('susceptible', label='Susceptible'),
            *classes,
            ss.State('cleared',        label='Cleared (from latent)'),
            ss.State('recovered',      label='Recovered (from non-infectious)'),
            ss.State('non_infectious', label='Non-infectious TB'),
            ss.State('asymptomatic',   label='Asymptomatic TB'),
            ss.State('symptomatic',    label='Symptomatic TB'),
        )
        classes = [getattr(self, f'infection_{m}') for m in range(len(weights))]
        foi = p.beta*(p.trans_asymp*self.asymptomatic + self.symptomatic)
        split = [None] + list(weights[1:]) # Class 0 takes the rest
        self.define_transitions(
            transmission = ss.Transition(self.susceptible, classes, rate=foi, p=split),
            reinfect_cle = ss.Transition(self.cleared,     classes, rate=p.rr_reinfection_cleared*foi, p=split),
            reinfect_rec = ss.Transition(self.recovered,   classes, rate=p.rr_reinfection_rec*foi, p=split),
        )
        for m, state in enumerate(classes):
            self.define_transitions(ss.Transition(state, self.cleared, rate=p.inf_cle, name=f'clearance_{m}'))
            if m > 0:
                self.define_transitions(
                    ss.Transition(state, self.non_infectious, rate=ss.peryear(rates[m]*frac_non),       name=f'latent_non_{m}'),
                    ss.Transition(state, self.asymptomatic,   rate=ss.peryear(rates[m]*(1 - frac_non)), name=f'latent_asy_{m}'),
                )
        self.define_transitions(
            recovery    = self.non_infectious >> self.recovered      @ p.non_rec,
            worsening   = self.non_infectious >> self.asymptomatic   @ p.non_asy,
            improvement = self.asymptomatic   >> self.non_infectious @ p.asy_non,
            symptoms    = self.asymptomatic   >> self.symptomatic    @ p.asy_sym,
            remission   = self.symptomatic    >> self.asymptomatic   @ p.sym_asy,
            tb_death    = ss.death(self.symptomatic, rate=p.sym_dead),
            births      = ss.birth(self.susceptible, rate=p.mu),
            replace     = ss.null >> self.susceptible @ p.sym_dead*self.symptomatic.n, # Births that replace TB deaths, as tbsim.TB
            deaths      = ss.death(rate=p.mu),
        )
        infection = classes[0]
        for state in classes[1:]:
            infection = infection + state
        self.define_results(infection=infection)
        self.n_classes = len(weights)
        return


def make_tb(parset='A', beta=None, init=None, model='tb', **kwargs):
    """ The TB module: `model='tb'` for `tbsim.TB`, `'classes'` for `TBClasses` (k = 5.846) """
    pars = sc.mergedicts(parsets[parset], kwargs)
    if beta is not None:
        pars['beta'] = beta
    if model == 'tb':
        return tbsim.TB(init=init, **pars)
    elif model == 'classes':
        return TBClasses(init=init, **pars)
    raise ValueError(f'Unknown model "{model}"')


def make_sim(tb, solver='abm', dt=ss.months(1), start=1900, stop=2000, n_agents=50_000, seed=0, network=None):
    """
    A sim of one TB module. As agents, people mix homogeneously (v4's default with no network: the ODE's mass action), unless
    `network` is the mean number of contacts K of an `ss.RandomNet` with new contacts every step (the module's `beta` is then
    per contact: the ODE's `beta`/K)
    """
    kw = dict(start=start, stop=stop, dt=dt, total_pop=n_agents, rand_seed=seed, verbose=0)
    if solver == 'ode':
        return ss.Sim(tb, solver='ode', **kw)
    if network:
        kw['networks'] = ss.RandomNet(n_contacts=ss.constant(network), dur=0)
    return ss.Sim(tb, solver=solver, n_agents=n_agents, **kw)


def end_shares(sim):
    """ The share of people in each state at the end of a run """
    r = sim.results.tb
    alive = sim.results.alive.values[-1]
    return {s.name: float(r[s.name].values[-1]/alive) for s in sim.diseases.tb.graph.people_states}


def equilibrium(parset='A', beta=None, model='tb', years=2000, tol=1e-9, maxiter=50, strict=True, **kwargs):
    """
    The ODE's endemic equilibrium: the share of people in each state, found by running the ODE (yearly steps) until the
    shares change by less than `tol` per century (with `strict=False`, whatever the shares are after `maxiter` tries)
    """
    init = dict(symptomatic=0.01)
    for i in range(maxiter):
        sim = make_sim(make_tb(parset, beta, init=init, model=model, **kwargs), solver='ode', dt=ss.years(1), start=0, stop=years, n_agents=1)
        sim.run()
        shares = end_shares(sim)
        sim2 = make_sim(make_tb(parset, beta, init=shares, model=model, **kwargs), solver='ode', dt=ss.years(1), start=0, stop=100, n_agents=1)
        sim2.run()
        new = end_shares(sim2)
        change = max(abs(new[k] - shares[k]) for k in shares)
        if change < tol:
            return new
        init = new
    if not strict:
        return new
    raise RuntimeError(f'ODE did not reach equilibrium (last change {change:.2e})')


def prevalence(shares):
    return shares['asymptomatic'] + shares['symptomatic']


def calibrate(parset='A', model='tb', target=0.005, bounds=(0.5, 50), **kwargs):
    """ The ODE `beta` (per year) whose equilibrium active (asymptomatic + symptomatic) prevalence is `target` """
    gap = lambda beta: prevalence(equilibrium(parset, beta, model=model, years=3000, tol=1e-7, maxiter=3, strict=False, **kwargs)) - target
    beta = brentq(gap, *bounds, xtol=1e-6)
    eq = equilibrium(parset, beta, model=model, **kwargs)
    return sc.objdict(beta=beta, equilibrium=eq, prevalence=prevalence(eq))


def get_calib(key, parset, model='tb', **kwargs):
    """ The calibrated `beta` and equilibrium for an arm, computed once and kept in `results/calibration.json` """
    fn = folder / 'results' / 'calibration.json'
    calib = sc.loadjson(fn) if fn.exists() else {}
    if key not in calib:
        out = calibrate(parset, model=model, **kwargs)
        calib[key] = dict(parset=parset, model=model, kwargs=kwargs, beta=out.beta, equilibrium=out.equilibrium)
        sc.savejson(fn, calib, indent=2)
    return sc.objdict(calib[key])


def outputs(sim):
    """
    The outputs compared, as an `sc.objdict` of arrays: the share of people in each group of states, and per person per
    year, incidence (every entry to asymptomatic TB from latent infection or non-infectious TB, counted from the
    transitions) and TB deaths
    """
    r = sim.results.tb
    alive = sim.results.alive.values
    out = sc.objdict(year=np.array(sim.t.yearvec))
    tb = sim.diseases.tb
    names = set(r.keys())
    for key, states in groups.items():
        out[key] = sum(r[s].values for s in states if s in names)/alive
    dt = sim.t.dt_year
    into_asy = [k for k in names if k == 'latent_asy' or k.startswith('latent_asy_') or k == 'worsening']
    out.incidence = sum(r[k].values for k in into_asy)/alive/dt
    out.tb_deaths = r.tb_death.values/alive/dt
    return out


def progression(parset='B', model='classes', years=10, **kwargs):
    """
    ProgressionY1 and ProgressionRatio (tb_LAI_TPT's targets 0.0523 / 0.663) from an ODE cohort, as tb_LAI_TPT's
    cohort_progression and `tbode.ode_cohort_progression`: everyone latent at t = 0, no transmission, no TB or background
    deaths. Y1 is the share that has ever reached asymptomatic TB after a year (asymptomatic is made absorbing, so its
    share is the share ever there); the ratio is Y1 over the same share after `years` (10)
    """
    pars = sc.mergedicts(parsets[parset], kwargs)
    if model == 'classes':
        w, _ = latent_classes(pars.get('inf_non', 0.16), pars.get('inf_asy', 0.06), k_tsi)
        init = {f'infection_{m}': float(w[m]) for m in range(len(w))}
    else:
        init = dict(infection=1.0)
    tb = make_tb(parset, beta=0, init=init, model=model, sym_dead=0, mu=0, asy_non=0, asy_sym=0, **kwargs)
    sim = make_sim(tb, solver='ode', dt=ss.days(1), start=0, stop=years, n_agents=1).run()
    ever = sim.results.tb.asymptomatic.values/sim.results.alive.values
    yrs = np.array(sim.t.yearvec)
    y1 = np.interp(1, yrs, ever)
    y10 = np.interp(years, yrs, ever)
    return sc.objdict(ProgressionY1=y1, ProgressionRatio=y1/y10)


def versions():
    """ The commits of tbsim and starsim (ssv4) used """
    out = {}
    for name, path in [('tbsim', '/home/cliffk/idm/tbsim'), ('starsim', '/home/cliffk/idm/ssv4')]:
        info = sc.gitinfo(path)
        out[name] = f"{info['branch']} @ {info['hash']}" + (' (+ uncommitted changes)' if sc.runcommand(f'git -C {path} status --porcelain -uno').strip() else '')
    return out
