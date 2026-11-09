"""
Treatment as declared states and transitions (exp 23; the design is CONVENTIONS.md 9), for tbsim v4's ODE and agents alike.

- `TBTreat`: tbsim.TB with exp 23's treatment, in place of tbsim's own (which stays inert: theta = 0, p_relapse = 0)
- `lai_tpt_treatment()`: tb_LAI_TPT's best-fit coverage and efficacy over time, as ss.timeseries
- `Incidence`: agents' incidence, the first entry to asymptomatic TB per episode plus relapses (CONVENTIONS.md 5)
- `TxChecks`: agents' checks of the treatment pitfalls (exp 23's README)
- `ode_incidence()`: the same incidence from a run of v4's ODE, which has one non-infectious state
"""
import numpy as np
import sciris as sc
import starsim as ss
import tbsim

__all__ = ['TBTreat', 'theta_from_coverage', 'lai_tpt_treatment', 'Incidence', 'TxChecks', 'ode_incidence', 'tx_outputs']

tx_state_names = ['treatment_success', 'treatment_fail', 'relapsing']


class TBTreat(tbsim.TB):
    """
    tbsim.TB with treatment as in exp 23: people with symptomatic TB start treatment at `tx_rate`; a share `tx_success` start a course
    that will cure them, the rest one that will fail. Courses last `tx_dur`: exactly (course='fixed') or exponentially (course='rate').
    Failing courses face TB death at sym_dead, then return to symptomatic TB. A share `tx_p_relapse` of cures relapse to symptomatic TB
    at rate 1/dur_relapse_tx; the rest go to `treated`. Nobody on treatment is infectious or reinfected.

    Args:
        tx_rate (float/ss.timeseries): the rate of starting treatment per year; a timeseries of ss.peryear values to vary in time
        tx_success (float/ss.timeseries): adherence × efficacy, the chance a course will cure, set when it starts
        tx_dur (float): the course length (years)
        course (str): 'fixed' (exact, a conveyor of whole steps in the ODE) or 'rate'
        stages (int): with course='fixed', run the ODE's course as this many Erlang stages instead (agents keep the exact course)
        tx_p_relapse (float): the share of cures that relapse
        dur_relapse_tx (float): the mean time from cure to relapse (years)
        kwargs (dict): passed to tbsim.TB
    """
    def __init__(self, tx_rate=0.0, tx_success=1.0, tx_dur=0.5, course='fixed', stages=None, tx_p_relapse=0.0, dur_relapse_tx=1.5, **kwargs):
        for key in ['theta', 'phi', 'p_relapse']:
            if key in kwargs:
                errormsg = f'TBTreat replaces tbsim.TB\'s own treatment, so {key}= is not used; give tx_rate, tx_success, tx_p_relapse instead'
                raise ValueError(errormsg)
        if course not in ['fixed', 'rate']:
            raise ValueError(f'course is "fixed" or "rate", not "{course}"')
        super().__init__(**kwargs)
        sym_dead = self.tb_death.rate # The sym_dead parameter itself, so a change to it applies to deaths on treatment too
        rate = tx_rate if isinstance(tx_rate, ss.timeseries) else ss.peryear(tx_rate)
        dur = ss.constant(ss.years(tx_dur)) if course == 'fixed' else ss.years(tx_dur) # A plain duration is the mean of an exponential
        staged = dict(stages=stages) if (stages and course == 'fixed') else {}
        self.define_states(
            ss.State('treatment_success', label='On treatment, will be cured'),
            ss.State('treatment_fail',    label='On treatment, will fail'),
            ss.State('relapsing',         label='Cured, will relapse'),
        )
        self.define_transitions(
            tx_start    = ss.Transition(self.symptomatic, [self.treatment_success, self.treatment_fail], rate=rate, p=[tx_success, None]),
            tx_cure     = ss.Transition(self.treatment_success, [self.treated, self.relapsing], dur=dur, p=[None, tx_p_relapse], **staged),
            tx_fail     = ss.Transition(self.treatment_fail, self.symptomatic, dur=dur, **staged),
            tb_death_tx = ss.death(self.treatment_fail, rate=sym_dead),
            relapse     = ss.Transition(self.relapsing, self.symptomatic, dur=ss.years(dur_relapse_tx)),
        )
        if isinstance(getattr(self, 'replace', None), ss.Transition): # Births replace TB deaths on treatment too, as tbsim.TB's replace does those from symptomatic TB
            self.define_transitions(replace_tx=ss.null >> self.susceptible @ sym_dead*self.treatment_fail.n)
        self.tx = sc.objdict(tx_rate=tx_rate, tx_success=tx_success, tx_dur=tx_dur, course=course, stages=stages, tx_p_relapse=tx_p_relapse, dur_relapse_tx=dur_relapse_tx)
        return


def theta_from_coverage(c):
    """ The rate per year at which a share c of people would start treatment within a year: −ln(1 − c) """
    return -np.log1p(-np.asarray(c, dtype=float))


def lai_tpt_treatment(row=None, step=1/12, end=2024):
    """
    tb_LAI_TPT's best-fit treatment over time (row 0 of collapsed_samples_eps1.2.csv), as tb_LAI_TPT builds it: coverage 0 until
    the end of 1994, linear to tx_covg_2010 in 2010 and tx_covg_2024 in 2024, then constant; efficacy 0 until the end of 1994,
    efficacy_1995 from 1995, linear to efficacy_2024 in 2024, then constant; adherence 1.

    Returns tx_rate (ss.timeseries of ss.peryear, from coverage at every `step` years) and tx_success (ss.timeseries), plus the knots
    """
    r = sc.mergedicts(dict(tx_covg_2010=0.5478919423969836, tx_covg_2024=0.6772271561810697, efficacy_1995=0.4525494678080442,
                           efficacy_2024=0.5300106194827879, p_relapse=0.03394965858858717), row)
    t0 = 1994 + 364/365 # tb_LAI_TPT's first knot, 1994-12-31
    cov_t, cov_v = [t0, 2010, 2024], [0.0, r['tx_covg_2010'], r['tx_covg_2024']]
    times = np.round(np.arange(1995, end + step/2, step), 6)
    times = np.concatenate([[t0], times])
    covg = np.interp(times, cov_t, cov_v)
    tx_rate = ss.timeseries({float(t): ss.peryear(float(theta_from_coverage(c))) for t, c in zip(times, covg)}, mode='linear')
    tx_success = ss.timeseries({t0: 0.0, 1995.0: r['efficacy_1995'], 2024.0: r['efficacy_2024']}, mode='linear')
    coverage = lambda t: np.interp(t, cov_t, cov_v, left=0)
    efficacy = lambda t: np.interp(t, [t0, 1995, 2024], [0.0, r['efficacy_1995'], r['efficacy_2024']], left=0)
    return sc.objdict(tx_rate=tx_rate, tx_success=tx_success, coverage=coverage, efficacy=efficacy, p_relapse=r['p_relapse'])


class _Grow:
    """ Per-agent arrays that grow with the population """
    def _arr(self, name, dtype, fill):
        setattr(self, name, np.full(2*len(self.sim.people), fill, dtype=dtype))

    def _grow(self, uids):
        n = (int(np.max(uids)) + 1) if len(uids) else 0
        for name in self._arrays:
            a = getattr(self, name)
            if n > len(a):
                extra = np.full(max(len(a), n - len(a)), a.dtype.type(0) if a.dtype != float else np.nan, dtype=a.dtype)
                setattr(self, name, np.concatenate([a, extra]))
        return


def _uids(module, name):
    """ The agents who made a transition in this step (none if the model has no such transition) """
    tr = getattr(module, name, None)
    return np.asarray(tr.uids) if tr is not None else np.array([], dtype=int)


class Incidence(ss.Analyzer, _Grow):
    """
    Agents' incidence: the first entry to asymptomatic TB per episode, plus relapses (CONVENTIONS.md 5), counted from transitions

    Each agent has a flag for "has been asymptomatic in this episode". It is set on entering asymptomatic TB (from latent infection,
    non-infectious TB or symptomatic TB) and on relapse; it is cleared when a new episode can start: entering non-infectious TB from
    latent infection, recovery, and cure. A failed course keeps it, since the episode goes on.

    Args:
        share_flagged (float): the share of people starting in non-infectious TB who have been asymptomatic in this episode
    """
    def __init__(self, share_flagged=0.0, **kwargs):
        super().__init__(**kwargs)
        self.share_flagged = share_flagged
        return

    def init_post(self):
        super().init_post()
        tb = self.sim.diseases.tb
        self._arrays = ['flag']
        self._arr('flag', bool, False)
        flagged = ['asymptomatic', 'symptomatic', 'treatment_success', 'treatment_fail'] # Episodes under way at the start
        for state in flagged:
            if hasattr(tb, state):
                u = np.asarray(getattr(tb, state).uids)
                self._grow(u)
                self.flag[u] = True
        non = np.asarray(tb.non_infectious.uids)
        rng = np.random.default_rng(int(self.sim.pars.rand_seed) + 7)
        self.flag[non[rng.random(len(non)) < self.share_flagged]] = True
        npts = self.sim.t.npts
        self.first = np.zeros(npts) # New episodes in each step, relapses included
        self.relapses = np.zeros(npts)
        return

    def step(self):
        tb = self.sim.diseases.tb
        get = lambda name: _uids(tb, name)
        new_latent_tb = np.concatenate([get('latent_non'), get('remote_non')]) # Non-infectious TB from latent infection: not yet asymptomatic
        self._grow(new_latent_tb)
        self.flag[new_latent_tb] = False
        from_latent = np.concatenate([get('latent_asy'), get('remote_asy')]) # Always a first entry
        worsening, relapse = get('worsening'), get('relapse')
        self._grow(np.concatenate([from_latent, worsening, relapse, [0]]).astype(int))
        self.first[self.ti] = len(from_latent) + (~self.flag[worsening]).sum() + len(relapse)
        self.relapses[self.ti] = len(relapse)
        on = np.concatenate([from_latent, worsening, get('remission'), relapse]).astype(int)
        off = np.concatenate([get('recovery'), get('tx_cure')]).astype(int) # The episode ends
        self._grow(np.concatenate([on, off, [0]]).astype(int))
        self.flag[on] = True
        self.flag[off] = False
        return


class TxChecks(ss.Analyzer, _Grow):
    """
    Agents' checks of the treatment pitfalls in exp 23's README, from each step's transitions and the states at the end of the last step.
    Counts are summed over the run; `by_year` holds the courses started each year and the share that will cure.
    """
    def init_post(self):
        super().init_post()
        self._arrays = ['ti_start', 'prev_tx', 'prev_succ', 'prev_sym', 'prev_rel']
        self._arr('ti_start', float, np.nan)
        for name in self._arrays[1:]:
            self._arr(name, bool, False)
        self.c = sc.objdict({k: 0 for k in ['starts', 'fails', 'cures', 'restart_same_step', 'start_while_on_tx', 'start_not_symptomatic',
                                            'relapse_not_from_success', 'tb_death_from_success', 'bg_death_from_success', 'success_person_steps',
                                            'symptomatic_person_steps', 'courses_timed', 'course_steps_sum', 'course_steps_sq']})
        self.course_steps = sc.ddict(int) # Course lengths in steps, of courses started during the run
        self.by_year = sc.ddict(lambda: np.zeros(2)) # Courses started each year, and those that will cure
        self.alive = np.zeros(self.sim.t.npts)
        self._remember()
        return

    def _remember(self):
        tb = self.sim.diseases.tb
        for name, states in [('prev_tx', ['treatment_success', 'treatment_fail']), ('prev_succ', ['treatment_success']),
                             ('prev_sym', ['symptomatic']), ('prev_rel', ['relapsing'])]:
            arr = getattr(self, name)
            arr[:] = False
            for state in states:
                u = np.asarray(getattr(tb, state).uids)
                self._grow(u)
                getattr(self, name)[u] = True
        return

    def step(self):
        tb, c = self.sim.diseases.tb, self.c
        get = lambda name: _uids(tb, name)
        starts, fails, cures = get('tx_start'), get('tx_fail'), get('tx_cure')
        allu = np.concatenate([starts, fails, cures, get('tb_death'), get('deaths'), get('symptoms'), get('relapse'), [0]]).astype(int)
        self._grow(allu)
        c.starts += len(starts); c.fails += len(fails); c.cures += len(cures)
        c.restart_same_step += len(np.intersect1d(starts, fails))
        c.start_while_on_tx += len(np.setdiff1d(starts[self.prev_tx[starts]], fails)) # On treatment at the step's start, and didn't fail first
        entered_sym = np.concatenate([get('symptoms'), fails, get('relapse')])
        c.start_not_symptomatic += len(np.setdiff1d(starts[~self.prev_sym[starts]], entered_sym))
        rel_now = np.asarray(tb.relapsing.uids)
        new_rel = rel_now[~self.prev_rel[rel_now]]
        c.relapse_not_from_success += len(np.setdiff1d(new_rel[~self.prev_succ[new_rel]], starts)) # (a start and cure in one step is allowed)
        c.tb_death_from_success += int(self.prev_succ[get('tb_death')].sum() + self.prev_succ[get('tb_death_tx')].sum())
        c.bg_death_from_success += int(self.prev_succ[get('deaths')].sum())
        c.success_person_steps += int(self.prev_succ.sum())
        c.symptomatic_person_steps += int(self.prev_sym.sum())

        # Course lengths, and the share of courses that will cure by year of starting
        ended = np.concatenate([cures, fails]).astype(int)
        timed = ended[~np.isnan(self.ti_start[ended])]
        steps = (self.ti - self.ti_start[timed]).astype(int)
        for s in steps:
            self.course_steps[int(s)] += 1
        c.courses_timed += len(steps); c.course_steps_sum += int(steps.sum()); c.course_steps_sq += int((steps**2).sum())
        self.ti_start[ended] = np.nan
        self.ti_start[starts] = self.ti
        in_succ, in_fail = tb.treatment_success[ss.uids(starts)], tb.treatment_fail[ss.uids(starts)]
        year = int(np.floor(self.sim.t.yearvec[self.ti]))
        self.by_year[year] += [int(in_succ.sum() + in_fail.sum()), int(in_succ.sum())] # Those who already left in the step are left out
        self.alive[self.ti] = self.sim.people.alive.sum()
        self._remember()
        return


def tx_outputs(sim, incidence=None):
    """
    Per step, per person per year (as shares of the living): prevalence, incidence (with relapses), notifications, the share on
    treatment, TB mortality (with deaths on treatment), and their parts. `incidence` is the agents' Incidence analyzer, or for the
    ODE the output of ode_incidence()
    """
    r = sim.results.tb
    alive = np.asarray(sim.results.alive.values, dtype=float)
    dt = sim.t.dt_year
    v = lambda k: np.asarray(r[k].values, dtype=float) if k in r else np.zeros(len(alive))
    out = sc.objdict(
        year         = np.asarray(sim.t.yearvec, dtype=float),
        prevalence   = (v('asymptomatic') + v('symptomatic'))/alive,
        notifications= v('tx_start')/alive/dt,
        on_treatment = (v('treatment_success') + v('treatment_fail'))/alive,
        mortality    = (v('tb_death') + v('tb_death_tx'))/alive/dt,
        mortality_tx = v('tb_death_tx')/alive/dt,
        relapse      = v('relapse')/alive/dt,
    )
    if incidence is not None:
        first = incidence.first if hasattr(incidence, 'first') else incidence
        out.incidence = np.asarray(first, dtype=float)/alive/dt
    return out


def ode_incidence(sim, non0):
    """
    Incidence (first entry to asymptomatic TB per episode, plus relapses) from a run of v4's ODE of TBTreat, per step (people)

    v4's model has one non-infectious state, so its first entries are rebuilt: NON0, the non-infectious people not yet asymptomatic in
    this episode, has the same constant exit rate r = non_rec + non_asy + mu as all non-infectious people, and is entered only from
    latent infection. Each step's inflow (the step's latent_non + remote_non) is taken as uniform within the step, and NON0 is carried
    forward exactly; worsening from NON0 is non_asy × NON0's person-time in the step.

    Args:
        sim: a run of TBTreat with solver='ode'
        non0 (float): NON0 at the start, in people
    """
    tb = sim.diseases.tb
    r = sim.results.tb
    get = lambda k: np.asarray(r[k].values, dtype=float) if k in r else np.zeros(len(r['symptomatic'].values))
    rate =lambda name: float(getattr(tb, name).rate.raw.value) # Constant rates per year
    non_rec, non_asy = rate('recovery'), rate('worsening')
    mu = rate('deaths') if 'deaths' in tb.transitions else 0.0
    rr, dt = non_rec + non_asy + mu, sim.t.dt_year
    inflow = get('latent_non') + get('remote_non')
    decay = np.exp(-rr*dt)
    a = (1 - decay)/rr # Person-time per person present at the step's start
    b = (dt - a)/(rr*dt) # Person-time per person entering uniformly over the step
    n0 = np.zeros(len(inflow)); first_non = np.zeros(len(inflow))
    n0[0] = non0
    for i in range(1, len(inflow)): # Index i holds the step ending at i
        first_non[i] = non_asy*(n0[i-1]*a + inflow[i]*b)
        n0[i] = n0[i-1]*decay + inflow[i]*(1 - decay)/(rr*dt)
    first = get('latent_asy') + get('remote_asy') + first_non + get('relapse')
    return sc.objdict(first=first, non0=n0, first_non=first_non)
