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

__all__ = ['TBTreat', 'phase_lengths', 'outcome_shares', 'carve_out', 'theta_from_coverage', 'lai_tpt_treatment', 'Incidence', 'TxChecks', 'ode_incidence', 'tx_outputs']

tx_state_names = ['treatment_success', 'treatment_fail', 'relapsing']


tb_disease_states = ['non_infectious', 'asymptomatic', 'symptomatic'] # Starting states with TB disease: a course can cure them
latent_states = ['infection', 'infection_remote'] # Starting states with TB infection: a course can clear it (to treated, without relapse)


def _is_zero(x):
    return isinstance(x, (int, float, np.integer, np.floating)) and float(x) == 0.0


def _rate(x, scale=1.0):
    """
    A rate per year for a transition: a number, a dict of knots {year: value} (linear between them, constant outside), or a Starsim rate
    or parameter (only unscaled). `scale` multiplies it (e.g. q·h for loss to follow-up with cure)
    """
    if isinstance(x, dict):
        return ss.peryear(ss.timeseries({float(t): scale*float(v) for t, v in x.items()}, interpolation='linear'))
    if isinstance(x, (int, float, np.integer, np.floating)):
        return ss.peryear(scale*float(x))
    if scale != 1.0:
        errormsg = f'A rate given as {type(x).__name__} can\'t be scaled here; give a number or a dict of knots {{year: value}}'
        raise TypeError(errormsg)
    return x


def _prob(x):
    """ A probability: a number, a dict of knots {year: value} (linear), or an ss.timeseries """
    if isinstance(x, dict):
        return ss.timeseries({float(t): float(v) for t, v in x.items()}, interpolation='linear')
    return x


class TBTreat(tbsim.TB):
    """
    tbsim.TB with treatment as declared states and transitions (exp 23; CONVENTIONS.md 9), with death and loss to follow-up on treatment (exp 35).

    People start treatment from each state in `tx_from` at its rate. From a state with TB disease (non-infectious, asymptomatic or symptomatic
    TB), a course is drawn at its start onto the cure track (probability `tx_success`, the cure share among people who complete it) or the
    fail track. The course lasts `tx_dur`, in phases of `tx_phase_dur` (one phase by default), each a state of fixed length. In each phase:

    - **Death on treatment from TB:** at `tx_death_fail` on the fail track (default `sym_dead`, as exp 23) and `tx_death_cure` on the cure track
    - **Loss to follow-up:** at `tx_ltfu` on both tracks; on the cure track, a share `tx_ltfu_cure` of those lost (by phase) is cured
    - **End of the phase:** the next phase; after the last, the cure track is cured (to `treated`, or `relapsing` with share `tx_p_relapse`,
      then symptomatic TB at rate 1/`dur_relapse_tx`) and the fail track returns to its starting state

    Everyone else leaving treatment returns to the state they started from. TB death on treatment applies by starting state: `tx_death_fail` and
    `tx_death_cure` are symptomatic TB's, unless given as {starting state: rate}; other states with TB disease have none unless given.

    From latent infection (`infection`, `infection_remote`), a course is drawn onto the clearance track (probability `tx_clear_latent`, a
    regimen's efficacy as in a trial) or not. At the end of the course the clearance track goes to `treated` (as previously treated, but never
    to `relapsing`), the other back to its latent state; loss to follow-up at `tx_ltfu`, the clearance track's losses cleared with the shares
    `tx_ltfu_cure`; no TB death. From any other state without TB disease (e.g. susceptible: a false positive), a course has one track and no
    cure: no TB death, loss to follow-up at `tx_ltfu`, and everyone returns to their starting state. Nobody on treatment transmits or is
    reinfected, and natural history pauses during the course. Background deaths (`mu`, or a demography module) apply in every state.

    With the defaults, this is exp 23's model: from symptomatic TB only, one phase, TB death at `sym_dead` on failing courses, no death on the
    cure track and no loss to follow-up. State and transition names are then exp 23's (`treatment_success`, `treatment_fail`, `tx_start`,
    `tx_cure`, `tx_fail`, `tb_death_tx`); otherwise they are generated, and listed in `self.tx_info`.

    Args:
        tx_rate (float/ss.timeseries): the rate of starting treatment from symptomatic TB per year (a timeseries of ss.peryear values to vary in time); used if tx_from is None
        tx_success (float/ss.timeseries/dict): the chance a course is on the cure track, set when it starts (a dict: knots {year: value})
        tx_dur (float): the course length (years)
        course (str): 'fixed' (exact; a conveyor of whole steps in the ODE) or 'rate' (one phase only)
        stages (int): with course='fixed' and one phase, run the ODE's course as this many Erlang stages instead (agents keep the exact course)
        tx_p_relapse (float): the share of cures that relapse
        dur_relapse_tx (float): the mean time from cure to relapse (years)
        tx_phase_dur (float/list): the length of each phase (years), tx_dur being a whole number of them; or a list of the phases' lengths, whose sum
            is then the course (e.g. whole steps, from phase_lengths(): v4's ODE runs each fixed length as whole steps). None: one phase
        tx_from (dict): {starting state: rate of starting treatment}, each rate as tx_rate; None: {'symptomatic': tx_rate}
        tx_death_fail (rate/list/dict): the TB death rate on the fail track from symptomatic TB: a rate for every phase or a list by phase (each a
            number, a dict of knots {year: value}, or a Starsim rate); None: sym_dead. Or {starting state: one of these} for each state with TB
            disease (those not given have none)
        tx_death_cure (rate/list/dict): the TB death rate on the cure track, as tx_death_fail; default 0
        tx_clear_latent (float/dict): the chance a course from latent infection clears it (a number, or knots {year: value}); needed if tx_from has a latent state
        tx_ltfu (rate): the loss-to-follow-up rate, both tracks and every phase (a number or a dict of knots); default 0
        tx_ltfu_cure (float/list/func): the chance that someone on the cure track lost in a phase is cured: a number, a list by phase, or a
            function of the fraction of the course completed by the middle of the phase
        kwargs (dict): passed to tbsim.TB
    """
    def __init__(self, tx_rate=0.0, tx_success=1.0, tx_dur=0.5, course='fixed', stages=None, tx_p_relapse=0.0, dur_relapse_tx=1.5, tx_phase_dur=None,
                 tx_from=None, tx_death_fail=None, tx_death_cure=0.0, tx_ltfu=0.0, tx_ltfu_cure=0.0, tx_clear_latent=None, **kwargs):
        for key in ['theta', 'phi', 'p_relapse']:
            if key in kwargs:
                errormsg = f'TBTreat replaces tbsim.TB\'s own treatment, so {key}= is not used; give tx_rate, tx_success, tx_p_relapse instead'
                raise ValueError(errormsg)
        if course not in ['fixed', 'rate']:
            raise ValueError(f'course is "fixed" or "rate", not "{course}"')
        if isinstance(tx_phase_dur, (list, tuple, np.ndarray)): # The phases' lengths, e.g. whole steps from phase_lengths(); the course is their sum
            widths = [float(x) for x in tx_phase_dur]
            tx_dur = sum(widths)
        else:
            n = 1 if tx_phase_dur is None else int(round(tx_dur/tx_phase_dur))
            if n < 1 or (tx_phase_dur is not None and abs(n*tx_phase_dur - tx_dur) > 1e-9):
                raise ValueError(f'tx_dur ({tx_dur}) must be a whole number of phases of tx_phase_dur ({tx_phase_dur}); or give the phases\' lengths as a list')
            widths = [tx_dur/n]*n
        M = len(widths)
        if M > 1 and course == 'rate':
            raise ValueError('Phases need a fixed course (course="fixed")')
        if M > 1 and stages:
            raise ValueError('Erlang stages (stages=) apply to a course of one phase only')
        super().__init__(**kwargs)
        sym_dead = self.tb_death.rate # The sym_dead parameter itself, so a change to it applies to deaths on treatment too
        ends = np.cumsum(widths)
        mid = [float((e - wk/2)/tx_dur) for e, wk in zip(ends, widths)] # The fraction of the course completed by the middle of each phase

        def per_phase(x, name):
            if callable(x) and not isinstance(x, (ss.timeseries,)) and not hasattr(x, 'raw') and not isinstance(x, dict):
                return [x(f) for f in mid]
            if isinstance(x, (list, tuple)):
                if len(x) != M:
                    raise ValueError(f'{name} has {len(x)} values, for {M} phases')
                return list(x)
            return [x]*M

        by_origin = lambda x: isinstance(x, dict) and len(x) > 0 and all(isinstance(k, str) for k in x) # {starting state: rate}, not knots {year: value}
        def for_origin(x, origin, default_sym):
            """ A TB death rate for courses from `origin`, by phase: from {starting state: rate}, or symptomatic TB's (none for other states) """
            if by_origin(x):
                return per_phase(x.get(origin, 0.0), f'rate for {origin}')
            return per_phase(default_sym if x is None else x, 'rate') if origin == 'symptomatic' else [0.0]*M
        ltfu_cure = per_phase(tx_ltfu_cure, 'tx_ltfu_cure')
        origins = {'symptomatic': tx_rate} if tx_from is None else dict(tx_from)
        if any(o in latent_states for o in origins) and tx_clear_latent is None:
            raise ValueError('Treatment from latent infection needs tx_clear_latent, the chance a course clears it')
        dur_of = lambda wk: ss.constant(ss.years(wk)) if course == 'fixed' else ss.years(wk) # A plain duration is the mean of an exponential
        staged = dict(stages=stages) if (stages and course == 'fixed') else {}
        info = sc.objdict(phases=M, phase_dur=widths, phase_mid=mid, origins=list(origins), cure=[], fail=[], notb=[], clear_track=[], noclear_track=[], origin_of={}, start=[],
                          phase=[], tx_cure=[], tx_fail=[], tx_end=[], tx_clear=[], tb_death_fail=[], tb_death_cure=[], ltfu_cure=[], ltfu_clear=[], ltfu_back=[])
        T = self.define_transitions
        deaths = [] # (name, rate, state) of each death on treatment, for births replacing them
        rel = [] if 'relapsing' in self.states else [ss.State(name='relapsing', label='Cured, will relapse')] # Declared with the first treatment states, as exp 23
        for origin, rate in origins.items():
            src = getattr(self, origin)
            legacy = origin == 'symptomatic' and M == 1 # Exp 23's names
            suffix = '' if origin == 'symptomatic' else f'_{origin}'
            ph = lambda k: '' if M == 1 else f'_{k + 1}'
            start_name = 'tx_start' if origin == 'symptomatic' else f'tx_start{suffix}'
            r = rate if (isinstance(rate, ss.timeseries) or not isinstance(rate, (int, float, dict))) else _rate(rate)
            if origin in tb_disease_states:
                death_fail = for_origin(tx_death_fail, origin, sym_dead)
                death_cure = for_origin(tx_death_cure, origin, 0.0)
                cure = [ss.State(name=f'treatment_success{suffix}{ph(k)}', label=f'On treatment from {origin}, will be cured, phase {k + 1}') for k in range(M)]
                fail = [ss.State(name=f'treatment_fail{suffix}{ph(k)}', label=f'On treatment from {origin}, will fail, phase {k + 1}') for k in range(M)]
                self.define_states(*cure, *fail, *rel)
                rel = []
                cure = [getattr(self, s.name) for s in cure]
                fail = [getattr(self, s.name) for s in fail]
                info.cure += [s.name for s in cure]; info.fail += [s.name for s in fail]
                info.origin_of.update({s.name: origin for s in cure + fail})
                T(**{start_name: ss.Transition(src, {cure[0]: _prob(tx_success), fail[0]: None}, rate=r)})
                info.start.append(start_name)
                for track in [cure, fail]:
                    for k in range(M - 1):
                        name = f'tx_next_{track[k].name}'
                        T(**{name: ss.Transition(track[k], track[k + 1], dur=dur_of(widths[k]))})
                        info.phase.append(name)
                name = 'tx_cure' if legacy else f'tx_cure{suffix}'
                T(**{name: ss.Transition(cure[-1], {self.treated: None, self.relapsing: tx_p_relapse}, dur=dur_of(widths[-1]), **staged)})
                info.tx_cure.append(name)
                name = 'tx_fail' if legacy else f'tx_fail{suffix}'
                T(**{name: ss.Transition(fail[-1], src, dur=dur_of(widths[-1]), **staged)})
                info.tx_fail.append(name)
                for k in range(M):
                    if not _is_zero(death_fail[k]):
                        name = 'tb_death_tx' if legacy else f'tb_death_{fail[k].name}'
                        T(**{name: ss.death(fail[k], rate=_rate(death_fail[k]))})
                        info.tb_death_fail.append(name)
                        deaths.append((name, _rate(death_fail[k]), fail[k]))
                    if not _is_zero(death_cure[k]):
                        name = f'tb_death_{cure[k].name}'
                        T(**{name: ss.death(cure[k], rate=_rate(death_cure[k]))})
                        info.tb_death_cure.append(name)
                        deaths.append((name, _rate(death_cure[k]), cure[k]))
                if not _is_zero(tx_ltfu):
                    for k in range(M):
                        q = float(ltfu_cure[k])
                        if q > 0: # Lost and cured, at q·h; lost and back, at (1 − q)·h: the same as one exit at h split q : 1 − q
                            name = f'tx_ltfu_cured_{cure[k].name}'
                            T(**{name: ss.Transition(cure[k], {self.treated: None, self.relapsing: tx_p_relapse}, rate=_rate(tx_ltfu, q))})
                            info.ltfu_cure.append(name)
                        if q < 1:
                            name = f'tx_ltfu_{cure[k].name}'
                            T(**{name: ss.Transition(cure[k], src, rate=_rate(tx_ltfu, 1 - q))})
                            info.ltfu_back.append(name)
                        name = f'tx_ltfu_{fail[k].name}'
                        T(**{name: ss.Transition(fail[k], src, rate=_rate(tx_ltfu))})
                        info.ltfu_back.append(name)
            elif origin in latent_states: # Latent infection: cleared (to treated, never relapsing) or not; no TB death
                clear = [ss.State(name=f'treatment_clear{suffix}{ph(k)}', label=f'On treatment from {origin}, will clear, phase {k + 1}') for k in range(M)]
                noclear = [ss.State(name=f'treatment_noclear{suffix}{ph(k)}', label=f'On treatment from {origin}, will not clear, phase {k + 1}') for k in range(M)]
                self.define_states(*clear, *noclear)
                clear = [getattr(self, x.name) for x in clear]
                noclear = [getattr(self, x.name) for x in noclear]
                info.clear_track += [x.name for x in clear]; info.noclear_track += [x.name for x in noclear]
                info.origin_of.update({x.name: origin for x in clear + noclear})
                T(**{start_name: ss.Transition(src, {clear[0]: _prob(tx_clear_latent), noclear[0]: None}, rate=r)})
                info.start.append(start_name)
                for track in [clear, noclear]:
                    for k in range(M - 1):
                        name = f'tx_next_{track[k].name}'
                        T(**{name: ss.Transition(track[k], track[k + 1], dur=dur_of(widths[k]))})
                        info.phase.append(name)
                name = f'tx_clear{suffix}'
                T(**{name: ss.Transition(clear[-1], self.treated, dur=dur_of(widths[-1]), **staged)})
                info.tx_clear.append(name)
                name = f'tx_noclear{suffix}'
                T(**{name: ss.Transition(noclear[-1], src, dur=dur_of(widths[-1]), **staged)})
                info.tx_end.append(name)
                if not _is_zero(tx_ltfu):
                    for k in range(M):
                        q = float(ltfu_cure[k])
                        if q > 0:
                            name = f'tx_ltfu_cleared_{clear[k].name}'
                            T(**{name: ss.Transition(clear[k], self.treated, rate=_rate(tx_ltfu, q))})
                            info.ltfu_clear.append(name)
                        if q < 1:
                            name = f'tx_ltfu_{clear[k].name}'
                            T(**{name: ss.Transition(clear[k], src, rate=_rate(tx_ltfu, 1 - q))})
                            info.ltfu_back.append(name)
                        name = f'tx_ltfu_{noclear[k].name}'
                        T(**{name: ss.Transition(noclear[k], src, rate=_rate(tx_ltfu))})
                        info.ltfu_back.append(name)
            else: # No TB disease or infection: one track, no cure, no TB death; everyone returns to the starting state
                notb = [ss.State(name=f'treatment_notb{suffix}{ph(k)}', label=f'On treatment from {origin} (no TB disease), phase {k + 1}') for k in range(M)]
                self.define_states(*notb)
                notb = [getattr(self, s.name) for s in notb]
                info.notb += [s.name for s in notb]
                info.origin_of.update({s.name: origin for s in notb})
                T(**{start_name: ss.Transition(src, notb[0], rate=r)})
                info.start.append(start_name)
                for k in range(M - 1):
                    name = f'tx_next_{notb[k].name}'
                    T(**{name: ss.Transition(notb[k], notb[k + 1], dur=dur_of(widths[k]))})
                    info.phase.append(name)
                name = f'tx_end{suffix}'
                T(**{name: ss.Transition(notb[-1], src, dur=dur_of(widths[-1]), **staged)})
                info.tx_end.append(name)
                if not _is_zero(tx_ltfu):
                    for k in range(M):
                        name = f'tx_ltfu_{notb[k].name}'
                        T(**{name: ss.Transition(notb[k], src, rate=_rate(tx_ltfu))})
                        info.ltfu_back.append(name)
        if rel:
            self.define_states(*rel)
        self.define_transitions(
            relapse = ss.Transition(self.relapsing, self.symptomatic, dur=ss.years(dur_relapse_tx)),
        )
        if isinstance(getattr(self, 'replace', None), ss.Transition): # Births replace TB deaths on treatment too, as tbsim.TB's replace does those from symptomatic TB
            for name, rate, state in deaths:
                rname = 'replace_tx' if name == 'tb_death_tx' else f'replace_{name}'
                T(**{rname: ss.unborn >> self.susceptible @ rate*state.n})
        info.states = info.cure + info.fail + info.clear_track + info.noclear_track + info.notb
        info.deaths = info.tb_death_fail + info.tb_death_cure
        self.tx_info = info
        self.tx = sc.objdict(tx_rate=tx_rate, tx_success=tx_success, tx_dur=tx_dur, course=course, stages=stages, tx_p_relapse=tx_p_relapse,
                             dur_relapse_tx=dur_relapse_tx, tx_phase_dur=tx_phase_dur, tx_from=tx_from, tx_death_fail=tx_death_fail,
                             tx_death_cure=tx_death_cure, tx_ltfu=tx_ltfu, tx_ltfu_cure=tx_ltfu_cure, tx_clear_latent=tx_clear_latent)
        return


def phase_lengths(dur, width, dt):
    """
    A course of `dur` years in phases of about `width` years, each a whole number of steps of `dt` years, adding up to the course in whole steps
    (nearest): e.g. 6 months in 2-month phases at weekly steps is 9 + 9 + 8 weeks (26 weeks, as v4's ODE runs a 6-month course). The phases
    differ by at most a step, the longer first. v4's ODE runs each fixed length as whole steps, so phases not in whole steps would each be rounded
    (3 × 9 weeks, 27 in all); give agents the same lengths.

    Returns:
        list of the phases' lengths in years
    """
    steps = int(round(dur/dt))
    n = max(1, int(round(dur/width)))
    base, extra = divmod(steps, n)
    return [(base + (1 if k < extra else 0))*dt for k in range(n)]


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
    tx_rate = ss.timeseries({float(t): ss.peryear(float(theta_from_coverage(c))) for t, c in zip(times, covg)}, interpolation='linear')
    tx_success = ss.timeseries({t0: 0.0, 1995.0: r['efficacy_1995'], 2024.0: r['efficacy_2024']}, interpolation='linear')
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


def _tx(module, key, legacy):
    """ Names of the treatment states or transitions of a kind from module.tx_info (TBTreat), or the legacy names """
    info = getattr(module, 'tx_info', None)
    return list(info[key]) if info is not None else list(legacy)


def _uids_all(module, names):
    """ The agents who made any of these transitions in this step """
    parts = [_uids(module, n) for n in names]
    return np.concatenate(parts).astype(int) if parts else np.array([], dtype=int)


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
        flagged = ['asymptomatic', 'symptomatic'] + _tx(tb, 'cure', ['treatment_success']) + _tx(tb, 'fail', ['treatment_fail']) # Episodes under way at the start
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
        off = np.concatenate([get('recovery'), _uids_all(tb, _tx(tb, 'tx_cure', ['tx_cure']) + _tx(tb, 'ltfu_cure', []))]).astype(int) # The episode ends (at any cure)
        self._grow(np.concatenate([on, off, [0]]).astype(int))
        self.flag[on] = True
        self.flag[off] = False
        return


class TxChecks(ss.Analyzer, _Grow):
    """
    Agents' checks of the treatment pitfalls in exp 23's README, and of exp 35's death and loss to follow-up, from each step's transitions and the
    states at the end of the last step. Counts are summed over the run; `by_year` holds the courses started each year and the share that will cure.

    Exp 23's checks (each should be 0, or rare where noted in its README): a start while on treatment (start_while_on_tx), a start from symptomatic
    TB by someone not symptomatic (start_not_symptomatic), relapsing entered other than from the cure track (relapse_not_from_success), and TB
    death at sym_dead or on the fail track among people on the cure track (tb_death_from_success). Exp 35's: anyone on treatment from a state
    without TB disease cured (notb_cured) or dying of TB (tb_death_notb); and of those leaving treatment other than by cure or death, how many
    are in their starting state at the end of the step (left_to_origin) or elsewhere (left_elsewhere: should be rare, from a further move in the same step).
    """
    def init_post(self):
        super().init_post()
        tb = self.sim.diseases.tb
        self.tx_states = _tx(tb, 'states', ['treatment_success', 'treatment_fail'])
        self.cure_states = _tx(tb, 'cure', ['treatment_success'])
        self.notb_states = _tx(tb, 'notb', [])
        self.latent_tx_states = _tx(tb, 'clear_track', []) + _tx(tb, 'noclear_track', [])
        info = getattr(tb, 'tx_info', None)
        self.origins = list(info.origins) if info is not None else ['symptomatic']
        self.origin_of = {k: self.origins.index(v) for k, v in info.origin_of.items()} if info is not None else {'treatment_success': 0, 'treatment_fail': 0}
        self._arrays = ['ti_start', 'prev_tx', 'prev_succ', 'prev_sym', 'prev_rel', 'prev_notb', 'prev_lat', 'origin']
        self._arr('ti_start', float, np.nan)
        for name in ['prev_tx', 'prev_succ', 'prev_sym', 'prev_rel', 'prev_notb', 'prev_lat']:
            self._arr(name, bool, False)
        self._arr('origin', int, -1)
        self.c = sc.objdict({k: 0 for k in ['starts', 'fails', 'cures', 'ends', 'restart_same_step', 'start_while_on_tx', 'start_not_symptomatic',
                                            'relapse_not_from_success', 'tb_death_from_success', 'bg_death_from_success', 'success_person_steps',
                                            'symptomatic_person_steps', 'courses_timed', 'course_steps_sum', 'course_steps_sq', 'tx_deaths_fail',
                                            'tx_deaths_cure', 'ltfu', 'ltfu_cured', 'notb_cured', 'tb_death_notb', 'left_to_origin', 'left_elsewhere', 'latent_cleared',
                                            'latent_to_relapsing', 'tb_death_latent']})
        self.course_steps = sc.ddict(int) # Course lengths in steps, of courses started during the run and completed
        self.by_year = sc.ddict(lambda: np.zeros(2)) # Courses started each year, and those that will cure
        self.alive = np.zeros(self.sim.t.npts)
        self._remember()
        return

    def _remember(self):
        tb = self.sim.diseases.tb
        for name, states in [('prev_tx', self.tx_states), ('prev_succ', self.cure_states), ('prev_sym', ['symptomatic']), ('prev_rel', ['relapsing']),
                             ('prev_notb', self.notb_states), ('prev_lat', self.latent_tx_states)]:
            arr = getattr(self, name)
            arr[:] = False
            for state in states:
                u = np.asarray(getattr(tb, state).uids)
                self._grow(u)
                arr[u] = True
        self.origin[:] = -1
        for state in self.tx_states:
            u = np.asarray(getattr(tb, state).uids)
            self.origin[u] = self.origin_of[state]
        return

    def step(self):
        tb, c = self.sim.diseases.tb, self.c
        get = lambda name: _uids(tb, name)
        many = lambda key, legacy: _uids_all(tb, _tx(tb, key, legacy))
        starts, fails, cures, ends = many('start', ['tx_start']), many('tx_fail', ['tx_fail']), many('tx_cure', ['tx_cure']), many('tx_end', [])
        d_fail, d_cure = many('tb_death_fail', ['tb_death_tx']), many('tb_death_cure', [])
        ltfu_cured, ltfu_back, ltfu_clear, clears = many('ltfu_cure', []), many('ltfu_back', []), many('ltfu_clear', []), many('tx_clear', [])
        allu = np.concatenate([starts, fails, cures, ends, d_fail, d_cure, ltfu_cured, ltfu_back, ltfu_clear, clears, get('tb_death'), get('deaths'), get('symptoms'),
                               get('relapse'), [0]]).astype(int)
        self._grow(allu)
        c.starts += len(starts); c.fails += len(fails); c.cures += len(cures); c.ends += len(ends)
        c.tx_deaths_fail += len(d_fail); c.tx_deaths_cure += len(d_cure); c.ltfu += len(ltfu_cured) + len(ltfu_back) + len(ltfu_clear); c.ltfu_cured += len(ltfu_cured)
        c.latent_cleared += len(clears) + len(ltfu_clear)
        returned = np.concatenate([fails, ends, ltfu_back]).astype(int)
        c.restart_same_step += len(np.intersect1d(starts, returned))
        c.start_while_on_tx += len(np.setdiff1d(starts[self.prev_tx[starts]], returned)) # On treatment at the step's start, and didn't leave first
        sym_starts = get('tx_start')
        entered_sym = np.concatenate([get('symptoms'), fails, ltfu_back, get('relapse')]).astype(int)
        c.start_not_symptomatic += len(np.setdiff1d(sym_starts[~self.prev_sym[sym_starts]], entered_sym))
        rel_now = np.asarray(tb.relapsing.uids)
        self._grow(rel_now)
        new_rel = rel_now[~self.prev_rel[rel_now]]
        c.relapse_not_from_success += len(np.setdiff1d(new_rel[~self.prev_succ[new_rel]], starts)) # (a start and cure in one step is allowed)
        c.tb_death_from_success += int(self.prev_succ[get('tb_death')].sum() + self.prev_succ[d_fail].sum())
        c.bg_death_from_success += int(self.prev_succ[get('deaths')].sum())
        c.success_person_steps += int(self.prev_succ.sum())
        c.symptomatic_person_steps += int(self.prev_sym.sum())
        cured_now = np.concatenate([np.asarray(tb.treated.uids), rel_now]).astype(int)
        self._grow(cured_now)
        c.notb_cured += int(self.prev_notb[cured_now].sum()) # On treatment without TB disease at the step's start, cured (treated or relapsing) at its end
        c.tb_death_notb += int(self.prev_notb[np.concatenate([get('tb_death'), d_fail, d_cure]).astype(int)].sum())
        c.tb_death_latent += int(self.prev_lat[np.concatenate([get('tb_death'), d_fail, d_cure]).astype(int)].sum())
        c.latent_to_relapsing += int(self.prev_lat[new_rel].sum()) # On treatment for latent infection at the step's start, relapsing at its end: never

        # Where those leaving treatment other than by cure or death are at the end of the step
        for k, origin in enumerate(self.origins):
            left = returned[self.origin[returned] == k]
            left = left[np.asarray(self.sim.people.alive.raw[left], dtype=bool)]
            at = np.asarray(getattr(tb, origin)[ss.uids(left)], dtype=bool) if len(left) else np.array([], dtype=bool)
            on_tx = np.zeros(len(left), dtype=bool)
            for state in self.tx_states: # Restarted in the same step: allowed
                if len(left):
                    on_tx |= np.asarray(getattr(tb, state)[ss.uids(left)], dtype=bool)
            c.left_to_origin += int(at.sum() + (on_tx & ~at).sum())
            c.left_elsewhere += int((~at & ~on_tx).sum())

        # Course lengths (completed courses), and the share of courses that will cure by year of starting
        ended = np.concatenate([cures, fails, ends, clears]).astype(int)
        timed = ended[~np.isnan(self.ti_start[ended])]
        steps = (self.ti - self.ti_start[timed]).astype(int)
        for s in steps:
            self.course_steps[int(s)] += 1
        c.courses_timed += len(steps); c.course_steps_sum += int(steps.sum()); c.course_steps_sq += int((steps**2).sum())
        self.ti_start[np.concatenate([ended, clears, d_fail, d_cure, ltfu_cured, ltfu_back, ltfu_clear]).astype(int)] = np.nan
        self.ti_start[starts] = self.ti
        on_cure = np.zeros(len(starts), dtype=bool); on_any = np.zeros(len(starts), dtype=bool)
        for state in self.tx_states:
            v = np.asarray(getattr(tb, state)[ss.uids(starts)], dtype=bool) if len(starts) else np.array([], dtype=bool)
            on_any |= v
            if state in self.cure_states:
                on_cure |= v
        year = int(np.floor(self.sim.t.yearvec[self.ti]))
        self.by_year[year] += [int(on_any.sum()), int(on_cure.sum())] # Those who already left in the step are left out
        self.alive[self.ti] = self.sim.people.alive.sum()
        self._remember()
        return


def tx_outputs(sim, incidence=None):
    """
    Per step, per person per year (as shares of the living): prevalence, incidence (with relapses), notifications (every start of treatment),
    the share on treatment, TB mortality (with deaths on treatment), and their parts: deaths on treatment on the fail and cure tracks, losses to
    follow-up, and those cured on being lost. `incidence` is the agents' Incidence analyzer, or for the ODE the output of ode_incidence().
    Stocks at index i are at t_i; flows at index i are the step from t_i to t_(i+1) (see ode_incidence())
    """
    r = sim.results.tb
    tb = sim.diseases.tb
    alive = np.asarray(sim.results.alive.values, dtype=float)
    dt = sim.t.dt_year
    v = lambda k: np.asarray(r[k].values, dtype=float) if k in r else np.zeros(len(alive))
    total = lambda names: sum((v(k) for k in names), np.zeros(len(alive)))
    d_fail, d_cure = total(_tx(tb, 'tb_death_fail', ['tb_death_tx'])), total(_tx(tb, 'tb_death_cure', []))
    ltfu_cured, ltfu_back, ltfu_clear = total(_tx(tb, 'ltfu_cure', [])), total(_tx(tb, 'ltfu_back', [])), total(_tx(tb, 'ltfu_clear', []))
    out = sc.objdict(
        year         = np.asarray(sim.t.yearvec, dtype=float),
        prevalence   = (v('asymptomatic') + v('symptomatic'))/alive,
        notifications= total(_tx(tb, 'start', ['tx_start']))/alive/dt,
        on_treatment = total(_tx(tb, 'states', ['treatment_success', 'treatment_fail']))/alive,
        mortality    = (v('tb_death') + d_fail + d_cure)/alive/dt,
        mortality_tx = (d_fail + d_cure)/alive/dt,
        mortality_tx_cure = d_cure/alive/dt,
        ltfu         = (ltfu_cured + ltfu_back + ltfu_clear)/alive/dt,
        latent_cleared = (total(_tx(tb, 'tx_clear', [])) + ltfu_clear)/alive/dt,
        ltfu_cured   = ltfu_cured/alive/dt,
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

    Starsim's flow results at index i are the step from t_i to t_(i+1), as are those of agents and the Incidence analyzer (starsim d83f4f87,
    exp 26). At the earlier pins (starsim 80c81f53, exps 19–25), index i held the step ending at t_i, and this function followed that.

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
    for i in range(len(inflow)): # Index i holds the step from t_i to t_(i+1); n0[i] is NON0 at t_i
        first_non[i] = non_asy*(n0[i]*a + inflow[i]*b)
        if i + 1 < len(inflow):
            n0[i+1] = n0[i]*decay + inflow[i]*(1 - decay)/(rr*dt)
    first = get('latent_asy') + get('remote_asy') + first_non + get('relapse')
    return sc.objdict(first=first, non0=n0, first_non=first_non)


#%% Outcome shares of a course, and the carve-out (exp 35)

def outcome_shares(c, mu_f, mu_c, h, q, widths, bg=0.0):
    """
    The shares of people starting a course (from a state with TB disease) by outcome, exactly in continuous time: phases of fixed lengths
    `widths` (years), in each competing exits by TB death (mu_f on the fail track, mu_c on the cure track), background death (bg), and loss to
    follow-up (h); on the cure track, a share q of those lost is cured. mu_f, mu_c and q are numbers or lists by phase.

    Returns:
        objdict: cured (completing on the cure track, or lost and cured), cured_ltfu, died (TB and background), died_tb, died_tb_cure (TB deaths
        on the cure track), died_bg, lost (all lost to follow-up, cured or not: WHO's category), failed (completing on the fail track), and
        omega (the cure track's share of TB deaths)
    """
    M = len(widths)
    per = lambda x: list(x) if isinstance(x, (list, tuple, np.ndarray)) else [x]*M
    def track(mu, qs):
        reach, out = 1.0, dict(died_tb=0.0, died_bg=0.0, lost=0.0, lost_cured=0.0)
        for mk, w, qk in zip(per(mu), widths, qs):
            lam = mk + h + bg
            leave = reach*(-np.expm1(-lam*w))
            if lam > 0:
                out['died_tb'] += leave*mk/lam; out['died_bg'] += leave*bg/lam; out['lost'] += leave*h/lam; out['lost_cured'] += leave*h/lam*qk
            reach *= np.exp(-lam*w)
        out['complete'] = reach
        return out
    cu, fa = track(mu_c, per(q)), track(mu_f, [0.0]*M)
    died_tb = c*cu['died_tb'] + (1 - c)*fa['died_tb']
    return sc.objdict(cured=c*(cu['complete'] + cu['lost_cured']), cured_ltfu=c*cu['lost_cured'], died=died_tb + c*cu['died_bg'] + (1 - c)*fa['died_bg'],
                      died_tb=died_tb, died_tb_cure=c*cu['died_tb'], died_bg=c*cu['died_bg'] + (1 - c)*fa['died_bg'],
                      lost=c*cu['lost'] + (1 - c)*fa['lost'], failed=(1 - c)*fa['complete'], omega=c*cu['died_tb']/died_tb if died_tb > 0 else 0.0)


def carve_out(s, d, ell, widths, q=0.0, omega=0.0, bg=0.0):
    """
    Exp 35's carve-out: the cure track's share c, the fail track's TB death rate mu_f (the same in every phase), the cure track's TB death rate
    in the first phase mu_c (0 if omega = 0) and the loss-to-follow-up rate h, such that a course started from a state with TB disease has overall
    cure s (with partial cures), a died share d (TB and background deaths: WHO's "died"), a lost share ell (WHO's "lost to follow-up"), and the
    cure track a share omega of TB deaths on treatment. Solved exactly (outcome_shares()).

    If the solution has c > 1, c is capped at 1 and mu_c and h are solved for d and ell alone (overall cure then falls short of s); `c_raw`
    is the uncapped value and `capped` says so.

    Returns:
        objdict: c, c_raw, capped, mu_f, mu_c (by phase), h, q (by phase), shares (outcome_shares() at the solution), and the residual
    """
    from scipy.optimize import least_squares
    M = len(widths)
    qs = list(q) if isinstance(q, (list, tuple, np.ndarray)) else [q]*M
    muc = lambda m1: [m1] + [0.0]*(M - 1)
    def shares(x):
        c, lmf, lmc, lh = x
        return outcome_shares(c, np.exp(lmf), muc(np.exp(lmc) if omega > 0 else 0.0), np.exp(lh), qs, widths, bg)
    def resid(x):
        o = shares(x)
        r = [o.cured - s, o.died - d, o.lost - ell]
        if omega > 0:
            r.append(o.omega - omega)
        return r
    x0 = [s, np.log(0.3), np.log(0.3), np.log(max(-np.log1p(-ell)/sum(widths), 1e-6))]
    if omega > 0:
        sol = least_squares(resid, x0, bounds=([0, -20, -20, -20], [5, 5, 5, 5]), xtol=1e-14, ftol=1e-14, gtol=1e-14)
        x = sol.x
    else:
        sol = least_squares(lambda y: resid([y[0], y[1], 0.0, y[2]]), [x0[0], x0[1], x0[3]], bounds=([0, -20, -20], [5, 5, 5]), xtol=1e-14, ftol=1e-14, gtol=1e-14)
        x = [sol.x[0], sol.x[1], -np.inf, sol.x[2]]
    c_raw = float(x[0])
    capped = c_raw > 1
    if capped: # The fail track is empty: deaths only on the cure track, all in the first phase whatever omega
        def r2(y):
            o = outcome_shares(1.0, 0.0, muc(np.exp(y[0])), np.exp(y[1]), qs, widths, bg)
            return [o.died - d, o.lost - ell]
        y = least_squares(r2, [np.log(0.3), x0[3]], bounds=([-20, -20], [5, 5]), xtol=1e-14, ftol=1e-14, gtol=1e-14).x
        c, mu_f, mu_c, h = 1.0, 0.0, muc(float(np.exp(y[0]))), float(np.exp(y[1]))
    else:
        c, mu_f, mu_c, h = c_raw, float(np.exp(x[1])), muc(float(np.exp(x[2])) if omega > 0 else 0.0), float(np.exp(x[3]))
    o = outcome_shares(c, mu_f, mu_c, h, qs, widths, bg)
    res = float(np.max(np.abs([o.cured - s, o.died - d, o.lost - ell] + ([o.omega - omega] if omega > 0 and not capped else []))))
    return sc.objdict(c=c, c_raw=c_raw, capped=capped, mu_f=mu_f, mu_c=mu_c, h=h, q=qs, shares=o, residual=res)
