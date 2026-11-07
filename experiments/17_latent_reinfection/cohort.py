"""
Experiment 17, question 2: the per-person cohort outputs, defined once and used by both environments (numpy and sciris only).

Each model's agents report, every step, which agents had which events; `Tracker` turns them into the outputs. The events reported
at analyzer step ti happen at time ti + offset (in steps): v4's analyzers see the starting state at ti = 0 and the step ending at
time ti after that (offset 0); v0.10's see the step ending at time ti + 1 (offset 1, as tb_LAI_TPT's FirstPassage). The cohort is
infected at time 0.

Events (arrays of agent indices):
- new:      a new infection of someone no longer infected (cleared or recovered), starting a latent spell
- reset:    a reinfection of a latent person that restarts their latent course (v0.10: any latent reinfection, which resets the
            clock; v4: remote -> recent, since reinfecting a recent infection changes nothing)
- progress: latent -> non-infectious or asymptomatic TB (an episode of TB)
- clear:    latent -> cleared
- asy, sym: entering asymptomatic or symptomatic TB

An agent can have more than one event in a step. They are applied in the only order possible: resets (only while latent), new
infections of people who weren't latent at the start of the step, exits from latent infection, then new infections of people who
were latent at the start (they must have left first).
"""
import numpy as np
import sciris as sc

__all__ = ['Tracker', 'event_names']

event_names = ['new', 'reset', 'progress', 'clear', 'asy', 'sym']
max_count = 20 # Counts per person (resets, episodes) are binned up to this, with the last bin holding everything above


class Tracker:
    """
    Per-person cohort outputs over the run.

    Args:
        n (int): agents, all latent at time 0
        npts (int): steps
        dt_year (float): the step, in years
        offset (int): the time, in steps, of the events reported at analyzer step 0 (v4: 0; v0.10: 1)
        infect_at_start (bool): whether new infections and resets happen at the start of their step rather than the end (v0.10's
            TBConstantFOI infects where transmission is, before the step's transitions: True; v4 dates them at the end of the step: False)
    """
    def __init__(self, n, npts, dt_year, offset, infect_at_start=False):
        self.n, self.npts, self.dt_year, self.offset = n, npts, dt_year, offset
        self.infect_at_start = infect_at_start
        self.latent = np.ones(n, bool)
        self.spell_start = np.zeros(n, int)  # When the current latent spell began
        self.last_reset = np.zeros(n, int)   # When the current latent course began (infection or reset)
        self.episodes = np.zeros(n, int)     # Progression events
        self.resets_pre = np.zeros(n, int)   # Resets before the first progression
        self.everA = np.zeros(n, bool)       # Ever progressed (non-infectious TB or worse)
        self.everB = np.zeros(n, bool)       # Ever asymptomatic (or symptomatic) TB
        self.everC = np.zeros(n, bool)       # Ever symptomatic TB
        self.time_latent = np.zeros(n, int)  # Total time latent over completed spells, in steps (open spells are added in results())
        self.n_spells = np.ones(n, int)      # Latent spells begun (the first at time 0)
        self.t = (np.arange(npts) + offset)*dt_year
        self.curves = sc.objdict({k: np.zeros(npts) for k in ['A', 'B', 'C', 'latent']})
        nb = npts + offset + 1
        self.hists = sc.objdict(
            reset_to_first_prog = np.zeros(nb, int), # Time from the last reset (or infection) to the first progression, in steps
            reset_to_prog       = np.zeros(nb, int), # ... to every progression
            spell_prog          = np.zeros(nb, int), # Latent spells ending in progression, by duration in steps
            spell_clear         = np.zeros(nb, int), # ... ending in clearance
        )
        return

    def step(self, ti, ev):
        """ Apply step ti's events (dict of index arrays, keys from event_names; missing keys are empty) """
        ev = {k: np.asarray(ev.get(k, []), dtype=int) for k in event_names}
        now = ti + self.offset
        t_inf = now - 1 if self.infect_at_start else now # When this step's new infections and resets happened
        was_latent = self.latent.copy()

        # Resets: only while latent
        r = ev['reset']
        assert self.latent[r].all(), 'A reset of someone not latent'
        self.last_reset[r] = t_inf
        self.resets_pre[r[self.episodes[r] == 0]] += 1

        # New infections of people not latent at the start of the step, then exits, then new infections of people who just left
        new_first = ev['new'][~was_latent[ev['new']]]
        new_after = ev['new'][was_latent[ev['new']]]
        self._start(new_first, t_inf)
        for key, hist in [('progress', 'spell_prog'), ('clear', 'spell_clear')]:
            u = ev[key]
            assert self.latent[u].all(), f'A latent exit ({key}) of someone not latent'
            np.add.at(self.hists[hist], now - self.spell_start[u], 1)
            self.time_latent[u] += now - self.spell_start[u]
            self.latent[u] = False
        p = ev['progress']
        first = p[self.episodes[p] == 0]
        np.add.at(self.hists.reset_to_first_prog, now - self.last_reset[first], 1)
        np.add.at(self.hists.reset_to_prog, now - self.last_reset[p], 1)
        self.episodes[p] += 1
        self.everA[p] = True
        self._start(new_after, now) # Only possible when infection is dated at the end of the step (v4): they left latent infection first

        self.everB[ev['asy']] = True
        self.everC[ev['sym']] = True
        self.everA |= self.everB # Asymptomatic TB is reached through progression, possibly in the same step
        c = self.curves
        c.A[ti], c.B[ti], c.C[ti], c.latent[ti] = self.everA.mean(), self.everB.mean(), self.everC.mean(), self.latent.mean()
        return

    def _start(self, u, now):
        assert not self.latent[u].any(), 'A new infection of someone already latent'
        self.latent[u] = True
        self.spell_start[u] = now
        self.last_reset[u] = now
        self.n_spells[u] += 1
        return

    def results(self):
        """ The outputs: curves by step, histograms, and per-person counts """
        out = sc.objdict(t=self.t, **self.curves)
        out.hists = sc.objdict({k: v.copy() for k, v in self.hists.items()})
        open_ = np.flatnonzero(self.latent) # Spells still open at the end
        end = self.npts - 1 + self.offset
        out.hists.spell_open = np.bincount(end - self.spell_start[open_], minlength=end + 1)
        total = self.time_latent.copy() # Each person's total time latent over the run, open spells included up to the end
        total[open_] += end - self.spell_start[open_]
        out.hists.time_latent = np.bincount(total, minlength=end + 1)       # People, by total time latent in steps
        out.hists.time_latent_open = np.bincount(total[open_], minlength=end + 1) # ... of whom still latent at the end
        out.n_spells = np.bincount(np.minimum(self.n_spells, max_count), minlength=max_count + 1)
        out.dt_year = self.dt_year
        out.resets_pre = np.bincount(np.minimum(self.resets_pre, max_count), minlength=max_count + 1)
        out.episodes = np.bincount(np.minimum(self.episodes, max_count), minlength=max_count + 1)
        reset = self.resets_pre > 0
        out.split = sc.objdict(n_reset=int(reset.sum()), n_reset_prog=int((reset & self.everA).sum()),
                               n_none=int((~reset).sum()), n_none_prog=int((~reset & self.everA).sum()))
        out.n = self.n
        return out
