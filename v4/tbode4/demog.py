"""
Demography matched between v4's agents and its ODE (CONVENTIONS.md 10; exps 28–29).

- WPP's rates as bundled with Starsim (every 10 years), interpolated linearly by year
- Deaths by age and sex (both solvers), and the ODE's births by age (`ss.birth`)
- Agents' births: `ss.Pregnancy` with fertility shifted younger by its mean gestation, the steady-state share pregnant at the start
  (`PregnancyStart`), and births as female as the ODE's
- The ODE's start: WPP's split by sex within each age group (`joint_sex`), and fine age groups (`age_edges`)

Moved here from exps 28 and 29's run.py, which keep their own copies so they can be rerun as they were.
"""
import numpy as np
import pandas as pd
import starsim as ss

__all__ = ['wpp_edges', 'sexes', 'WPP', 'wpp_rate', 'death_rate', 'ode_births', 'ode_demographics', 'gestation', 'shifted_fertility', 'PregnancyStart',
           'agent_demographics', 'female_share', 'joint_sex', 'age_edges']

wpp_edges = [0, 1] + list(range(5, 101, 5)) # WPP's mortality groups: 0, 1–4, 5–9, ..., 95–99, 100+
sexes = ['female', 'male']


class WPP:
    """ A country's bundled WPP rates as arrays by decade [decade, age group], and their linear interpolation to any year """
    def __init__(self, country='India', decades=None):
        self.country = country
        self.decades = np.arange(1950, 2101, 10) if decades is None else np.asarray(decades)
        mx = ss.data.rates(country, 'deaths')
        fr = ss.data.rates(country, 'fertility')
        self.mx = {s: mx[mx.Sex == s.capitalize()].pivot(index='Time', columns='AgeGrpStart', values='mx').loc[self.decades, wpp_edges].values for s in sexes}
        self.fert_ages = sorted(fr.AgeGrpStart.unique()) # 10, 15, ..., 50
        self.asfr = fr.pivot(index='Time', columns='AgeGrpStart', values='ASFR').loc[self.decades, self.fert_ages].values/1000 # Per woman per year
        return

    def at(self, table, year):
        return np.array([np.interp(year, self.decades, table[:, j]) for j in range(table.shape[1])])

    def mx_at(self, year):
        """ Death rates [sex, WPP age group] in `year` """
        return np.array([self.at(self.mx[s], year) for s in sexes])

    def asfr_at(self, year):
        """ Fertility rates by 5-year group of age at delivery (10–14, ..., 50–54) in `year` """
        return self.at(self.asfr, year)

    def asfr_by_age(self, values, ages):
        """ One year's rates (`values`, by WPP's groups) at each age in `ages`, as a step function; 0 outside 10–54 """
        fa = np.array(self.fert_ages, dtype=float)
        j = np.clip(np.digitize(ages, fa) - 1, 0, len(fa) - 1)
        return np.where((ages >= 10) & (ages < 55), values[j], 0.0)


def wpp_rate(wpp, values, frozen=False, frozen_year=None):
    """ A rate per year over time: a timeseries on WPP's decades (linear between them), or frozen at its value in frozen_year (default the first decade) """
    if frozen:
        y = wpp.decades[0] if frozen_year is None else frozen_year
        return ss.peryear(float(np.interp(y, wpp.decades, values)))
    return ss.peryear(ss.timeseries({float(y): float(v) for y, v in zip(wpp.decades, values)}, interpolation='linear'))


def death_rate(wpp, frozen=False):
    """ ss.death of everyone, by WPP's age groups and sex (agents and the ODE) """
    rows = [dict(age=a, sex=s, value=wpp_rate(wpp, wpp.mx[s][:, j], frozen)) for s in sexes for j, a in enumerate(wpp_edges)]
    return ss.death(rate=ss.by(pd.DataFrame(rows)))


def ode_births(wpp, frozen=False):
    """ The ODE's births: each woman's age-specific fertility, by her (WPP) age group; none for men. Compartmental only """
    zero = ss.peryear(0.0)
    fert = {a: wpp_rate(wpp, wpp.asfr[:, j], frozen) for j, a in enumerate(wpp.fert_ages)}
    rows = [dict(age=a, sex='female', value=fert.get(a, zero)) for a in wpp_edges] + [dict(age=a, sex='male', value=zero) for a in wpp_edges]
    return ss.birth(rate=ss.by(pd.DataFrame(rows)))


def ode_demographics(wpp, frozen=False):
    """ The ODE's births and deaths. Stratify the sim by age (e.g. age_edges(0.25)) and sex (joint_sex) """
    return [ode_births(wpp, frozen), death_rate(wpp, frozen)]


def gestation():
    """ ss.Pregnancy's mean duration of pregnancy, in years (0.769 at the default dur_pregnancy) """
    sim = ss.Sim(diseases=ss.SIS(beta=0, init_prev=0), demographics=[ss.Pregnancy(fertility_rate=30, rate_per=1000)], n_agents=100, dt=ss.months(1), start=2000, stop=2001, verbose=0)
    sim.init()
    return float(sim.demographics.pregnancy.mean_dur_pregnancy)


def shifted_fertility(wpp, g=None, grid=0.25, years=None):
    """ A fertility table for ss.Pregnancy (per 1,000 women): conception at age x takes WPP's rate for age at delivery x + g, on a grid of `grid` years,
    one row per year. Returns the table and the ages at conception it covers (min_age, max_age) """
    g = gestation() if g is None else g
    years = np.arange(int(wpp.decades[0]), int(wpp.decades[-1]) + 1) if years is None else years
    ages = np.round(np.arange(10 - g, 55 - g, grid), 6)
    rows = []
    for y in years:
        vals = wpp.asfr_by_age(wpp.asfr_at(y), ages + g)
        rows += [dict(Time=y, AgeGrp=a, ASFR=1000*v) for a, v in zip(ages, vals)]
    return pd.DataFrame(rows), float(ages[0]), float(ages[-1] + grid)


class PregnancyStart(ss.Pregnancy):
    """ ss.Pregnancy with the steady-state share pregnant at the start: each woman pregnant with probability h·d (her conception rate × the
    mean gestation, unrounded), instead of burnin_plan()'s 1 − exp(−h·steps·dt), with steps·dt the mean gestation rounded to whole steps """
    def burnin_plan(self, engine):
        plan = super().burnin_plan(engine)
        if plan is None:
            return plan
        orig, steps = plan.share, plan.steps
        factor = ss.solvers.convert(ss.years(self.mean_dur_pregnancy), engine.sim.unit_base)/engine.dt/steps # d ÷ (steps·dt)
        plan.share = lambda: np.minimum(-np.log1p(-np.minimum(orig(), 1 - 1e-12))*factor, 1.0) # orig = 1 − exp(−h·steps·dt), so this is h·d
        return plan


def female_share(country='India', year=1950):
    """ The bundled population in `year`: its groups' lower edges, the female share in each, and overall """
    pyr = ss.data.population(country, year=year)
    lows = np.array(ss.data.age_edges(pyr.index), dtype=float)
    return lows, pyr.female.values/(pyr.female.values + pyr.male.values), float(pyr.female.sum()/(pyr.female.sum() + pyr.male.sum()))


def agent_demographics(wpp, year=1950, p_female=None):
    """ Agents' births and deaths: PregnancyStart with the shifted table, births as female as the ODE's started with joint_sex (the population's
    female share overall in `year`, since the ODE's births take its starting share), and the ODE's deaths """
    tab, lo, hi = shifted_fertility(wpp)
    p_female = female_share(wpp.country, year)[2] if p_female is None else p_female
    return [PregnancyStart(fertility_rate=tab, min_age=lo, max_age=hi, p_female=ss.bernoulli(p_female), name='pregnancy'), death_rate(wpp)]


def joint_sex(edges, country='India', year=1950):
    """ The ODE's starting split by sex within each of its age groups (`edges`), from the bundled population; each group takes the share of the
    population group holding it. Pass as sim.stratify(age=edges, sex=joint_sex(edges)). The ODE's births then take its starting female share overall """
    lows, share, _ = female_share(country, year)
    f = share[np.clip(np.digitize(edges, lows) - 1, 0, len(lows) - 1)]
    return ss.levels(sexes, fractions=ss.by(pd.DataFrame({'female': f, 'male': 1 - f}, index=list(edges)), dims=['age', 'sex']))


def age_edges(width, top=100):
    """ Age groups of `width` years from 0 to `top` (the last, top+, open) """
    return [float(a) for a in np.round(np.arange(0, top + 1e-9, width), 6)]
