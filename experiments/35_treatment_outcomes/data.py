"""
Experiment 35, stage 0: India's treatment outcomes from WHO's TB database, as shares of evaluated cases, and the died and lost-to-follow-up
shares the model uses (fitted_schedule(): all new cases; died constant, lost linear between fitted knots). See README.md. Run from the repo root:

    v4/.venv/bin/python -I experiments/35_treatment_outcomes/data.py           # uses data/TB_outcomes.csv if present
    v4/.venv/bin/python -I experiments/35_treatment_outcomes/data.py download  # fetch WHO's CSVs again first
    v4/.venv/bin/python -I experiments/35_treatment_outcomes/data.py options   # yearly and smoothed options (outcome_options.png), not used
    v4/.venv/bin/python -I experiments/35_treatment_outcomes/data.py fitted    # the schedules used (outcome_schedule.csv, .png)

Cohorts: new smear-positive cases to 2011 (new_sp_*: cured + completed = success; defaulted = lost to follow-up), new and relapse cases from 2012
(newrel_*). Not evaluated (the cohort less every outcome) is left out of the denominator (user, 2026-11-23). These are notified, mostly public-sector
cohorts; private-sector outcomes are largely missing (README).

Knot values are means of the yearly shares over a window around each knot (`windows`): 1995 from 1999–2001, because before 1999 most of each cohort
wasn't evaluated during the rollout of DOTS and the died shares are implausibly small; 2010 from 2008–2012; 2024 from 2021–2023 (the latest cohorts).
"""
import sys
import urllib.request
import numpy as np
import pandas as pd
import sciris as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

here = sc.thispath(__file__)
datadir = here / 'data'
urls = dict(outcomes='https://extranet.who.int/tme/generateCSV.asp?ds=outcomes', dictionary='https://extranet.who.int/tme/generateCSV.asp?ds=dictionary')
windows = sc.objdict({1995: (1999, 2001), 2010: (2008, 2012), 2024: (2021, 2023)}) # Knot: years averaged (inclusive)


def download():
    datadir.mkdir(exist_ok=True)
    for name, url in urls.items():
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req) as r, open(datadir / f'TB_{name if name == "outcomes" else "data_dictionary"}.csv', 'wb') as fh:
            fh.write(r.read())
    return


def yearly_outcomes(iso3='IND'):
    """ Each cohort year's outcomes: counts, the evaluated share of the cohort, and success, failed, died and lost as shares of those evaluated """
    d = pd.read_csv(datadir / 'TB_outcomes.csv')
    c = d[d.iso3 == iso3].set_index('year').sort_index()
    rows = []
    for y, r in c.iterrows():
        if pd.notna(r.get('newrel_coh')):
            n = dict(cohort=r.newrel_coh, success=r.newrel_succ, failed=r.newrel_fail, died=r.newrel_died, lost=r.newrel_lost, cohort_type='new and relapse')
        elif pd.notna(r.get('new_sp_coh')):
            n = dict(cohort=r.new_sp_coh, success=r.new_sp_cur + r.new_sp_cmplt, failed=r.new_sp_fail, died=r.new_sp_died, lost=r.new_sp_def,
                     cohort_type='new smear-positive')
        else:
            continue
        ev = n['success'] + n['failed'] + n['died'] + n['lost']
        rows.append(dict(year=int(y), **n, evaluated=ev, share_evaluated=ev/n['cohort'], **{f'{k}_share': n[k]/ev for k in ['success', 'failed', 'died', 'lost']}))
    return pd.DataFrame(rows).set_index('year')


def knots(df):
    """ d (died) and ℓ (lost) at each knot: the mean of the yearly shares of evaluated cases over the knot's window; success and failed too, for reference """
    rows = []
    for k, (a, b) in windows.items():
        w = df.loc[a:b]
        rows.append(dict(knot=k, years=f'{a}–{b}', **{s: float(w[f'{s}_share'].mean()) for s in ['success', 'failed', 'died', 'lost']},
                         share_evaluated=float(w.share_evaluated.mean()), cohort_types=', '.join(sorted(set(w.cohort_type)))))
    return pd.DataFrame(rows).set_index('knot')


def schedule(kn, key):
    """ A share over calendar time: linear between the knots, constant outside them (as tb_LAI_TPT's coverage and efficacy) """
    x, y = np.array(kn.index, dtype=float), kn[key].values
    return lambda t: np.interp(t, x, y)


def plot(df, kn):
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
    ax = axes[0]
    for s, c in [('died', '#d62728'), ('lost', '#1f77b4'), ('failed', '#7f7f7f')]:
        ax.plot(df.index, 100*df[f'{s}_share'], 'o-', ms=3, color=c, label=f'{s} (share of evaluated)')
    ax.axvline(2011.5, color='k', lw=0.6, ls=':')
    ax.text(2011.7, ax.get_ylim()[1]*0.95, 'new smear-positive | new and relapse', fontsize=7, va='top')
    ax.set_ylabel('% of evaluated cases'); ax.set_title('India: deaths, losses and failures on treatment (WHO)', fontsize=10)
    ax.legend(fontsize=7, frameon=False)
    ax = axes[1]
    ax.plot(df.index, 100*df.share_evaluated, 'o-', ms=3, color='k', label='Evaluated, % of cohort')
    ax.plot(df.index, 100*df.success_share, 'o-', ms=3, color='#2ca02c', label='Success, % of evaluated')
    ax.set_ylim(0, 105); ax.set_title('Evaluated share of each cohort, and treatment success', fontsize=10); ax.legend(fontsize=7, frameon=False)
    for a in axes:
        a.set_xlabel('Cohort year')
    fig.tight_layout()
    fig.savefig(here / 'outcomes.png', dpi=120)
    plt.close(fig)
    return


def main():
    if not (datadir / 'TB_outcomes.csv').exists():
        download()
    df = yearly_outcomes()
    kn = knots(df)
    df.to_csv(here / 'outcomes.csv')
    plot(df, kn)
    with pd.option_context('display.width', 200, 'display.max_columns', 20):
        text = ('India, WHO treatment outcomes by cohort year (shares of evaluated cases; not evaluated left out of the denominator):\n' +
                df[['cohort_type', 'cohort', 'share_evaluated', 'success_share', 'failed_share', 'died_share', 'lost_share']].round(4).to_string() +
                '\n\n(Smear-positive cohorts only to 2011; the schedules used are from all new cases: fitted_schedule(), outcome_schedule.csv)')
    sc.savetext(here / 'outcomes.log', text)
    print(text)
    return df, kn


if __name__ == '__main__' and not (len(sys.argv) > 1 and sys.argv[1] in ['options', 'fitted']):
    if len(sys.argv) > 1 and sys.argv[1] == 'download':
        download()
    main()


#%% Options for the outcome schedules: every new case (smear-negative and extrapulmonary too), yearly or smoothed

def all_new_cases(iso3='IND', ratio_years=(2004, 2006)):
    """
    Died, lost and failed as shares of evaluated cases, for every new case each year: new smear-positive plus new smear-negative/extrapulmonary/
    smear unknown cohorts (new_snep_*) for 2004–2011, new and relapse cohorts from 2012. Before 2004 only smear-positive cohorts were reported;
    they are scaled by the ratio of all new cases' share to smear-positive cases' over `ratio_years` (adjusted=True)
    """
    d = pd.read_csv(datadir / 'TB_outcomes.csv')
    c = d[d.iso3 == iso3].set_index('year').sort_index()
    rows = []
    for y, r in c.iterrows():
        sp = None
        if pd.notna(r.get('new_sp_coh')):
            sp = dict(ev=r.new_sp_cur + r.new_sp_cmplt + r.new_sp_died + r.new_sp_fail + r.new_sp_def, died=r.new_sp_died, lost=r.new_sp_def, failed=r.new_sp_fail)
        if pd.notna(r.get('newrel_coh')):
            n = dict(ev=r.newrel_succ + r.newrel_fail + r.newrel_died + r.newrel_lost, died=r.newrel_died, lost=r.newrel_lost, failed=r.newrel_fail)
            kind = 'new and relapse'
        elif pd.notna(r.get('new_snep_coh')) and sp is not None:
            sn = dict(ev=r.new_snep_cmplt + r.new_snep_died + r.new_snep_fail + r.new_snep_def, died=r.new_snep_died, lost=r.new_snep_def, failed=r.new_snep_fail)
            n = {k: sp[k] + sn[k] for k in sp}
            kind = 'all new'
        elif sp is not None:
            n, kind = sp, 'smear-positive only'
        else:
            continue
        rows.append(dict(year=int(y), kind=kind, **{k: n[k]/n['ev'] for k in ['died', 'lost', 'failed']},
                         **({f'sp_{k}': sp[k]/sp['ev'] for k in ['died', 'lost', 'failed']} if sp is not None else {})))
    df = pd.DataFrame(rows).set_index('year')
    a, b = ratio_years
    for k in ['died', 'lost', 'failed']:
        ratio = (df.loc[a:b, k]/df.loc[a:b, f'sp_{k}']).mean()
        df[f'{k}_adj'] = np.where(df.kind == 'smear-positive only', df[k]*ratio, df[k])
        df.attrs[f'ratio_{k}'] = float(ratio)
    return df


def schedule_options(df, first=1999, hold=(1999, 2001), last_year=2023, lam=3.0):
    """ The died and lost schedules by year, 1995 on: 'knots' (the three knots so far), 'yearly' (the adjusted yearly values, linear between),
    'smooth' (held at its 1999 value before 1999; a smoothing spline through the adjusted values, 1999–2023, penalty lam: GCV's choice barely smooths, and lam = 3 removes year-to-year noise in deaths but keeps the fall in losses); each held at its 1999–2001 mean before 1999 and at its
    last value after 2023 """
    from scipy.interpolate import make_smoothing_spline
    years = np.arange(1995, 2041)
    out = {}
    kn = knots(yearly_outcomes())
    for k, key in [('died', 'died'), ('lost', 'lost')]:
        s = df.loc[first:last_year, f'{k}_adj']
        x, v = s.index.values.astype(float), s.values
        spl = make_smoothing_spline(x, v, lam=lam)
        smooth = spl(np.clip(years, first, last_year).astype(float))
        smooth[years < first] = spl(float(first)) # Held at its own first value, so the curve is continuous
        yearly = np.interp(years, x, v)
        yearly[years < first] = s.loc[hold[0]:hold[1]].mean()
        out[k] = pd.DataFrame(dict(knots=np.interp(years, kn.index.values.astype(float), kn[key].values), yearly=yearly, smooth=smooth), index=years)
    return out


def plot_options():
    df = all_new_cases()
    opts = schedule_options(df)
    raw = yearly_outcomes()
    fig, axes = plt.subplots(1, 3, figsize=(19, 5))
    for ax, k in zip(axes[:2], ['died', 'lost']):
        ax.plot(raw.index, 100*raw[f'{k}_share'], 's', ms=3, color='0.6', label='Smear-positive only to 2011 (used so far)')
        sn = df[df.kind == 'all new']
        ax.plot(sn.index, 100*sn[k], 'o', ms=4, color='#1f77b4', label='All new cases, 2004–2011 (smear-pos. + smear-neg./EP)')
        nr = df[df.kind == 'new and relapse']
        ax.plot(nr.index, 100*nr[k], 'o', ms=4, color='#2ca02c', label='New and relapse, 2012 on')
        sp = df[df.kind == 'smear-positive only']
        ax.plot(sp.index, 100*sp[f'{k}_adj'], 'o', ms=4, mfc='white', color='#1f77b4', label=f'1994–2003 scaled to all new cases (× {df.attrs[f"ratio_{k}"]:.2f}, 2004–06 ratio)')
        o = opts[k]
        ax.plot(o.index, 100*o.knots, '--', color='0.4', lw=1.2, label='Knots so far (1995, 2010, 2024)')
        ax.plot(o.index, 100*o.yearly, '-', color='#ff7f0e', lw=1.2, label='Yearly values (adjusted), linear')
        ax.plot(o.index, 100*o.smooth, '-', color='#d62728', lw=2, label='Smoothing spline, 1999–2023 (penalty 3)')
        ax.axvspan(1994.5, 1998.5, color='0.92', zorder=0)
        ax.text(1995, ax.get_ylim()[1]*0.97, 'held before 1999\n(cohorts mostly\nnot evaluated)', fontsize=7, va='top')
        ax.set_title(f'{"Died" if k == "died" else "Lost to follow-up"}, % of evaluated new cases (India, WHO)', fontsize=10)
        ax.set_xlabel('Cohort year'); ax.set_xlim(1994, 2041)
    axes[0].legend(fontsize=7, frameon=True, loc='lower right')
    ax = axes[2]
    ax.plot(raw.index, 100*raw.share_evaluated, 'o-', ms=3, color='k', label='Evaluated, % of cohort (smear-pos. to 2011)')
    ax.plot(raw.index, 100*raw.success_share, 'o-', ms=3, color='#2ca02c', label='Success, % of evaluated')
    ax.plot(df.index, 100*df.failed_adj, 'o-', ms=3, color='#9467bd', label='Failed, % of evaluated (all new, adjusted)')
    ax.set_ylim(0, 105); ax.legend(fontsize=7, frameon=False); ax.set_title('Evaluated share, success and failure', fontsize=10); ax.set_xlabel('Cohort year')
    fig.tight_layout(); fig.savefig(here / 'outcome_options.png', dpi=120); plt.close(fig)
    rows = pd.concat({k: v for k, v in opts.items()}, axis=1)
    rows.to_csv(here / 'outcome_options.csv')
    df.to_csv(here / 'outcomes_all_new.csv')
    return df, opts


if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'options':
    plot_options()


#%% The schedules used: deaths constant, losses linear between knots fitted to the data

lost_knots = [1999, 2003, 2014, 2023] # Where the decline in losses changes pace: the fast fall after DOTS's rollout, the plateau, the later fall


def fitted_schedule(df=None, first=1999, last_year=2023):
    """
    The died and lost shares used by the model (user, 2026-11-23): died constant, the mean of all new cases' adjusted shares 1999–2023
    (its trend is noise-dominated); lost linear between `lost_knots`, the continuous piecewise-linear curve fitted by least squares to the adjusted
    shares 1999–2023, constant before the first knot and after the last. Returns the knots {year: value} for each, and the fit's residual SD
    """
    df = all_new_cases() if df is None else df
    s = df.loc[first:last_year]
    x = s.index.values.astype(float)
    kn = np.array(lost_knots, dtype=float)
    basis = np.column_stack([np.interp(x, kn, np.eye(len(kn))[j]) for j in range(len(kn))]) # Hat functions: the curve's value at each knot
    coef, *_ = np.linalg.lstsq(basis, s['lost_adj'].values, rcond=None)
    lost = {int(k): float(v) for k, v in zip(kn, coef)}
    died_mean = float(s['died_adj'].mean())
    return sc.objdict(died={first: died_mean}, lost=lost, died_sd=float(s['died_adj'].std()), lost_resid_sd=float(np.std(s['lost_adj'].values - basis @ coef)))


def plot_fitted():
    df = all_new_cases()
    fs = fitted_schedule(df)
    years = np.arange(1995, 2041)
    raw = yearly_outcomes()
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
    for ax, k in zip(axes, ['died', 'lost']):
        ax.plot(raw.index, 100*raw[f'{k}_share'], 's', ms=3, color='0.6', label='Smear-positive only to 2011 (used in the first runs)')
        for kind, c, mfc in [('all new', '#1f77b4', None), ('new and relapse', '#2ca02c', None), ('smear-positive only', '#1f77b4', 'white')]:
            d = df[df.kind == kind]
            lab = dict(**{'all new': 'All new cases, 2004–2011', 'new and relapse': 'New and relapse, 2012 on', 'smear-positive only': f'1994–2003, scaled to all new cases (× {df.attrs[f"ratio_{k}"]:.2f})'})[kind]
            ax.plot(d.index, 100*d[f'{k}_adj'], 'o', ms=4, color=c, mfc=mfc or c, label=lab)
        kn = fs[k]
        xs = np.array(sorted(kn), dtype=float)
        ax.plot(years, 100*np.interp(years, xs, [kn[int(x)] for x in xs]), '-', color='#d62728', lw=2, label='Used: ' + ('constant (mean 1999–2023)' if k == 'died' else 'linear between fitted knots'))
        if k == 'lost':
            ax.plot(xs, [100*kn[int(x)] for x in xs], 'D', color='#d62728', ms=6)
        ax.axvspan(1994.5, 1998.5, color='0.92', zorder=0)
        ax.set_title(f'{"Died" if k == "died" else "Lost to follow-up"}, % of evaluated new cases (India, WHO)', fontsize=10)
        ax.set_xlabel('Cohort year'); ax.set_xlim(1994, 2041)
        ax.legend(fontsize=7, frameon=True, loc='lower left' if k == 'died' else 'upper right')
    fig.tight_layout(); fig.savefig(here / 'outcome_schedule.png', dpi=120); plt.close(fig)
    rows = [dict(year=y, died=float(np.interp(y, sorted(fs.died), [fs.died[k] for k in sorted(fs.died)])),
                 lost=float(np.interp(y, sorted(fs.lost), [fs.lost[k] for k in sorted(fs.lost)]))) for y in years]
    pd.DataFrame(rows).to_csv(here / 'outcome_schedule.csv', index=False)
    sc.savetext(here / 'outcome_schedule.log', f'Died: constant {fs.died} (SD of yearly values {fs.died_sd:.4f}); lost knots {fs.lost} (residual SD {fs.lost_resid_sd:.4f})')
    print(f'Died: constant {fs.died} (SD of yearly values {fs.died_sd:.4f}); lost knots ' + ', '.join(f'{k}: {v:.4f}' for k, v in fs.lost.items()) + f' (residual SD {fs.lost_resid_sd:.4f})')
    return fs


if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'fitted':
    plot_fitted()
