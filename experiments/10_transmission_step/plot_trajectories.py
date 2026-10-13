"""
Experiment 10: prevalence and incidence trajectories per arm (mean across reps ± 2 SE) against the ODE equilibrium, from outputs.obj.
Saves trajectories.png and the plotted series (trajectories.csv).
"""

import os
import numpy as np
import pandas as pd
import sciris as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import tbode
from step_map import settings

outdir = os.path.dirname(os.path.abspath(__file__))
smooth_years = 5  # Moving average for incidence (annual values are noisy)

panels = {  # Column: (title, setting, arms)
    'A': ('Set A: tbsim defaults + k (near threshold), monthly dt', 'A', ['1: A, network', '2: A, mass action exp', '3: A, mass action linear']),
    'B': ('Set B: tb_LAI_TPT natural history + k, monthly dt', 'B', ['4: B, network', '5: B, mass action linear']),
    'k0': ('k = 0, tbsim defaults (near threshold)', 'k0', ['6: k=0 monthly start, network', '7: k=0 weekly within, mass action exp']),
}
styles = {
    '1: A, network': ('Network (= 1 − e^(−x))', 'C3', '-'),
    '2: A, mass action exp': ('Mass action, 1 − e^(−x)', 'C1', '--'),
    '3: A, mass action linear': ('Mass action, linear x', 'C0', '-'),
    '4: B, network': ('Network (= 1 − e^(−x))', 'C3', '-'),
    '5: B, mass action linear': ('Mass action, linear x', 'C0', '-'),
    '6: k=0 monthly start, network': ("Network, monthly, 'start' (exp 05 arm)", 'C3', '-'),
    '7: k=0 weekly within, mass action exp': ("Mass action 1 − e^(−x), weekly, 'within'", 'C1', '--'),
}
outputs = {'prev_active': 'Prevalence (ASY + SYM) per 100k', 'incidence_py': 'Incidence per 100k per year'}


def yearly(o, k):
    """Calendar-year values of one rep's output (annualize() already assigns each year's mean to every step in it)."""
    yr = np.floor(o.time + 1e-9)
    years = np.unique(yr)
    return years, np.array([np.nanmean(np.asarray(o[k])[yr == y]) for y in years])


if __name__ == '__main__':
    d = sc.load(os.path.join(outdir, 'outputs.obj'))
    abms = d['abms']
    sets = settings()
    rows = []
    fig, axs = plt.subplots(2, 3, figsize=(17, 9), sharex=True)
    for j, (title, setting, arm_list) in enumerate(panels.values()):
        ode = tbode.TBODE(**sets[setting])
        eq = ode.outputs(ode.equilibrium())
        for i, (k, ylabel) in enumerate(outputs.items()):
            ax = axs[i, j]
            ode_val = eq[k][0] * 1e5
            ax.axhline(ode_val, color='k', lw=1.5, ls=':', label=f'ODE equilibrium ({ode_val:.0f})')
            for arm in arm_list:
                oo = abms[arm]
                years = yearly(oo[0], k)[0]
                Y = np.vstack([yearly(o, k)[1] for o in oo]) * 1e5
                mean, se = Y.mean(axis=0), Y.std(axis=0, ddof=1) / np.sqrt(len(Y))
                if k == 'incidence_py':  # Smooth the mean and SE with a centred moving average
                    kern = np.ones(smooth_years) / smooth_years
                    mean, se = [np.convolve(v, kern, mode='same') / np.convolve(np.ones_like(v), kern, mode='same') for v in (mean, se)]
                last50 = np.mean([np.nanmean(o[k][o.time >= 1950]) for o in oo]) * 1e5
                lab, col, ls = styles[arm]
                ax.fill_between(years, mean - 2 * se, mean + 2 * se, color=col, alpha=0.15, lw=0)
                ax.plot(years, mean, color=col, ls=ls, lw=1.8, label=f'{lab}: {last50:.0f} ({last50 / ode_val:.2f}×, {len(oo)} reps)')
                rows += [dict(panel=setting, arm=arm, output=k, year=y, mean=m, se=s, ode=ode_val) for y, m, s in zip(years, mean, se)]
            ax.set_ylim(0, ode_val * 1.25)
            ax.set_ylabel(ylabel)
            if i == 0:
                ax.set_title(title, fontsize=10)
            if i == 1:
                ax.set_xlabel('Year')
            ax.legend(fontsize=7.5, frameon=False, loc='lower left')
    fig.suptitle('TBExact transmission options vs the ODE: mean across reps (band: ±2 SE), started at the ODE equilibrium.\n'
                 'Legend values: mean over 1950–2000 (as in the summary tables). Incidence: 5-year moving average.', fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, 'trajectories.png'), dpi=110)
    pd.DataFrame(rows).to_csv(os.path.join(outdir, 'trajectories.csv'), index=False)
    print('Saved trajectories.png and trajectories.csv')
