"""
Put ODE and ABM outputs on a common footing and compare them.
"""

import numpy as np
import sciris as sc
import matplotlib.pyplot as plt

__all__ = ['compare_states', 'ode_outputs', 'abm_outputs', 'annualize', 'summarize', 'plot_comparison']

# States compared, using the ABM's names. The ODE's CLEARED/RECOVERED/TREATED are summed into CLEARED
compare_states = ['SUSCEPTIBLE', 'INFECTION', 'CLEARED', 'NON_INFECTIOUS', 'ASYMPTOMATIC', 'SYMPTOMATIC']


def ode_outputs(res):
    """Proportions in each compared state plus prev_active, frac_asymp, incidence_py and tb_deaths_py, from TBODE results."""
    out = sc.objdict(time=res.time, N=res.N)
    for s in compare_states:
        num = res.CLEARED + res.RECOVERED + res.TREATED if s == 'CLEARED' else res[s]
        out[s] = num / res.N
    for k in ['prev_active', 'frac_asymp', 'incidence_py', 'tb_deaths_py']:
        out[k] = res[k]
    out.tb_deaths_ppy = res.tb_deaths_py
    return out


def abm_outputs(sim):
    """
    The same outputs as ode_outputs, from a run tbsim.Sim with a StateCounts analyzer. Per-step event rates
    (incidence_py, tb_deaths_py) are noisy; use annualize() before comparing them.
    """
    r = sim.results.statecounts
    n = np.asarray(r.n_alive, dtype=float)
    dt = sim.t.dt_year
    out = sc.objdict(time=np.asarray(sim.t.yearvec), N=n, dt=dt)
    with np.errstate(divide='ignore', invalid='ignore'):
        for s in compare_states:
            out[s] = np.asarray(r[s]) / n
        asy, sym = np.asarray(r.ASYMPTOMATIC, float), np.asarray(r.SYMPTOMATIC, float)
        out.prev_active = (asy + sym) / n
        out.frac_asymp = asy / (asy + sym)
        out.incidence_py = np.asarray(r.new_incident) / (n * dt)
        out.tb_deaths_py = np.asarray(sim.get_tb().results.new_deaths) / (n * dt)
    out.tb_deaths_ppy = out.tb_deaths_py
    return out


def annualize(out, keys=('incidence_py', 'tb_deaths_py')):
    """Return a copy of ``out`` with ``keys`` replaced by calendar-year averages (assigned to every step in that year)."""
    out = sc.dcp(out)
    yr = np.floor(out.time + 1e-9)
    for k in keys:
        v = np.asarray(out[k], dtype=float)
        avg = np.full_like(v, np.nan)
        for y in np.unique(yr):
            m = yr == y
            avg[m] = np.nanmean(v[m])
        out[k] = avg
    out.tb_deaths_ppy = out.tb_deaths_py
    return out


def summarize(outs, quantiles=(0.1, 0.5, 0.9)):
    """Quantiles across a list of abm_outputs (one per seed). Returns {key: array of shape (len(quantiles), n_t)}."""
    keys = [k for k in outs[0] if k not in ('time', 'dt')]
    summ = sc.objdict(time=outs[0].time)
    for k in keys:
        summ[k] = np.nanquantile(np.vstack([o[k] for o in outs]), quantiles, axis=0)
    return summ


def plot_comparison(odes, abms, keys=None, xlim=None, title=None, scale=None):
    """
    Plot ODE lines against ABM median and 10-90% bands.

    Args:
        odes (dict): {label: ode_outputs}
        abms (dict): {label: summarize output}
        keys (list): outputs to plot (default: compared states plus prev_active and tb_deaths_ppy)
        scale (dict): {key: multiplier} for display (e.g. {'prev_active': 1e5})
    """
    keys = keys or compare_states + ['prev_active', 'tb_deaths_ppy']
    scale = scale or {}
    fig, axs = plt.subplots(2, int(np.ceil(len(keys) / 2)), figsize=(16, 7), squeeze=False)
    for ax, k in zip(axs.flat, keys):
        m = scale.get(k, 1)
        for i, (lab, a) in enumerate(abms.items()):
            col = f'C{i}'
            ax.fill_between(a.time, m * a[k][0], m * a[k][-1], color=col, alpha=0.25, lw=0)
            ax.plot(a.time, m * a[k][1], color=col, lw=1.5, label=f'ABM: {lab}')
        for j, (lab, o) in enumerate(odes.items()):
            ax.plot(o.time, m * o[k], color='k', ls=['--', ':', '-.'][j % 3], lw=1.5, label=f'ODE: {lab}')
        ax.set_title(k + (f' (×{m:g})' if m != 1 else ''))
        ax.set_ylim(bottom=0)
        if xlim:
            ax.set_xlim(xlim)
    for ax in axs.flat[len(keys):]:
        ax.set_visible(False)
    axs.flat[0].legend(fontsize=8, frameon=False)
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    return fig
