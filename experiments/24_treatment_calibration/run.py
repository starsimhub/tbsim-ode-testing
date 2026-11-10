"""
Experiment 24, v4 environment: calibrate the ODE with treatment, time-varying sym_dead and over-detection to tb_LAI_TPT's applicable
targets, with IMIS. See README.md. Run from the repo root, in order:

    PYTHONPATH=v4 v4/.venv/bin/python experiments/24_treatment_calibration/run.py coverage    # prior predictive check
    PYTHONPATH=v4 v4/.venv/bin/python experiments/24_treatment_calibration/run.py reid        # re-identification from a known point
    PYTHONPATH=v4 v4/.venv/bin/python experiments/24_treatment_calibration/run.py real        # the real targets
    PYTHONPATH=v4 v4/.venv/bin/python experiments/24_treatment_calibration/run.py real_seed2  # the same with another seed (convergence)
    PYTHONPATH=v4 v4/.venv/bin/python experiments/24_treatment_calibration/run.py plots       # summaries and figures from the saved runs
"""
import os
os.environ.setdefault('NUMBA_NUM_THREADS', '1') # One thread per process (Starsim's kernels otherwise use every core); before importing starsim
import sys
import numpy as np
import pandas as pd
import sciris as sc
import matplotlib.pyplot as plt

outdir = sc.thispath(__file__)
sys.path.insert(0, str(outdir))
import calib as cb # noqa: E402
import tbode4 # noqa: E402
from tbode4 import fastode as fo # noqa: E402

fn = sc.objdict(coverage=outdir / 'coverage.obj', reid=outdir / 'reid.obj', real=outdir / 'real.obj', real_seed2=outdir / 'real_seed2.obj')

# The re-identification truth: exp 21's set 1 natural history, with tb_LAI_TPT's best-fit treatment, sym_dead trend and over-detection
lai = dict(sym_dead_2024=0.14261830792417335, tx_covg_2024=0.6772271561810697, tx_covg_2010=0.5478919423969836, efficacy_1995=0.4525494678080442,
           efficacy_2024=0.5300106194827879, p_relapse=0.03394965858858717, over_detection=0.18310258381507125)
truth = dict(beta=22.10, trans_asymp=0.7488, pi=0.3518, inf_cle=0.6638, inf_non=0.3153, inf_asy=0.3735, k_latent=7.517, phi=0.01194, non_rec=0.2207,
             non_asy=0.2776, asy_non=1.612, asy_sym=0.7560, sym_asy=0.7593, sym_dead_1995=0.3959,
             sym_dead_frac=(0.3959 - lai['sym_dead_2024'])/(0.3959 - cb.sym_dead_floor), tx_covg_2024=lai['tx_covg_2024'],
             tx_covg_2010_frac=lai['tx_covg_2010']/lai['tx_covg_2024'], efficacy_1995=lai['efficacy_1995'],
             efficacy_frac=(lai['efficacy_2024'] - lai['efficacy_1995'])/(cb.efficacy_ceiling - lai['efficacy_1995']), p_relapse=lai['p_relapse'],
             over_detection=lai['over_detection'])
labels = dict(prevalence='Prevalence 2020 (per 100k)', prop_asymptomatic='Share asymptomatic', prev_to_inc='Prevalence ÷ incidence (y)',
              mort_to_inc='Mortality ÷ incidence', prev_to_notif='Prevalence ÷ notifications (y)', inc_ratio='Incidence 2024 ÷ 2015',
              mort_ratio='Mortality 2024 ÷ 2015', prev_2040_ratio='Prevalence 2040 ÷ 2020', ProgressionY1='ProgressionY1', ProgressionRatio='ProgressionRatio')
scale = dict(prevalence=1e5)
# Data over time, per 100,000 (tb_LAI_TPT's components): WHO incidence and mortality, notifications, the survey's prevalence
data = sc.objdict(incidence={2015: 237.197, 2020: 197.482, 2024: 187.31}, mortality={2015: 32.03, 2020: 22.918, 2024: 21.073},
                  notifications={2015: 126.222, 2020: 116.666, 2024: 172.908}, prevalence={2020: 310.296})


def target_interval(t):
    """ The target's mean and 95% interval (Jeffreys for a Binomial, ±1.96σ in the log for a log-Normal) """
    if t['kind'] == 'binomial':
        from scipy.stats import beta
        return t['x']/t['n'], *beta.ppf([0.025, 0.975], t['x'] + 0.5, t['n'] - t['x'] + 0.5)
    return t['mean'], t['mean']*np.exp(-1.96*t['sigma']), t['mean']*np.exp(1.96*t['sigma'])


def synthetic_targets():
    """ Targets at the truth's own outputs, with the real targets' sample sizes and sigmas """
    out = cb.simulate(np.array([truth[k] for k in cb.names]))
    tg = sc.dcp(cb.targets)
    for k, t in tg.items():
        if t['kind'] == 'binomial':
            t['x'] = int(round(out[k]*t['n']))
        else:
            t['mean'] = float(out[k])
    return tg, out


def coverage(n=5000, seed=2):
    """ Prior predictive check: the share of prior draws with an endemic equilibrium, and each output's spread against its target """
    Z = np.random.default_rng(seed).random((n, cb.d))
    outs, ll = cb.evaluate(Z)
    df = pd.DataFrame([dict(o) for o in outs])
    sc.save(fn.coverage, sc.objdict(Z=Z, theta=cb.from_unit(Z), outs=df, ll=ll))
    print(f'{n} prior draws: no endemic equilibrium in {df.prevalence.isna().mean():.1%}')
    rows = []
    for k, t in cb.targets.items():
        m, lo, hi = target_interval(t)
        v = df[k].dropna().values
        rows.append(dict(target=k, target_mean=m, target_lo=lo, target_hi=hi, prior_p05=np.percentile(v, 5), prior_p50=np.percentile(v, 50),
                         prior_p95=np.percentile(v, 95), share_in_interval=np.mean((v >= lo) & (v <= hi)), share_below=np.mean(v < lo)))
    cov = pd.DataFrame(rows)
    cov.to_csv(outdir / 'coverage.csv', index=False)
    every = np.ones(len(df), bool)
    for k, t in cb.targets.items():
        _, lo, hi = target_interval(t)
        every &= (df[k] >= lo).values & (df[k] <= hi).values
    with pd.option_context('display.width', 220):
        print(cov.round(4).to_string(index=False))
    print(f'Prior draws inside every target interval at once: {every.sum()} of {n}')
    return


def run_imis(part, tgts, **kwargs):
    res = cb.imis(tgts, **kwargs)
    res.targets = tgts
    sc.save(fn[part], res)
    tbode4.save_versions(outdir, 'v4')
    return res


def weighted_quantile(x, w, q):
    i = np.argsort(x)
    x, w = x[i], w[i]
    cw = (np.cumsum(w) - 0.5*w)/w.sum()
    return np.interp(q, cw, x)


def summarise(part):
    """ Posterior table (weighted mean and 2.5/50/97.5% per parameter, against the prior) and the fit to each target """
    res = sc.load(fn[part])
    w = res.weights
    rows = []
    for j, k in enumerate(cb.names):
        x = res.theta[:, j]
        lo, hi, _ = cb.priors[k]
        q = weighted_quantile(x, w, [0.025, 0.5, 0.975])
        row = dict(parameter=k, prior_lo=lo, prior_hi=hi, mean=float(w @ x), q025=q[0], q50=q[1], q975=q[2], ci_width_share_of_prior=(q[2] - q[0])/(hi - lo))
        if part == 'reid':
            row.update(truth=truth[k], truth_in_95=bool(q[0] <= truth[k] <= q[2]))
        rows.append(row)
    pars = pd.DataFrame(rows)
    pars.to_csv(outdir / f'{part}_posterior.csv', index=False)
    df = pd.DataFrame([dict(o) for o in res.outs])
    rows = []
    for k, t in res.targets.items():
        m, lo, hi = target_interval(t)
        v = df[k].values
        ok = np.isfinite(v)
        q = weighted_quantile(v[ok], w[ok], [0.025, 0.5, 0.975])
        rows.append(dict(target=k, target_mean=m, target_lo=lo, target_hi=hi, post_q025=q[0], post_q50=q[1], post_q975=q[2]))
    for k in cb.extra:
        v = df[k].values
        ok = np.isfinite(v)
        q = weighted_quantile(v[ok], w[ok], [0.025, 0.5, 0.975])*1e5
        rows.append(dict(target=k + ' (per 100k)', post_q025=q[0], post_q50=q[1], post_q975=q[2]))
    fit = pd.DataFrame(rows)
    fit.to_csv(outdir / f'{part}_fit.csv', index=False)
    last = res.log[-1]
    print(f'{part}: {last["n"]:,} model runs, {len(res.comps)} IMIS components, ESS {last["ess"]:.0f}, expected unique {last["expected_unique"]:.0f} of {res.j}')
    with pd.option_context('display.width', 220):
        print(pars.round(4).to_string(index=False))
        print(fit.round(5).to_string(index=False))
    return pars, fit


def posterior_draws(part, n=200, seed=24):
    """ A weighted resample of parameter sets (with replacement), and their trajectories 1995–2040 per 100,000 """
    res = sc.load(fn[part])
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(res.weights), size=n, replace=True, p=res.weights)
    years = np.arange(1995, 2041, 1.0)
    t_eval = np.concatenate([[cb.t0], years])
    rows = []
    for i in idx:
        r = cb.trajectories(res.theta[i], t_eval=t_eval)
        prev = (r.X[:, fo.A] + r.X[:, fo.Y])*1e5
        od = r.trend.over_detection
        for j, y in enumerate(years[:-1]):
            rows.append(dict(index=int(i), year=int(y), prevalence=prev[j + 1], incidence=(r.cum_inc[j + 2] - r.cum_inc[j + 1])*1e5,
                             notifications=(1 + od)*(r.cum_notif[j + 2] - r.cum_notif[j + 1])*1e5, mortality=(r.cum_deaths[j + 2] - r.cum_deaths[j + 1])*1e5))
    df = pd.DataFrame(rows)
    df.to_csv(outdir / f'{part}_trajectories.csv', index=False)
    return df


def plot(part):
    """ Marginals against the prior, the posterior predictive of each target, and trajectories against the data """
    res = sc.load(fn[part])
    w = res.weights
    fig, axes = plt.subplots(4, 6, figsize=(18, 11))
    for ax, (j, k) in zip(axes.flat, enumerate(cb.names)):
        lo, hi, lg = cb.priors[k]
        bins = np.geomspace(lo, hi, 40) if lg else np.linspace(lo, hi, 40)
        ax.hist(res.theta[:, j], bins=bins, weights=w, density=True, color='#3b6ea8', alpha=0.8)
        ax.axhline(1/(np.log(hi) - np.log(lo))/np.mean(bins) if lg else 1/(hi - lo), color='0.5', lw=0.8, ls='--')
        if lg:
            ax.set_xscale('log')
        if part == 'reid':
            ax.axvline(truth[k], color='#c4620b', lw=1.5)
        ax.set_title(k, fontsize=9)
        ax.set_yticks([])
    for ax in axes.flat[len(cb.names):]:
        ax.axis('off')
    fig.suptitle(f'Exp 24, {part}: posterior marginals (blue; dashed: prior, approximate for log-uniform{"; orange: truth" if part == "reid" else ""})', fontsize=10)
    fig.tight_layout()
    fig.savefig(outdir / f'{part}_marginals.png', dpi=130)
    plt.close(fig)

    df = pd.DataFrame([dict(o) for o in res.outs])
    fig, axes = plt.subplots(2, 5, figsize=(17, 6))
    for ax, (k, t) in zip(axes.flat, res.targets.items()):
        m, lo, hi = target_interval(t)
        s = scale.get(k, 1)
        v = df[k].values
        ok = np.isfinite(v)
        q = weighted_quantile(v[ok], w[ok], [0.001, 0.999])
        edges = np.linspace(min(q[0], lo), max(q[1], hi), 41)*s
        ax.hist(v[ok]*s, bins=edges, weights=w[ok], density=True, color='#3b6ea8', alpha=0.8)
        ax.set_xlim(edges[0], edges[-1])
        ax.axvspan(lo*s, hi*s, color='0.85', zorder=0)
        ax.axvline(m*s, color='0.2', lw=1)
        ax.set_title(labels[k], fontsize=9)
        ax.set_yticks([])
    fig.suptitle(f'Exp 24, {part}: posterior predictive of each target (blue) against the target mean and 95% interval (grey)', fontsize=10)
    fig.tight_layout()
    fig.savefig(outdir / f'{part}_fit.png', dpi=130)
    plt.close(fig)

    tr = posterior_draws(part)
    fig, axes = plt.subplots(1, 4, figsize=(17, 3.8))
    for ax, key in zip(axes, ['prevalence', 'incidence', 'notifications', 'mortality']):
        g = tr.groupby('year')[key]
        ax.fill_between(g.quantile(0.025).index, g.quantile(0.025), g.quantile(0.975), color='#3b6ea8', alpha=0.25, lw=0, label='Posterior 95%')
        ax.plot(g.median().index, g.median(), color='#1f4e8c', lw=1.5, label='Posterior median')
        if key in data:
            yrs = list(data[key])
            ax.plot(yrs, [data[key][y] for y in yrs], 'o', color='k', ms=4, label='Data (targets use ratios)')
        ax.set_title({'prevalence': 'Prevalence (1 Jan)', 'incidence': 'Incidence', 'notifications': 'Reported notifications', 'mortality': 'TB mortality'}[key] + ' per 100,000', fontsize=10)
        ax.set_xlabel('Year')
        ax.set_ylim(bottom=0)
    axes[0].legend(fontsize=8, frameon=False)
    fig.suptitle(f'Exp 24, {part}: trajectories of 200 posterior draws (all ages; the prevalence survey is 15+)', fontsize=10)
    fig.tight_layout()
    fig.savefig(outdir / f'{part}_trajectories.png', dpi=130)
    plt.close(fig)
    print(f'Saved {part}_marginals.png, {part}_fit.png, {part}_trajectories.png')
    return


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'coverage'
    if which == 'coverage':
        coverage()
    elif which == 'reid':
        tg, out = synthetic_targets()
        print('Synthetic targets (the truth\'s outputs):', {k: round(float(v), 5) for k, v in out.items()})
        run_imis('reid', tg)
        summarise('reid')
        plot('reid')
    elif which == 'real':
        run_imis('real', cb.targets)
        summarise('real')
        plot('real')
    elif which == 'real_seed2':
        run_imis('real_seed2', cb.targets, seed=2)
        summarise('real_seed2')
    elif which == 'plots':
        for part in ['reid', 'real']:
            if fn[part].exists():
                summarise(part)
                plot(part)
