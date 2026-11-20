"""
Experiment 32, v4 environment: recalibrate the ODE with ages and India's demography (tbode4.fastode_age), with exp 24's IMIS protocol. See
README.md. Run from the repo root, in order (on zebra, with EXP32_NCPUS=160):

    PYTHONPATH=v4 v4/.venv/bin/python experiments/32_demography_calibration/run.py coverage    # prior predictive check
    PYTHONPATH=v4 v4/.venv/bin/python experiments/32_demography_calibration/run.py reid        # re-identification from a known point
    PYTHONPATH=v4 v4/.venv/bin/python experiments/32_demography_calibration/run.py real        # the real targets, IMIS seed 1
    PYTHONPATH=v4 v4/.venv/bin/python experiments/32_demography_calibration/run.py real_seed2  # IMIS seed 2 (convergence)
    PYTHONPATH=v4 v4/.venv/bin/python experiments/32_demography_calibration/run.py summary     # tables and figures from the saved runs

`all` runs them in order.
"""
import os
for _k in ['NUMBA_NUM_THREADS', 'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']:
    os.environ.setdefault(_k, '1')
import sys
import numpy as np
import pandas as pd
import sciris as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

outdir = sc.thispath(__file__)
sys.path.insert(0, str(outdir))
import calib as cb # noqa: E402
import tbode4 # noqa: E402

fn = sc.objdict({k: outdir / f'{k}.obj' for k in ['coverage', 'reid', 'real', 'real_seed2']})
exp24 = cb.root / 'experiments' / '24_treatment_calibration'

# The re-identification truth: exp 24's (exp 21's set 1 natural history, with tb_LAI_TPT's best-fit treatment, sym_dead trend and over-detection)
lai = dict(sym_dead_2024=0.14261830792417335, tx_covg_2024=0.6772271561810697, tx_covg_2010=0.5478919423969836, efficacy_1995=0.4525494678080442,
           efficacy_2024=0.5300106194827879, p_relapse=0.03394965858858717, over_detection=0.18310258381507125)
truth = dict(beta=22.10, trans_asymp=0.7488, pi=0.3518, inf_cle=0.6638, inf_non=0.3153, inf_asy=0.3735, k_latent=7.517, phi=0.01194, non_rec=0.2207,
             non_asy=0.2776, asy_non=1.612, asy_sym=0.7560, sym_asy=0.7593, sym_dead_1995=0.3959,
             sym_dead_frac=(0.3959 - lai['sym_dead_2024'])/(0.3959 - cb.cb24.sym_dead_floor),
             tx_covg_2024=lai['tx_covg_2024'], tx_covg_2010_frac=lai['tx_covg_2010']/lai['tx_covg_2024'], efficacy_1995=lai['efficacy_1995'],
             efficacy_frac=(lai['efficacy_2024'] - lai['efficacy_1995'])/(cb.cb24.efficacy_ceiling - lai['efficacy_1995']), p_relapse=lai['p_relapse'],
             over_detection=lai['over_detection'])
labels = dict(prevalence='Prevalence 15+, 2020 (per 100k)', prop_asymptomatic='Share asymptomatic, 15+', prev_to_inc='Prevalence 15+ ÷ incidence (y)',
              mort_to_inc='Mortality ÷ incidence', prev_to_notif='Prevalence 15+ ÷ notifications (y)', inc_ratio='Incidence 2024 ÷ 2015',
              mort_ratio='Mortality 2024 ÷ 2015', prev_2040_ratio='Prevalence 15+, 2040 ÷ 2020', ProgressionY1='ProgressionY1', ProgressionRatio='ProgressionRatio')
scale = dict(prevalence=1e5)
# Data over time, per 100,000 (tb_LAI_TPT's components): WHO incidence and mortality, notifications, the survey's 15+ prevalence
data = sc.objdict(incidence={2015: 237.197, 2020: 197.482, 2024: 187.31}, mortality={2015: 32.03, 2020: 22.918, 2024: 21.073},
                  notifications={2015: 126.222, 2020: 116.666, 2024: 172.908}, prev15={2020: 310.296})
colors = sc.objdict(real='#1f4e8c', real_seed2='#c4620b', exp24='0.45', exp24_seed2='0.7')


def target_interval(t):
    """ The target's mean and 95% interval (Jeffreys for a Binomial, ±1.96σ in the log for a log-Normal) """
    if t['kind'] == 'binomial':
        from scipy.stats import beta
        return t['x']/t['n'], *beta.ppf([0.025, 0.975], t['x'] + 0.5, t['n'] - t['x'] + 0.5)
    return t['mean'], t['mean']*np.exp(-1.96*t['sigma']), t['mean']*np.exp(1.96*t['sigma'])


def weighted_quantile(x, w, q):
    i = np.argsort(x)
    x, w = x[i], w[i]
    cw = (np.cumsum(w) - 0.5*w)/w.sum()
    return np.interp(q, cw, x)


def truth_theta():
    return np.array([truth[k] for k in cb.names])


def synthetic_targets():
    """ Targets at the truth's own outputs, with the real targets' sample sizes and sigmas """
    out = cb.simulate(truth_theta())
    tg = sc.dcp(cb.targets)
    for k, t in tg.items():
        if t['kind'] == 'binomial':
            t['x'] = int(round(out[k]*t['n']))
        else:
            t['mean'] = float(out[k])
    return tg, out


def outs_df(res):
    return pd.DataFrame(res.outs, columns=res.out_keys)


#%% Runs

def coverage(n=5000, seed=2):
    """ Prior predictive check: the share of prior draws with an endemic state, and each output's spread against its target """
    Z = np.random.default_rng(seed).random((n, cb.d))
    T = sc.timer()
    outs, ll = cb.evaluate(Z)
    secs = T.toc(output=True)
    sc.save(fn.coverage, sc.objdict(Z=Z, theta=cb.from_unit(Z), outs=outs, out_keys=cb.out_keys, ll=ll, secs=secs, ncpus=cb.ncpus))
    tbode4.save_versions(outdir, 'v4')
    df = pd.DataFrame(outs, columns=cb.out_keys)
    lines = [f'{n} prior draws in {secs:.0f} s on {cb.ncpus} workers ({secs*cb.ncpus/n:.2f} s per draw per worker); no endemic state in {df.prevalence.isna().mean():.1%}; '
             f'errors {int(df.error.sum())}']
    rows = []
    every = np.ones(len(df), bool)
    for k, t in cb.targets.items():
        m, lo, hi = target_interval(t)
        v = df[k].dropna().values
        rows.append(dict(target=k, target_mean=m, target_lo=lo, target_hi=hi, prior_p05=np.percentile(v, 5), prior_p50=np.percentile(v, 50),
                         prior_p95=np.percentile(v, 95), share_in_interval=np.mean((v >= lo) & (v <= hi)), share_below=np.mean(v < lo)))
        every &= (df[k] >= lo).values & (df[k] <= hi).values
    cov = pd.DataFrame(rows)
    cov.to_csv(outdir / 'coverage.csv', index=False)
    with pd.option_context('display.width', 220):
        lines.append(cov.round(4).to_string(index=False))
    lines.append(f'Prior draws inside every target interval at once: {every.sum()} of {n}')
    write_log('coverage', lines)
    return


def run_imis(part, tgts, **kwargs):
    res = cb.imis(tgts, **kwargs)
    res.targets = tgts
    res.ncpus = cb.ncpus
    sc.save(fn[part], res)
    tbode4.save_versions(outdir, 'v4')
    return res


def write_log(part, lines):
    text = '\n'.join(lines)
    print(text, flush=True)
    sc.savetext(outdir / f'{part}.log', text)
    return


#%% Summaries

def summarise(part):
    """ Posterior table (weighted quantiles per parameter against the prior; truth checks for re-identification) and the fit to each target """
    res = sc.load(fn[part])
    w = res.weights
    df = outs_df(res)
    last = res.log[-1]
    lines = [f'{part}: {last["n"]:,} model runs, {len(res.comps)} IMIS components, {last["iteration"]} iterations, ESS {last["ess"]:.0f}, '
             f'expected unique {last["expected_unique"]:.0f} of {res.j} (stop at {(1 - np.exp(-1))*res.j:.0f}), {last["time"]/3600:.2f} h on {res.ncpus} workers; '
             f'errors {int(np.nansum(df.error))}; mixture check {res.mixture_check:.1e}']
    rows = []
    for j, k in enumerate(cb.names):
        x = res.theta[:, j]
        lo, hi, _ = cb.priors[k]
        q = weighted_quantile(x, w, [0.025, 0.25, 0.5, 0.75, 0.975])
        row = dict(parameter=k, prior_lo=lo, prior_hi=hi, mean=float(w @ x), q025=q[0], q25=q[1], q50=q[2], q75=q[3], q975=q[4], ci_width_share_of_prior=(q[4] - q[0])/(hi - lo))
        if part == 'reid':
            t = truth[k]
            far = np.abs(x - t) >= abs(q[2] - t) # Draws at least as far from the truth as the posterior median
            row.update(truth=t, truth_in_95=bool(q[0] <= t <= q[4]), truth_position_in_prior=(t - lo)/(hi - lo), max_ll_far=float(np.max(res.ll[far])))
        rows.append(row)
    pars = pd.DataFrame(rows)
    if part == 'reid':
        ll_truth = cb.loglik(cb.simulate(truth_theta()), res.targets)
        pars['max_ll_far_minus_truth'] = pars.max_ll_far - ll_truth
        lines.append(f'Log-likelihood at the truth {ll_truth:.3f}; best point found {np.max(res.ll):.3f}; truth inside its 95% interval for '
                     f'{pars.truth_in_95.sum()} of {len(pars)} parameters')
    pars.to_csv(outdir / f'{part}_posterior.csv', index=False)
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
        q = weighted_quantile(v[ok], w[ok], [0.025, 0.5, 0.975])*(1 if k == 'population_2020' else 1e5)
        rows.append(dict(target=k + ('' if k == 'population_2020' else ' (per 100k)'), post_q025=q[0], post_q50=q[1], post_q975=q[2]))
    fit = pd.DataFrame(rows)
    fit.to_csv(outdir / f'{part}_fit.csv', index=False)
    with pd.option_context('display.width', 250, 'display.max_columns', 30):
        lines += [pars.round(4).to_string(index=False), fit.round(5).to_string(index=False)]
    write_log(f'{part}_summary', lines)
    return pars, fit


def posterior_draws(part, n=200, seed=32):
    """ A weighted resample of parameter sets (with replacement), saved for the agents check, and their trajectories 1950–2040 by calendar year """
    res = sc.load(fn[part])
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(res.weights), size=n, replace=True, p=res.weights)
    draws = pd.DataFrame(res.theta[idx], columns=cb.names)
    draws.insert(0, 'index', idx)
    draws.to_csv(outdir / f'{part}_draws.csv', index=False)
    cals = cb.pool().map(cb.trajectory, list(res.theta[idx]), chunksize=1)
    rows = []
    for i, cal in zip(idx, cals):
        c = cal.loc[1950:2040].reset_index()
        c.insert(0, 'index', int(i))
        rows.append(c)
    tr = pd.concat(rows, ignore_index=True)
    tr.to_csv(outdir / f'{part}_trajectories.csv', index=False)
    return tr


def plot_part(part):
    """ Marginals against the prior, the posterior predictive of each target, and trajectories against the data """
    res = sc.load(fn[part])
    w = res.weights
    fig, axes = plt.subplots(4, 6, figsize=(18, 11))
    for ax, (j, k) in zip(axes.flat, enumerate(cb.names)):
        lo, hi, lg = cb.priors[k]
        bins = np.geomspace(lo, hi, 40) if lg else np.linspace(lo, hi, 40)
        ax.hist(res.theta[:, j], bins=bins, weights=w, density=True, color='#3b6ea8', alpha=0.8)
        if not lg:
            ax.axhline(1/(hi - lo), color='0.5', lw=0.8, ls='--')
        else:
            ax.set_xscale('log')
        if part == 'reid':
            ax.axvline(truth[k], color='#c4620b', lw=1.5)
        ax.set_title(k, fontsize=9)
        ax.set_yticks([])
    for ax in axes.flat[len(cb.names):]:
        ax.axis('off')
    fig.suptitle(f'Exp 32, {part}: posterior marginals (blue; dashed: the uniform prior{"; orange: truth" if part == "reid" else ""})', fontsize=10)
    fig.tight_layout()
    fig.savefig(outdir / f'{part}_marginals.png', dpi=130)
    plt.close(fig)

    df = outs_df(res)
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
    fig.suptitle(f'Exp 32, {part}: posterior predictive of each target (blue) against the target mean and 95% interval (grey)', fontsize=10)
    fig.tight_layout()
    fig.savefig(outdir / f'{part}_fit.png', dpi=130)
    plt.close(fig)

    tr = pd.read_csv(outdir / f'{part}_trajectories.csv')
    panels = sc.objdict({'prev_15+': 'Prevalence 15+ (1 Jan)', 'prevalence': 'Prevalence, all ages (1 Jan)', 'incidence': 'Incidence',
                         'notifications': 'Reported notifications', 'mortality': 'TB mortality', 'alive': 'Population ÷ 1950'})
    fig, axes = plt.subplots(2, 3, figsize=(15, 7.5))
    for ax, (key, title) in zip(axes.flat, panels.items()):
        g = tr.groupby('year')[key]
        lo, md, hi = g.quantile(0.025), g.median(), g.quantile(0.975)
        ax.fill_between(lo.index, lo, hi, color='#3b6ea8', alpha=0.25, lw=0, label='Posterior 95%')
        ax.plot(md.index, md, color='#1f4e8c', lw=1.5, label='Posterior median')
        dk = 'prev15' if key == 'prev_15+' else key
        if dk in data:
            yrs = list(data[dk])
            ax.plot(yrs, [data[dk][y] for y in yrs], 'o', color='k', ms=4, label='Data (targets use ratios)')
        ax.set_title(title + ('' if key == 'alive' else ' per 100,000'), fontsize=10)
        ax.set_xlabel('Year')
        ax.set_ylim(bottom=0)
    axes.flat[0].legend(fontsize=8, frameon=False)
    fig.suptitle(f'Exp 32, {part}: trajectories of 200 posterior draws by calendar year', fontsize=10)
    fig.tight_layout()
    fig.savefig(outdir / f'{part}_trajectories.png', dpi=130)
    plt.close(fig)
    return


def compare():
    """ Seed 1 against seed 2 on parameter marginals and predictions, and both against exp 24's two seeds """
    runs = sc.objdict(real=sc.load(fn.real), real_seed2=sc.load(fn.real_seed2))
    ext = sc.objdict()
    for k, f in [('exp24', exp24 / 'real.obj'), ('exp24_seed2', exp24 / 'real_seed2.obj')]:
        if f.exists():
            ext[k] = sc.load(f)
    allr = sc.mergedicts(runs, ext)
    rows = []
    for j, k in enumerate(cb.names):
        row = dict(parameter=k)
        for name, r in allr.items():
            q = weighted_quantile(r.theta[:, j], r.weights, [0.025, 0.25, 0.5, 0.75, 0.975])
            row.update({f'{name}_q025': q[0], f'{name}_q25': q[1], f'{name}_q50': q[2], f'{name}_q75': q[3], f'{name}_q975': q[4]})
        a, b = 'real', 'real_seed2'
        row['seeds_agree'] = bool(row[f'{b}_q25'] <= row[f'{a}_q50'] <= row[f'{b}_q75'] and row[f'{a}_q25'] <= row[f'{b}_q50'] <= row[f'{a}_q75'])
        rows.append(row)
    pars = pd.DataFrame(rows)
    pars.to_csv(outdir / 'seeds_parameters.csv', index=False)
    rows = []
    for k, t in cb.targets.items():
        row = dict(target=k)
        for name, r in runs.items():
            v = outs_df(r)[k].values
            ok = np.isfinite(v)
            q = weighted_quantile(v[ok], r.weights[ok], [0.025, 0.5, 0.975])
            row.update({f'{name}_q025': q[0], f'{name}_q50': q[1], f'{name}_q975': q[2]})
        rows.append(row)
    preds = pd.DataFrame(rows)
    preds.to_csv(outdir / 'seeds_targets.csv', index=False)
    corr = {}
    for name, r in allr.items():
        jk, ja = cb.names.index('k_latent'), cb.names.index('inf_asy')
        x, y, w = r.theta[:, jk], np.log(r.theta[:, ja]), r.weights
        mx, my = w @ x, w @ y
        corr[name] = float(w @ ((x - mx)*(y - my))/np.sqrt((w @ (x - mx)**2)*(w @ (y - my)**2)))
    lines = ['Seeds compared on parameter marginals (agree: each median inside the other\'s 25–75%):',
             f'{pars.seeds_agree.sum()} of {len(pars)} agree; disagree: {", ".join(pars.parameter[~pars.seeds_agree])}',
             'Weighted correlation of k_latent with log(inf_asy): ' + ', '.join(f'{k} {v:.2f}' for k, v in corr.items())]
    show = ['parameter'] + [f'{n}_{s}' for n in allr for s in ['q025', 'q50', 'q975']] + ['seeds_agree']
    with pd.option_context('display.width', 300, 'display.max_columns', 40):
        lines += [pars[show].round(4).to_string(index=False), preds.round(5).to_string(index=False)]
    write_log('compare', lines)

    fig, axes = plt.subplots(4, 6, figsize=(18, 11))
    for ax, (j, k) in zip(axes.flat, enumerate(cb.names)):
        lo, hi, lg = cb.priors[k]
        bins = np.geomspace(lo, hi, 40) if lg else np.linspace(lo, hi, 40)
        for name, r in allr.items():
            ax.hist(r.theta[:, j], bins=bins, weights=r.weights, density=True, histtype='step', lw=1.3, color=colors[name],
                    ls='-' if name.startswith('real') else '--', label=name)
        if lg:
            ax.set_xscale('log')
        ax.set_title(k, fontsize=9)
        ax.set_yticks([])
    for ax in axes.flat[len(cb.names):]:
        ax.axis('off')
    axes.flat[0].legend(fontsize=7, frameon=False)
    fig.suptitle('Exp 32 posterior marginals, seed 1 (blue) and seed 2 (orange), against exp 24\'s seeds without demography (grey, dashed)', fontsize=10)
    fig.tight_layout()
    fig.savefig(outdir / 'compare_marginals.png', dpi=130)
    plt.close(fig)
    return


def summary():
    for part in ['reid', 'real', 'real_seed2']:
        if fn[part].exists():
            summarise(part)
            if not (outdir / f'{part}_trajectories.csv').exists():
                posterior_draws(part)
            plot_part(part)
    if fn.real.exists() and fn.real_seed2.exists():
        compare()
    return


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'coverage'
    parts = ['coverage', 'reid', 'real', 'real_seed2', 'summary'] if which == 'all' else [which]
    for part in parts:
        T = sc.timer()
        print(f'=== {part} ({cb.ncpus} workers) ===', flush=True)
        if part == 'coverage':
            coverage()
        elif part == 'reid':
            tg, out = synthetic_targets()
            print('Synthetic targets (the truth\'s outputs):', {k: round(float(v), 5) for k, v in out.items()}, flush=True)
            run_imis('reid', tg)
        elif part == 'real':
            run_imis('real', cb.targets)
        elif part == 'real_seed2':
            run_imis('real_seed2', cb.targets, seed=2)
        elif part == 'summary':
            summary()
        print(f'=== {part} done in {T.toc(output=True)/60:.1f} min ===', flush=True)
    if cb._pool is not None: # Closed here, not at interpreter exit, where its teardown raises once modules are gone
        cb._pool.close()
        cb._pool.join()
