"""
Experiment 20, v4 environment: calibrate the ODE without treatment to tb_LAI_TPT's six applicable targets, with IMIS. See README.md.
Run from the repo root, in order:

    PYTHONPATH=v4 v4/.venv/bin/python experiments/20_ode_calibration/run.py coverage   # prior predictive: do the prior's outputs bracket the targets?
    PYTHONPATH=v4 v4/.venv/bin/python experiments/20_ode_calibration/run.py reid       # re-identification: IMIS on synthetic targets from a known point
    PYTHONPATH=v4 v4/.venv/bin/python experiments/20_ode_calibration/run.py real       # IMIS on the real targets
    PYTHONPATH=v4 v4/.venv/bin/python experiments/20_ode_calibration/run.py real_seed2 # the same with another seed, to check convergence
    PYTHONPATH=v4 v4/.venv/bin/python experiments/20_ode_calibration/run.py plots      # summaries and figures from the saved runs

Each part saves its raw results (`*.obj`) before summarising.
"""
import sys
import numpy as np
import pandas as pd
import sciris as sc
import matplotlib.pyplot as plt

outdir = sc.thispath(__file__)
sys.path.insert(0, str(outdir))
import calib as cb # noqa: E402
import tbode4 # noqa: E402

fn = sc.objdict(coverage=outdir / 'coverage.obj', reid=outdir / 'reid.obj', real=outdir / 'real.obj', real_seed2=outdir / 'real_seed2.obj')

# The re-identification truth: set B with remote progression at phi = 0.01 (exp 17's model, beta calibrated to 500/100k)
truth = dict(beta=21.14343, trans_asymp=0.9407760328175468, pi=0.2415635040650483, inf_cle=0.5365965386150265, inf_non=0.35069644069791794,
             inf_asy=0.23431100544176175, k_latent=4.631151620621326, phi=0.01, non_rec=0.1706711906227213, non_asy=0.2426759489829796,
             asy_non=1.738594655537684, asy_sym=0.7435451323482238, sym_asy=0.6630323475395543, sym_dead=0.2626559974664169)
labels = dict(prevalence='Prevalence (per 100k)', prop_asymptomatic='Share asymptomatic', mort_to_inc='Mortality ÷ incidence', prev_to_inc='Prevalence ÷ incidence (y)',
              ProgressionY1='ProgressionY1', ProgressionRatio='ProgressionRatio')
scale = dict(prevalence=1e5)


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
    data = sc.objdict(Z=Z, theta=cb.from_unit(Z), outs=df, ll=ll)
    sc.save(fn.coverage, data)
    print(f'{n} prior draws: no endemic equilibrium in {df.prevalence.isna().mean():.1%}')
    rows = []
    for k, t in {**cb.targets, **cb.reported}.items():
        m, lo, hi = target_interval(t)
        v = df[k].dropna().values
        rows.append(dict(target=k, in_likelihood=k in cb.targets, target_mean=m, target_lo=lo, target_hi=hi, prior_p05=np.percentile(v, 5), prior_p50=np.percentile(v, 50),
                         prior_p95=np.percentile(v, 95), share_in_interval=np.mean((v >= lo) & (v <= hi)), share_below=np.mean(v < lo)))
    cov = pd.DataFrame(rows)
    cov.to_csv(outdir / 'coverage.csv', index=False)
    both = np.ones(len(df), bool)
    for k, t in cb.targets.items():
        _, lo, hi = target_interval(t)
        both &= (df[k] >= lo).values & (df[k] <= hi).values
    print(cov.round(4).to_string(index=False))
    print(f'Prior draws inside every target interval at once: {both.sum()} of {n}')
    return data


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
        row = dict(parameter=k, prior_lo=lo, prior_hi=hi, mean=float(w @ x), q025=q[0], q50=q[1], q975=q[2],
                   ci_width_share_of_prior=(q[2] - q[0])/(hi - lo))
        if part == 'reid':
            row.update(truth=truth[k], truth_in_95=bool(q[0] <= truth[k] <= q[2]))
        rows.append(row)
    pars = pd.DataFrame(rows)
    pars.to_csv(outdir / f'{part}_posterior.csv', index=False)
    rows = []
    df = pd.DataFrame([dict(o) for o in res.outs])
    for k, t in {**res.targets, **cb.reported}.items():
        m, lo, hi = target_interval(t)
        v = df[k].values
        ok = np.isfinite(v)
        q = weighted_quantile(v[ok], w[ok], [0.025, 0.5, 0.975])
        rows.append(dict(target=k, in_likelihood=k in res.targets, target_mean=m, target_lo=lo, target_hi=hi, post_q025=q[0], post_q50=q[1], post_q975=q[2]))
    fit = pd.DataFrame(rows)
    fit.to_csv(outdir / f'{part}_fit.csv', index=False)
    last = res.log[-1]
    print(f'{part}: {last["n"]:,} model runs, {len(res.comps)} IMIS components, ESS {last["ess"]:.0f}, expected unique {last["expected_unique"]:.0f} of {res.j}')
    with pd.option_context('display.width', 220):
        print(pars.round(4).to_string(index=False))
        print(fit.round(5).to_string(index=False))
    return pars, fit


def plot(part):
    """ Marginals against the prior, and the posterior predictive of each target against its interval """
    res = sc.load(fn[part])
    w = res.weights
    fig, axes = plt.subplots(3, 5, figsize=(16, 9))
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
    axes.flat[-1].axis('off')
    fig.suptitle(f'Exp 20, {part}: posterior marginals (blue; dashed: prior, approximate for log-uniform{"; orange: truth" if part == "reid" else ""})', fontsize=10)
    fig.tight_layout()
    fig.savefig(outdir / f'{part}_marginals.png', dpi=150)
    df = pd.DataFrame([dict(o) for o in res.outs])
    shown = {**res.targets, **cb.reported}
    fig, axes = plt.subplots(1, len(shown), figsize=(3.2*len(shown), 3))
    for ax, (k, t) in zip(axes, shown.items()):
        m, lo, hi = target_interval(t)
        s = scale.get(k, 1)
        v = df[k].values
        ok = np.isfinite(v)
        q = weighted_quantile(v[ok], w[ok], [0.001, 0.999]) # Range from the posterior (and the target interval), not from every point sampled
        edges = np.linspace(min(q[0], lo), max(q[1], hi), 41)*s
        ax.hist(v[ok]*s, bins=edges, weights=w[ok], density=True, color='#3b6ea8', alpha=0.8)
        ax.set_xlim(edges[0], edges[-1])
        ax.axvspan(lo*s, hi*s, color='0.85', zorder=0)
        ax.axvline(m*s, color='0.2', lw=1)
        ax.set_title(labels[k] + ('' if k in res.targets else ' (not fitted)'), fontsize=9)
        ax.set_yticks([])
    fig.suptitle(f'Exp 20, {part}: posterior predictive of each target (blue) against the target mean and 95% interval (grey)', fontsize=10)
    fig.tight_layout()
    fig.savefig(outdir / f'{part}_fit.png', dpi=150)
    print(f'Saved {part}_marginals.png, {part}_fit.png')
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
        run_imis('real', cb.targets, max_iter=300)
        summarise('real')
        plot('real')
    elif which == 'real_seed2': # The same calibration with another seed: are the posteriors the same?
        run_imis('real_seed2', cb.targets, max_iter=300, seed=2)
        summarise('real_seed2')
    elif which == 'plots':
        for part in ['reid', 'real']:
            if fn[part].exists():
                summarise(part)
                plot(part)
