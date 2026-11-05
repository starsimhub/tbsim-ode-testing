"""
Latent-progression structures on one person's continuous-time Markov chain, solved exactly (expm): for cohort comparisons of
progression by time since infection (exps 14 on). numpy/scipy only.

Structures:
- 'exact':     tb_LAI_TPT / v0.10's hazard r·exp(−k·τ) (k_asy = k_non = k), as Poisson-mixture classes (exact). Reinfection of a
               latent person reassigns them across the classes with the Poisson weights, which is an exact clock reset.
- 'two-state': v4's recent and remote latent states (`k_latent`). Recent infections progress at inf_non / inf_asy and become remote at
               k_latent; remote ones progress at inf_non_remote / inf_asy_remote and clear at inf_cle, like recent ones. Reinfection
               moves remote → recent and leaves recent unchanged.

States: the structure's latent states, then CLE (cleared), NON (non-infectious TB), REC (recovered), ASY (asymptomatic TB), SYM
(symptomatic TB), DEAD. By default (cascade=False) ASY is absorbing, so its occupancy is 'ever reached asymptomatic TB' (tb_LAI_TPT's
curve B) and SYM is never reached. With cascade=True, active TB moves on (ASY -> NON, ASY -> SYM, SYM -> ASY, SYM -> DEAD), and
`absorbing` makes chosen states absorbing, for first-passage curves (see first_passage() and latent_exits()).
"""
import numpy as np
import sciris as sc
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply
from scipy.stats import poisson

__all__ = ['n_classes', 'latent_chain', 'cohort_curve', 'ever_asy', 'latent_hazard', 'progression_targets', 'first_passage', 'latent_exits',
           'target_means', 'target_sigmas']

# tb_LAI_TPT's progression targets (cohort_progression.py): ProgressionY1 = B(1), ProgressionRatio = B(1)/B(10); sigmas are log-space
target_means = dict(ProgressionY1=0.0523, ProgressionRatio=0.663)
target_sigmas = dict(ProgressionY1=0.094, ProgressionRatio=0.071)


def n_classes(r0, k, tol=1e-12):
    """ The number of Poisson-mixture classes (0..M-1) whose truncated weight is below tol """
    return max(2, int(poisson.isf(tol, r0/k)) + 2)


def latent_chain(p, structure='exact', lam=0.0, latent_reinf=True, mu=0.0, tol=1e-12, cascade=False, absorbing=None):
    """
    One person's generator Q, initial vector x0 (everyone just infected), and the state indices

    Args:
        p (dict): rates per year: inf_cle, inf_non, inf_asy, non_rec, non_asy, pi (reinfection risk ratio for everyone); and k (exact) or
            k_latent, inf_non_remote, inf_asy_remote (two-state; remote rates default 0); with cascade, also asy_non, asy_sym, sym_asy
            and sym_dead (default 0)
        structure (str): 'exact' or 'two-state'
        lam (float): force of infection, per year: people no longer infected are reinfected at pi·lam, and latent people at pi·lam by
            the structure's own rule (if latent_reinf)
        mu (float): background deaths, per year
        tol (float): truncated Poisson weight in the exact classes
        cascade (bool): whether asymptomatic TB moves on (False: ASY is absorbing, as for tb_LAI_TPT's curve B)
        absorbing (list): names of states to make absorbing ('LATENT' for all latent states, 'CLE', 'NON', 'REC', 'ASY', 'SYM'), e.g.
            ['NON', 'ASY'] for first progression

    Returns:
        Q (array), x0 (array), idx (objdict: latent (list of indices), CLE, NON, REC, ASY, SYM, DEAD, prog (progression rate of each
        latent state), and for the two-state structure recent and remote)
    """
    r0 = p['inf_non'] + p['inf_asy']
    if structure == 'exact':
        nL = n_classes(r0, p['k'], tol)
        entry = poisson.pmf(np.arange(nL), r0/p['k'])
        entry[-1] += 1 - entry.sum()
        non = np.arange(nL)*p['k']*p['inf_non']/r0 # Class m progresses at m·k, split between NON and ASY as inf_non : inf_asy
        asy = np.arange(nL)*p['k']*p['inf_asy']/r0
    elif structure == 'two-state':
        nL = 2
        entry = np.array([1.0, 0.0])
        non = np.array([p['inf_non'], p.get('inf_non_remote', 0.0)])
        asy = np.array([p['inf_asy'], p.get('inf_asy_remote', 0.0)])
    else:
        errormsg = f'structure is "exact" or "two-state", not {structure!r}'
        raise ValueError(errormsg)
    CLE, NON, REC, ASY, SYM, DEAD = range(nL, nL + 6)
    n = nL + 6
    Q = np.zeros((n, n))

    def add(i, j, rate):
        Q[i, j] += rate
        Q[i, i] -= rate

    for i in range(nL):
        add(i, CLE, p['inf_cle'])
        add(i, NON, non[i])
        add(i, ASY, asy[i])
    if structure == 'two-state':
        add(0, 1, p['k_latent']) # Recent -> remote
    add(NON, ASY, p['non_asy'])
    add(NON, REC, p['non_rec'])
    if cascade:
        add(ASY, NON, p['asy_non'])
        add(ASY, SYM, p['asy_sym'])
        add(SYM, ASY, p['sym_asy'])
        add(SYM, DEAD, p.get('sym_dead', 0.0))
    if lam > 0:
        for src in [CLE, REC]: # New infections of people no longer infected
            for j in range(nL):
                add(src, j, p['pi']*lam*entry[j])
        if latent_reinf:
            if structure == 'exact': # Clock reset: reassigned across the classes with the Poisson weights
                for i in range(nL):
                    for j in range(nL):
                        if i != j:
                            add(i, j, p['pi']*lam*entry[j])
            else: # Remote -> recent; reinfecting a recent infection does nothing
                add(1, 0, p['pi']*lam)
    if mu > 0:
        for i in range(n):
            if i != DEAD and (cascade or i != ASY):
                add(i, DEAD, mu)
    names = dict(CLE=[CLE], NON=[NON], REC=[REC], ASY=[ASY], SYM=[SYM], LATENT=list(range(nL)))
    for name in absorbing or []:
        for i in names[name]:
            Q[i, :] = 0.0 # No exits: the state is absorbing
    x0 = np.zeros(n)
    x0[:nL] = entry
    idx = sc.objdict(latent=list(range(nL)), CLE=CLE, NON=NON, REC=REC, ASY=ASY, SYM=SYM, DEAD=DEAD, prog=non + asy)
    if structure == 'two-state':
        idx.recent, idx.remote = 0, 1
    return Q, x0, idx


def cohort_curve(p, structure='exact', tmax=50, per_year=48, **kwargs):
    """ The state occupancies of a cohort on the grid 0, 1/per_year, ..., tmax years: (t, X (times × states), idx) """
    Q, x0, idx = latent_chain(p, structure, **kwargs)
    t = np.linspace(0, tmax, tmax*per_year + 1)
    X = expm_multiply(Q.T, x0, start=0, stop=tmax, num=len(t), endpoint=True)
    return t, X, idx


def ever_asy(p, structure='exact', years=(1, 10), **kwargs):
    """ The share ever reaching asymptomatic TB at the given years """
    Q, x0, idx = latent_chain(p, structure, **kwargs)
    return np.array([(x0 @ expm(Q*t))[idx.ASY] for t in years])


def latent_hazard(X, idx):
    """ The progression hazard among people still latent, at each time of a cohort_curve(), per year """
    L = X[:, idx.latent]
    tot = L.sum(axis=1)
    return np.divide(L @ idx.prog, tot, out=np.full(len(tot), np.nan), where=tot > 0)


def first_passage(p, structure='exact', curve='B', tmax=50, per_year=48):
    """
    tb_LAI_TPT's progression curves for a no-reinfection cohort (its ProgressionByTimeTracker), on the grid 0, 1/per_year, ..., tmax:
    'A' first progression to active TB (NON or ASY), 'B' ever reaching ASY, 'C' ever reaching SYM. Returns (t, curve)
    """
    targets = dict(A=['NON', 'ASY'], B=['ASY'], C=['SYM'])[curve]
    t, X, idx = cohort_curve(p, structure, tmax=tmax, per_year=per_year, cascade=True, absorbing=targets)
    return t, X[:, [idx[s] for s in targets]].sum(axis=1)


def latent_exits(p, structure='exact', tmax=50, per_year=48):
    """
    How a no-reinfection cohort first leaves latent infection (tb_LAI_TPT's InfectionOutcomeTracker): the cumulative shares that have
    left to CLE, NON and ASY, and the share still latent, on the grid 0, 1/per_year, ..., tmax. Returns (t, objdict of arrays)
    """
    t, X, idx = cohort_curve(p, structure, tmax=tmax, per_year=per_year, cascade=True, absorbing=['CLE', 'NON', 'ASY'])
    out = sc.objdict(cleared=X[:, idx.CLE], non=X[:, idx.NON], asy=X[:, idx.ASY], still_latent=X[:, idx.latent].sum(axis=1))
    return t, out


def progression_targets(p, structure='exact'):
    """ tb_LAI_TPT's progression observables: ProgressionY1 = B(1), ProgressionRatio = B(1)/B(10) (no reinfection, no deaths) """
    b1, b10 = ever_asy(p, structure, years=(1, 10))
    return dict(ProgressionY1=float(b1), ProgressionRatio=float(b1/b10))
