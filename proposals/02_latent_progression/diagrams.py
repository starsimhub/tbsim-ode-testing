"""
Diagrams of the latent-infection structures compared in options.md (A-D).

Each shows only transitions within, and exits from, latent infection: progression (up, to NON or ASY),
clearance (down), moves between latent states (across) and reinfection (dashed). All rates are per year.
Background deaths are not shown.

Run from the repo root: v4/.venv/bin/python proposals/02_latent_progression/diagrams.py
"""
import os
import textwrap
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

folder = os.path.dirname(os.path.abspath(__file__))

# Colours: fills carry the node type; all text stays in ink
ink = '#1f2328'
muted = '#656d76'
latent_fill = '#dbe9f6'
latent_edge = '#4a7fb5'
out_fill = '#eef0f2'
out_edge = '#8c959f'
reinf = '#c4620b'

w, h = 1.15, 0.8    # latent node size
y_lat = 3.0         # latent row
y_top = (4.7, 5.5)  # progression box
y_bot = (0.5, 1.3)  # cleared box
x_span = (1.9, 9.9)


def box(ax, x, y, text, fill=latent_fill, edge=latent_edge, width=w, height=h, fontsize=10, dashed=False):
    patch = FancyBboxPatch((x - width/2, y - height/2), width, height, boxstyle='round,pad=0.02,rounding_size=0.12',
                           fc=fill, ec=edge, lw=1.2, ls='--' if dashed else '-')
    ax.add_patch(patch)
    ax.text(x, y, text, ha='center', va='center', color=ink, fontsize=fontsize)


def arrow(ax, start, end, color=ink, dashed=False, rad=0.0, lw=1.2):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle='-|>', mutation_scale=12, color=color, lw=lw,
                                 ls='--' if dashed else '-', connectionstyle=f'arc3,rad={rad}', shrinkA=0, shrinkB=0))


def label(ax, x, y, text, ha='left', color=ink, fontsize=9):
    ax.text(x, y, text, ha=ha, va='center', color=color, fontsize=fontsize,
            bbox=dict(fc='white', ec='none', pad=1.0, alpha=0.9))


def frame(title, subtitle, note):
    """ Figure with the shared outcome boxes and the entry node """
    fig, ax = plt.subplots(figsize=(8.4, 4.6))
    ax.set_xlim(0, 10.2)
    ax.set_ylim(-0.4, 6.6)
    ax.axis('off')
    fig.text(0.02, 0.96, title, ha='left', va='top', fontsize=12.5, weight='bold', color=ink)
    fig.text(0.02, 0.905, textwrap.fill(subtitle, 120), ha='left', va='top', fontsize=9.5, color=muted)
    xm, wid = sum(x_span)/2, x_span[1] - x_span[0]
    box(ax, xm, sum(y_top)/2, 'Progression: to NON or ASY, split as inf_non : inf_asy', fill=out_fill, edge=out_edge, width=wid, height=y_top[1]-y_top[0])
    box(ax, xm, sum(y_bot)/2, 'Cleared', fill=out_fill, edge=out_edge, width=wid, height=y_bot[1]-y_bot[0])
    ax.text(0.75, y_lat, 'New\ninfection', ha='center', va='center', color=ink, fontsize=10)
    fig.text(0.02, 0.03, textwrap.fill(note, 135), ha='left', va='bottom', fontsize=8.5, color=muted)
    return fig, ax


def up(ax, x, text, left=False):
    arrow(ax, (x, y_lat + h/2), (x, y_top[0]))
    label(ax, x + (-0.08 if left else 0.08), (y_lat + h/2 + y_top[0])/2, text, ha='right' if left else 'left')


def down(ax, x, text='inf_cle'):
    arrow(ax, (x, y_lat - h/2), (x, y_bot[1]))
    label(ax, x + 0.08, (y_lat - h/2 + y_bot[1])/2, text)


def across(ax, x0, x1, text=None, below=False):
    arrow(ax, (x0 + w/2, y_lat), (x1 - w/2, y_lat))
    if text:
        label(ax, (x0 + x1)/2, y_lat + (-0.28 if below else 0.28), text, ha='center')


def save(fig, name):
    for ext in ['png', 'svg']:
        fig.savefig(os.path.join(folder, f'{name}.{ext}'), dpi=200)
    plt.close(fig)
    print(f'Saved {name}.png/.svg')


def option_a():
    fig, ax = frame('A. Poisson-mixture latent classes',
                    'Class fixed at infection; class m progresses at m·k. Across classes, the hazard is exactly r₀·e^(−kτ).',
                    'r₀ = inf_non + inf_asy; p_m = Poisson(m; r₀/k). Class 0 never progresses. M ≈ 8 for set B (truncation 1e-10; 4–5 classes give 1e-6). Background deaths not shown.')
    xs = [2.9, 4.5, 6.1, 7.4, 8.9]
    names = ['Class 0', 'Class 1', 'Class 2', None, 'Class M']
    ups = [None, 'k', '2k', None, 'M·k']
    gx0, gx1 = xs[0] - w/2 - 0.25, xs[-1] + w/2 + 0.25
    box(ax, (gx0 + gx1)/2, y_lat, '', fill='none', edge=muted, width=gx1 - gx0, height=h + 0.5, dashed=True)
    for x, name, u in zip(xs, names, ups):
        if name is None:
            ax.text(x, y_lat, '…', ha='center', va='center', color=ink, fontsize=14)
            continue
        box(ax, x, y_lat, name)
        down(ax, x)
        if u:
            up(ax, x, u)
        else:
            label(ax, x, y_lat + h/2 + 0.6, 'no progression', ha='center', color=muted)
    arrow(ax, (1.2, y_lat), (gx0, y_lat))
    label(ax, (1.2 + gx0)/2, y_lat + 0.3, 'p_m', ha='center')
    arrow(ax, (gx0, y_lat - 0.35), (0.75, y_lat - 0.45), color=reinf, dashed=True, rad=-0.6)
    label(ax, 0.0, y_lat - 1.8, 'Reinfection π·λ,\nfrom any class:\nreassigned with p_m\n(= clock reset)', color=reinf, fontsize=8)
    save(fig, 'diagram_A_classes')


def option_b():
    fig, ax = frame('B. Two latent states, with remote progression (v4\'s k_latent)',
                    'Recent infections become remote at k_latent; remote infections progress at a fraction φ of the recent rates.',
                    'φ = inf_non_remote / inf_non = inf_asy_remote / inf_asy. φ = 0 is the set-aside single-state version. Background deaths not shown.')
    xr, xm = 3.7, 7.6
    box(ax, xr, y_lat, 'Recent')
    box(ax, xm, y_lat, 'Remote')
    arrow(ax, (1.2, y_lat), (xr - w/2, y_lat))
    up(ax, xr, 'r₀ = inf_non + inf_asy', left=True)
    up(ax, xm, 'φ·r₀')
    down(ax, xr)
    down(ax, xm)
    across(ax, xr, xm, 'k_latent', below=True)
    arrow(ax, (xm - w/2 + 0.15, y_lat + h/2), (xr + w/2 - 0.15, y_lat + h/2), color=reinf, dashed=True, rad=0.45)
    label(ax, (xr + xm)/2, y_lat + 1.0, 'Reinfection π·λ', ha='center', color=reinf, fontsize=8.5)
    label(ax, xr - w/2 - 0.05, y_lat - h/2 - 0.3, 'Reinfection of recent: no change', ha='right', color=reinf, fontsize=8.5)
    save(fig, 'diagram_B_two_states')


def option_c():
    fig, ax = frame('C. A chain of recent stages',
                    'Recent stages with declining progression rates r₁ > r₂ > … > rₙ, then remote. Generalises B; approaches any shape with enough stages.',
                    'Stage rates fitted, or tied to k by a formula. φ·r₀ may be 0. Reinfection: a choice of which stage to return to (not drawn). Background deaths not shown.')
    xs = [2.75, 4.45, 5.95, 7.45, 9.15]
    names = ['Recent 1', 'Recent 2', None, 'Recent n', 'Remote']
    ups = ['r₁', 'r₂', None, 'rₙ', 'φ·r₀']
    arrow(ax, (1.2, y_lat), (xs[0] - w/2, y_lat))
    for i, (x, name, u) in enumerate(zip(xs, names, ups)):
        if name is None:
            ax.text(x, y_lat, '…', ha='center', va='center', color=ink, fontsize=14)
            continue
        box(ax, x, y_lat, name)
        up(ax, x, u)
        down(ax, x)
    across(ax, xs[0], xs[1], 'γ₁')
    arrow(ax, (xs[1] + w/2, y_lat), (xs[2] - 0.25, y_lat))
    arrow(ax, (xs[2] + 0.25, y_lat), (xs[3] - w/2, y_lat))
    across(ax, xs[3], xs[4], 'γₙ')
    save(fig, 'diagram_C_stage_chain')


def option_d():
    fig, ax = frame('D. A non-exponential dwell time for progression',
                    'One latent state. Progression happens after a dwell time drawn from a distribution F, competing with clearance.',
                    'Agents sample F exactly; the ODE approximates it with Erlang stages (or solver=\'ide\'), so the two solvers differ. Background deaths not shown.')
    x = 5.9
    box(ax, x, y_lat, 'Latent')
    arrow(ax, (1.2, y_lat), (x - w/2, y_lat))
    up(ax, x, 'after T ~ F (not exponential)')
    down(ax, x, 'inf_cle (whichever comes first)')
    save(fig, 'diagram_D_dwell_time')


if __name__ == '__main__':
    option_a()
    option_b()
    option_c()
    option_d()
