""" Figure for treatment_parameters: the suggested shapes of partial benefit after loss to follow-up, q(f) = f^a """
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

out = sys.argv[1] if len(sys.argv) > 1 else 'qf_shapes.png'
D = 6 # Course length, months
m = np.linspace(0, D, 301)
fig, ax = plt.subplots(figsize=(6.5, 3.6))
for a, style, label in [(0.3, '-', 'a = 0.3 (prior lower bound)'), (0.5, '-', 'a = 0.5'), (1.0, '--', 'a = 1 (exp 35 placeholder, q = f)')]:
    ax.plot(m, (m/D)**a, style, lw=2, label=label)

# What the model uses with 2-month phases: q at each phase's midpoint (1, 3, 5 months), for a = 0.5
edges, mids = np.array([0, 2, 4, 6]), np.array([1, 3, 5])
ax.stairs((mids/D)**0.5, edges, baseline=None, color='C1', alpha=0.35, lw=6, label='a = 0.5 as run (2-month phases)')
ax.set_xlabel('Months of treatment completed when lost to follow-up')
ax.set_ylabel('Chance of cure, relative to\ncompleting the course (q)')
ax.set_xlim(0, D)
ax.set_ylim(0, 1.02)
ax.set_xticks(range(D + 1))
ax.grid(alpha=0.3)
ax.legend(frameon=False, fontsize=8, loc='lower right')
fig.tight_layout()
fig.savefig(out, dpi=200)
print('saved', out)
