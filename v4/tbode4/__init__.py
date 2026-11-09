"""
The v4 track: tbsim v4 (on Starsim v4), where one declared TB model runs as an ODE or as agents.

Not importable alongside tbode (the v0.10 track): each track has its own environment. Results pass between them as saved files.
"""

import os
# Starsim's kernels use every core by default; one thread per process is faster here (exp 23). This only takes effect if set before
# starsim (numba) is imported, so scripts set it at their top too; set here, it reaches worker processes started later
os.environ.setdefault('NUMBA_NUM_THREADS', '1')

from .versions import * # noqa: E402
from .sims import *
from .latent import *
from .treat import *
from . import fastode

ncpus = 4 # Cores for sc.parallelize in experiment scripts; runs are on a laptop, so don't raise this
