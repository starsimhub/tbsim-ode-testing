"""
The v4 track: tbsim v4 (on Starsim v4), where one declared TB model runs as an ODE or as agents.

Not importable alongside tbode (the v0.10 track): each track has its own environment. Results pass between them as saved files.
"""

from .versions import *
from .sims import *
from .latent import *

ncpus = 4 # Cores for sc.parallelize in experiment scripts; runs are on a laptop, so don't raise this
