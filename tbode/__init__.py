"""
ODE versions of TBsim, plus helpers for building matched agent-based runs and comparing the two.
"""

from .ode import *
from .abm import *
from .exact import *
from .compare import *
from .flux import *

ncpus = 4 # Cores for sc.parallelize in experiment scripts; runs are on a laptop, so don't raise this
