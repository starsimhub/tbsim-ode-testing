"""
pytest setup for v4/tests: one numba thread per process, set before any test module imports starsim (as tbode4 does on import). Otherwise a test
module that imports starsim first starts numba with every core, and agent runs in later tests fail when tbode4 asks for one thread.
"""
import os
os.environ.setdefault('NUMBA_NUM_THREADS', '1')
