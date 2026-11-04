"""
The versions a run used, saved beside its outputs (as versions.json), so that each result can be traced to its code.

The same file is in both tracks (tbode/versions.py for v0.10, tbode4/versions.py for v4), since their environments can't import each other.
"""
import sys
import json
import datetime
import subprocess
import importlib.metadata as im
from pathlib import Path

__all__ = ['versions', 'save_versions']


def package_info(name):
    """ The version and, for a package installed from git (as both tracks pin tbsim and starsim), its commit """
    try:
        dist = im.distribution(name)
    except im.PackageNotFoundError:
        return None
    out = dict(version=dist.version)
    direct = dist.read_text('direct_url.json')
    if direct:
        d = json.loads(direct)
        out['source'] = d.get('url')
        vcs = d.get('vcs_info')
        if vcs:
            out['commit'] = vcs.get('commit_id')
        elif d.get('dir_info', {}).get('editable'):
            out['editable'] = True # No commit to record: a working copy
    return out


def repo_info():
    """ This repo's commit, and whether the working tree had changes (which the commit doesn't capture) """
    root = Path(__file__).resolve().parent
    git = lambda *args: subprocess.run(['git', *args], cwd=root, capture_output=True, text=True).stdout.strip()
    return dict(commit=git('rev-parse', 'HEAD'), dirty=bool(git('status', '--porcelain')))


def versions(track):
    """
    The versions in use: Python, tbsim and starsim (with their commits), and this repo

    Args:
        track (str): 'v0.10' or 'v4'
    """
    return dict(
        track   = track,
        date    = datetime.datetime.now().isoformat(timespec='seconds'),
        python  = sys.version.split()[0],
        tbsim   = package_info('tbsim'),
        starsim = package_info('starsim'),
        repo    = repo_info(),
    )


def save_versions(folder, track):
    """ Write versions.json into an experiment's folder, and return what was written """
    out = versions(track)
    path = Path(folder) / 'versions.json'
    path.write_text(json.dumps(out, indent=2) + '\n')
    return out
