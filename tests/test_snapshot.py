# Author: Michel Sanner, PhD. Copyright (c) 2026.
"""A snapshot of ADapt made with `git archive` (how an installer ships it,
without .git) still records which commit it is: adapt/_build.py is stamped by
export-subst, and env.provenance falls back to it. Before, adapt_commit was
None outside a git checkout."""
import os
import subprocess
import sys

from common import check, outdir

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)


def test_archived_snapshot_records_its_commit():
    head = subprocess.check_output(['git', '-C', REPO, 'rev-parse', 'HEAD'],
                                   universal_newlines=True).strip()
    dest = os.path.join(outdir('snapshot'), 'adapt')
    os.makedirs(dest)
    archive = subprocess.Popen(['git', '-C', REPO, 'archive', 'HEAD'], stdout=subprocess.PIPE)
    subprocess.check_call(['tar', '-x', '-C', dest], stdin=archive.stdout)
    check(archive.wait() == 0, 'git archive failed')
    check(not os.path.exists(os.path.join(dest, '.git')), 'the snapshot has a .git')
    code = ('import json, sys; from adapt import env; p = env.provenance(); '
            'print(json.dumps([p["adapt_commit"], p["adapt_source"]]))')
    out = subprocess.check_output([sys.executable, '-c', code], cwd=dest,
                                  env=dict(os.environ, PYTHONPATH=dest),
                                  universal_newlines=True, stderr=subprocess.DEVNULL)
    commit, source = __import__('json').loads(out.strip().splitlines()[-1])
    check(commit == head and source == 'snapshot', (commit, source, head))
