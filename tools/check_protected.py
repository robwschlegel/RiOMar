#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Refuse changes that delete code or outputs the manuscript may depend on.

    python tools/check_protected.py                    # uncommitted changes vs HEAD
    python tools/check_protected.py --range A..B       # a commit range
    python tools/check_protected.py --hook             # Claude Code PreToolUse hook

Blocks, unless the item is listed in metadata/approved_deletions.txt:
  1. deleting (or untracking) any tracked file under code/, func/, tools/,
     metadata/, figures/ or animations/, or any file named in
     metadata/paragraph_source_registry.csv / figure_table_registry.csv
     (whatever that row's status says);
  2. removing an R (`name <- function`) or Python (`def name`) function
     definition.
Moves and renames pass: a deleted file whose exact content still exists at
another path, or a removed function matched by an added one with the same
argument list, counts as moved, not deleted.

Why: manuscript/ is gitignored, so an agent cannot see what the text still
cites; "looks unused" is not evidence of "is unused" (see CLAUDE.md, "Rules
for automated changes"). Deletions need Robert's per-item approval, recorded
as a line in metadata/approved_deletions.txt.
"""

import argparse
import csv
import json
import os
import re
import subprocess
import sys

proj_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GUARDED_DIRS = ('code/', 'func/', 'tools/', 'metadata/', 'figures/', 'animations/')
APPROVED_PATH = os.path.join(proj_dir, 'metadata', 'approved_deletions.txt')
# group 1 = function name, group 2 = argument list (used to recognise renames)
R_DEF = re.compile(r'^\s*([A-Za-z_.][A-Za-z0-9_.]*)\s*(?:<-|=)\s*function\s*(\(.*)$')
PY_DEF = re.compile(r'^\s*def\s+([A-Za-z_][A-Za-z0-9_]*)\s*(\(.*)$')


def git(*args):
    return subprocess.run(['git', *args], cwd=proj_dir, capture_output=True, text=True, check=True).stdout


def approved_items():
    if not os.path.exists(APPROVED_PATH):
        return set()
    with open(APPROVED_PATH) as f:
        return {line.split('#')[0].strip() for line in f if line.split('#')[0].strip()}


def registry_paths():
    """Every repo path either registry names (source files, check paths, scripts)."""
    paths = set()
    para = os.path.join(proj_dir, 'metadata', 'paragraph_source_registry.csv')
    if os.path.exists(para):
        for row in csv.DictReader(open(para, newline='')):
            for entry in (row.get('source_files') or '').split(';'):
                entry = entry.strip()
                if entry.startswith(('EXTERNAL-CITATION:', 'UNVERIFIED:')):
                    continue
                paths.add(entry.split(':', 1)[1] if entry.startswith('PLOT-ONLY:') else entry)
    fig = os.path.join(proj_dir, 'metadata', 'figure_table_registry.csv')
    if os.path.exists(fig):
        for row in csv.DictReader(open(fig, newline='')):
            paths.add((row.get('check_path') or '').strip())
            for col in ('r_function', 'python_entry_point'):
                for name in re.findall(r'[\w./-]+\.(?:R|py)\b', row.get(col) or ''):
                    paths.add(name if '/' in name else f'func/{name}')
    return {p for p in paths if p}


def new_blobs(base, target):
    """Blob ids of every file in the post-change state (a commit, or the working tree)."""
    if target:
        return {line.split()[2] for line in git('ls-tree', '-r', target).splitlines()}
    blobs = {line.split()[1] for line in git('ls-files', '-s').splitlines()}
    changed = [p for p in git('diff', base, '--name-only', '--diff-filter=AMR').splitlines()
               if os.path.isfile(os.path.join(proj_dir, p))]
    if changed:
        blobs |= set(git('hash-object', '--', *changed).split())
    return blobs


def deleted_files(base, target):
    """Deleted paths whose content does not survive anywhere else (i.e. real deletions, not moves)."""
    diff_args = [base, target] if target else [base]
    deleted = [line.split('\t')[1] for line in git('diff', '--name-status', '-M', *diff_args).splitlines()
               if line.startswith('D\t')]
    if not deleted:
        return []
    surviving = new_blobs(base, target)
    return [p for p in deleted if git('rev-parse', f'{base}:{p}').strip() not in surviving]


def removed_functions(base, target):
    """Function names whose definition is removed and not matched by an added definition
    with the same argument list (a rename)."""
    diff_args = [base, target] if target else [base]
    removed, added = [], []
    for line in git('diff', '-M', '--unified=0', *diff_args, '--', '*.R', '*.py').splitlines():
        if line[:1] not in '+-' or line.startswith(('---', '+++')):
            continue
        for pattern in (R_DEF, PY_DEF):
            m = pattern.match(line[1:])
            if m:
                (added if line[0] == '+' else removed).append((m.group(1), m.group(2).strip()))
    added_names = {name for name, _ in added}
    spare_signatures = [sig for name, sig in added]
    result = []
    for name, sig in removed:
        if name in added_names:
            continue  # redefined (edited in place, or moved to another file)
        if sig in spare_signatures:
            spare_signatures.remove(sig)  # renamed: same arguments, new name
            continue
        result.append(name)
    return sorted(set(result))


def check(base, target=None):
    approved = approved_items()
    protected = registry_paths()
    problems = []
    for path in deleted_files(base, target):
        if path in approved:
            continue
        if path in protected:
            problems.append(f'deletes {path} -- named in a manuscript registry')
        elif path.startswith(GUARDED_DIRS):
            problems.append(f'deletes {path}')
    for name in removed_functions(base, target):
        if name not in approved:
            problems.append(f'removes function {name}()')
    return problems


def report(problems):
    print('BLOCKED by tools/check_protected.py -- this change deletes things that may be in use:', file=sys.stderr)
    for p in problems:
        print(f'  - {p}', file=sys.stderr)
    print('Moves/renames are fine. To delete, get Robert\'s explicit per-item approval and add each path or '
          'function name to metadata/approved_deletions.txt (see CLAUDE.md, "Rules for automated changes").',
          file=sys.stderr)


def hook_range(command):
    """(base, target) a git commit/push command would publish; None if it is neither or has nothing to push."""
    if re.search(r'\bgit\b[^;&|]*\bcommit\b', command):
        return ('HEAD', None)  # staged + unstaged vs HEAD, so `git add -A && git commit` is covered too
    if re.search(r'\bgit\b[^;&|]*\bpush\b', command):
        # Only the commits this push would publish: those on no remote branch yet
        unpushed = git('rev-list', 'HEAD', '--not', '--remotes').split()
        return (f'{unpushed[-1]}~1', 'HEAD') if unpushed else None
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--range', help='commit range A..B to check (default: uncommitted changes vs HEAD)')
    parser.add_argument('--hook', action='store_true', help='read a Claude Code PreToolUse payload on stdin')
    args = parser.parse_args()

    if args.hook:
        payload = json.load(sys.stdin)
        rng = hook_range(payload.get('tool_input', {}).get('command', ''))
        if rng is None:  # not a commit/push, or nothing unpushed
            return
        problems = check(*rng)
        if problems:
            report(problems)
            sys.exit(2)  # exit code 2 = block the tool call and show stderr to Claude
        return

    problems = check(*args.range.split('..', 1)) if args.range else check('HEAD')
    if problems:
        report(problems)
        sys.exit(1)
    print('OK: no unapproved deletions.')


if __name__ == '__main__':
    main()
