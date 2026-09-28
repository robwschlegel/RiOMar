#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check that every reference to an R/Python file in the repo points at a file
that exists -- a static substitute for re-running the pipeline after files
are moved.

    python tools/check_paths.py                  # list unresolved references
    python tools/check_paths.py --save base.txt  # record them (before a move)
    python tools/check_paths.py --against base.txt   # fail on any new ones (after)

Scans tracked code (code/, func/, tools/), metadata/*.R and both metadata
registries for:
  - path references containing '/' and ending in .R/.py (source("func/x.R"),
    'func/tools/run_x.R', registry source_files, comments too) -- resolved
    from the repo root, or from func/ (code joins 'analysis/x.R' onto func_dir);
  - bare 'name.R' / 'name.py' string literals in code lines (not comments)
    -- resolved against func/, since the pipeline joins them onto func_dir;
  - Python `import name` / `from name import` of a module that exists
    somewhere under func/ or tools/ -- the importing file must put that
    module's folder on sys.path.
Historical mentions of long-gone files (e.g. in comments) are expected, so
the useful signal is the difference before vs. after a change (--against).
"""

import argparse
import os
import re
import subprocess
import sys

proj_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCAN = ('code/', 'func/', 'tools/', 'metadata/')
SCAN_EXT = ('.R', '.py', '.csv')
PATH_REF = re.compile(r'(?<![\w./-])((?:[\w.-]+/)+[\w.-]+\.(?:R|py))\b')
BARE_LITERAL = re.compile(r'''["']([\w.-]+\.(?:R|py))["']''')
PY_IMPORT = re.compile(r'^\s*(?:from\s+([\w.]+)\s+import|import\s+([\w., ]+))')


def tracked_files():
    out = subprocess.run(['git', 'ls-files'], cwd=proj_dir, capture_output=True, text=True, check=True).stdout
    return out.splitlines()


def code_part(line, ext):
    """The line with any trailing comment removed (good enough: '#' inside strings is rare here)."""
    return line.split('#', 1)[0] if ext in ('.R', '.py') else line


def unresolved_references():
    files = tracked_files()
    existing = set(files)
    modules = {os.path.splitext(f)[0].split('/')[-1]: os.path.dirname(f)
               for f in files if f.endswith('.py') and f.startswith(('func/', 'tools/'))}
    problems = set()
    for rel in files:
        if not (rel.startswith(SCAN) and rel.endswith(SCAN_EXT)) or rel == 'tools/check_paths.py':
            continue
        ext = os.path.splitext(rel)[1]
        text = open(os.path.join(proj_dir, rel), encoding='utf-8', errors='replace').read()
        for n, line in enumerate(text.splitlines(), 1):
            for ref in PATH_REF.findall(line):
                if ref.startswith(('output/', 'figures/', 'data/', 'manuscript/', 'http', '~')):
                    continue
                if ref not in existing and f'func/{ref}' not in existing:  # code often joins onto func_dir
                    problems.add(f'{rel}: missing file {ref}')
            if ext in ('.R', '.py'):
                code = code_part(line, ext)
                for name in BARE_LITERAL.findall(code):
                    if re.search(r'''["'](analysis|tools)["']''', code):
                        continue  # os.path.join(func_dir, 'analysis', 'x.R') style: subfolder given
                    if f'func/{name}' not in existing and any(f.endswith('/' + name) for f in files):
                        problems.add(f'{rel}: bare "{name}" no longer in func/ (now in a subfolder)')
            if ext == '.py':
                m = PY_IMPORT.match(code_part(line, ext))
                names = [m.group(1)] if m and m.group(1) else (m.group(2).split(',') if m else [])
                for name in (x.strip().split(' ')[0].split('.')[0] for x in names):
                    folder = modules.get(name)
                    if folder and folder not in ('func', os.path.dirname(rel)) and \
                            os.path.basename(folder) not in text:
                        problems.add(f'{rel}: imports {name} from {folder}/ but never adds it to sys.path')
    return problems


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--save', help='write the unresolved references to this file')
    parser.add_argument('--against', help='fail if any reference is unresolved that was not in this saved file')
    args = parser.parse_args()

    problems = unresolved_references()
    if args.save:
        with open(args.save, 'w') as f:
            f.write('\n'.join(sorted(problems)) + '\n')
        print(f'Saved {len(problems)} unresolved references to {args.save}')
        return
    if args.against:
        baseline = {line.strip() for line in open(args.against) if line.strip()}
        new = sorted(problems - baseline)
        for p in new:
            print(f'NEW  {p}')
        if new:
            sys.exit(f'{len(new)} new unresolved references')
        print(f'OK: no new unresolved references ({len(problems)} pre-existing, {len(baseline - problems)} fixed)')
        return
    for p in sorted(problems):
        print(p)
    print(f'{len(problems)} unresolved references')


if __name__ == '__main__':
    main()
