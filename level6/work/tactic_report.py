# For each counted case: which tactic (solver folder) provides the shortest valid program per pattern.
import os, sys, glob, collections
L6 = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, L6); os.chdir(L6)
from sim6 import parse_input, geometry, simulate
from combine import read_programs

def label(path):
    parts = os.path.relpath(path, 'work').replace(os.sep, '/').split('/')
    if parts[0] == 'cand': return parts[1].replace('x_', '').split('_')[0]
    return parts[0]

srcs = [p for p in glob.glob('work/**/*.txt', recursive=True)
        if '_harvest' not in p and 'tmp_best' not in p and os.path.getsize(p) >= 1000]
for case in ['2-medium', '3-large', '4-enormous']:
    cases = parse_input('in_level-6_%s.txt' % case)
    files = []
    for p in srcs:
        try: progs = read_programs(p)
        except Exception: continue
        if len(progs) == len(cases): files.append((label(p), progs))
    geos, best = {}, []
    for i, (r, g, cnt) in enumerate(cases):
        if r not in geos: geos[r] = geometry(r)
        opts = sorted(((len(pr[i]), lab, k) for k, (lab, pr) in enumerate(files) if pr[i] is not None))
        for n, lab, k in opts:
            if simulate(r, g, files[k][1][i], geos[r]) >= cnt:
                ties = sorted({l for m, l, _ in opts if m == n}); best.append((n, ties)); break
    wins = collections.Counter()
    for n, ties in best:
        for t in ties: wins[t] += 1
    print('%s: best-of-all = %d lines; patterns won per tactic (ties count for each): %s'
          % (case, sum(n for n, _ in best), dict(wins.most_common())))
