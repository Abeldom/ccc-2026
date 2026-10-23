# Combine candidate solutions: for every pattern, re-validate each candidate program with the exact
# simulator (sim6.py) and keep the shortest valid one.
# Usage: python combine.py <case name, e.g. 2-medium> <output file> <candidate dir> [<candidate dir> ...]
#   Each candidate dir may contain out_level-6_<case>.txt; programs that are missing, malformed or
#   invalid are ignored for that pattern.
import os, sys
from collections import Counter
from sim6 import parse_input, geometry, simulate

def read_programs(path):
    """Split a file into programs (lists of 6-token lists); None for unparsable programs."""
    progs, cur, ok = [], [], True
    for line in open(path).read().split('\n'):
        if line.strip() == '':
            if cur: progs.append(cur if ok else None); cur, ok = [], True
            continue
        toks = line.split()
        if len(toks) != 6 or any(t not in ('F', 'B', 'T', '-') for t in toks): ok = False
        cur.append(toks)
    if cur: progs.append(cur if ok else None)
    return progs

def main(case, outp, dirs):
    cases = parse_input('in_level-6_%s.txt' % case)
    cands = {}
    for d in dirs:
        f = os.path.join(d, 'out_level-6_%s.txt' % case)
        if os.path.exists(f):
            progs = read_programs(f)
            if len(progs) == len(cases): cands[d] = progs
            else: print('skip %s: %d programs, expected %d' % (f, len(progs), len(cases)))
    geos, best, wins, total = {}, [], Counter(), 0
    for i, (r, g, cnt) in enumerate(cases):
        if r not in geos: geos[r] = geometry(r)
        options = sorted((len(p[i]), d) for d, p in cands.items() if p[i] is not None)
        chosen = None
        for n, d in options:                      # shortest first; first valid wins
            if simulate(r, g, cands[d][i], geos[r]) >= cnt:
                chosen = (n, d); break
        assert chosen, 'no valid program for case %d' % i
        best.append(cands[chosen[1]][i]); wins[chosen[1]] += 1; total += chosen[0]
    open(outp, 'w', newline='\n').write('\n'.join('\n'.join(' '.join(l) for l in p) + '\n' for p in best))
    print('%s: score %d; wins %s' % (case, total, dict(wins)))

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3:])
