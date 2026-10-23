# Merge solution files ("# idx ..." headers followed by program lines) into a full output file.
# For each pattern keep the shortest program that the exact simulator (sim6.py) accepts; baseline fallback.
# Usage: python merge.py <case> <out file> <sol file> [<sol file> ...]
import os, sys
L6 = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, L6)
from sim6 import parse_input, geometry, simulate


def read_sol(path):
    res = {}
    cur = None
    for line in open(path):
        line = line.rstrip('\n')
        if line.startswith('#'):
            cur = int(line.split()[1])
            res.setdefault(cur, []).append([])
        elif line.strip() and cur is not None:
            res[cur][-1].append(line.split())
    return res


def read_out(path):
    progs, cur = [], []
    for line in open(path).read().split('\n'):
        if not line.strip():
            if cur: progs.append(cur); cur = []
            continue
        cur.append(line.split())
    if cur: progs.append(cur)
    return progs


def main(case, outp, sols):
    cases = parse_input(os.path.join(L6, 'in_level-6_%s.txt' % case))
    base = read_out(os.path.join(L6, 'work', 'cand', 'baseline', 'out_level-6_%s.txt' % case))
    best = list(base)
    src = ['base'] * len(cases)
    if os.path.exists(outp):
        old = read_out(outp)
        if len(old) == len(cases):
            for i, p in enumerate(old):
                if len(p) < len(best[i]): best[i] = p; src[i] = 'old'
    geos = {}
    cand = {}
    for sp in sols:
        for i, plist in read_sol(sp).items():
            for p in plist: cand.setdefault(i, []).append((len(p), sp, p))
    for i, (r, g, cnt) in enumerate(cases):
        if r not in geos: geos[r] = geometry(r)
        # verify current best (old file may contain anything)
        if src[i] == 'old' and simulate(r, g, best[i], geos[r]) < cnt:
            best[i] = base[i]; src[i] = 'base'
        for n, sp, p in sorted(cand.get(i, []), key=lambda x: x[0]):
            if n >= len(best[i]): break
            if all(len(l) == 6 and all(c in 'FBT-' for c in l) for l in p) and simulate(r, g, p, geos[r]) >= cnt:
                best[i] = p; src[i] = sp
                break
    open(outp, 'w', newline='\n').write('\n'.join('\n'.join(' '.join(l) for l in p) + '\n' for p in best))
    from collections import Counter
    print(case, 'total', sum(len(p) for p in best), Counter(os.path.basename(s) for s in src))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3:])
