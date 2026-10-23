# Compare a candidate output against the published best, per pattern.
# Usage: python cmp.py <case> <candidate file> [reference file]
import sys, os
L6 = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, L6)
from sim6 import parse_input, geometry, simulate

def read_programs(path):
    progs, cur = [], []
    for line in open(path).read().split('\n'):
        if line.strip() == '':
            if cur: progs.append(cur); cur = []
            continue
        cur.append(line.split())
    if cur: progs.append(cur)
    return progs

def main():
    case, cand = sys.argv[1], sys.argv[2]
    ref = sys.argv[3] if len(sys.argv) > 3 else os.path.join(L6, 'out_level-6_%s.txt' % case)
    cases = parse_input(os.path.join(L6, 'in_level-6_%s.txt' % case))
    A, B = read_programs(cand), read_programs(ref)
    geos = {}
    win = tie = lose = inval = 0; gain = 0; tot = 0
    wl = []
    for i, (r, g, cnt) in enumerate(cases):
        if r not in geos: geos[r] = geometry(r)
        ok = i < len(A) and simulate(r, g, A[i], geos[r]) >= cnt
        if not ok: inval += 1; tot += len(B[i]); continue
        a, b = len(A[i]), len(B[i])
        tot += min(a, b)
        if a < b: win += 1; gain += b - a; wl.append('%d:%d<%d' % (i, a, b))
        elif a == b: tie += 1
        else: lose += 1
    print('%s: wins %d (gain %d) ties %d loses %d invalid %d | merged total %d' % (case, win, gain, tie, lose, inval, tot))
    if wl: print(' '.join(wl))

main()
