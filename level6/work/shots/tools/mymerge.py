# Merge: per pattern shortest VALID program among the published best + given candidate files
# -> writes work/cand/shots/out_level-6_<case>.txt (only if all valid) and reports wins vs published.
# usage: python mymerge.py <case> file1 [file2 ...]
import sys, os
L6 = r"C:\Users\domok\OneDrive\Desktop\ccc 2026\level6"
sys.path.insert(0, L6)
from sim6 import parse_input, parse_output, geometry, simulate

def main():
    case = sys.argv[1]
    files = sys.argv[2:]
    cases = parse_input(os.path.join(L6, 'in_level-6_%s.txt' % case))
    n = len(cases)
    pub = os.path.join(L6, 'out_level-6_%s.txt' % case)
    mine = os.path.join(L6, 'work', 'cand', 'shots', 'out_level-6_%s.txt' % case)
    srcs = [pub] + ([mine] if os.path.exists(mine) else []) + files
    geos = {}
    best = [None] * n
    pubLen = [None] * n
    for si, f in enumerate(srcs):
        try:
            progs = parse_output(f, n)
        except Exception as e:
            print('skip', f, e); continue
        if len(progs) != n: print('skip', f, 'count', len(progs)); continue
        for i, ((r, g, cnt), p) in enumerate(zip(cases, progs)):
            if si == 0: pubLen[i] = len(p)
            if best[i] is not None and len(p) >= len(best[i]): continue
            if r not in geos: geos[r] = geometry(r)
            if simulate(r, g, p, geos[r]) >= cnt: best[i] = p
    assert all(b is not None for b in best)
    tot = sum(len(b) for b in best)
    wins = sum(1 for i in range(n) if len(best[i]) < pubLen[i])
    gain = sum(pubLen[i] - len(best[i]) for i in range(n))
    tmp = mine + '.tmp'
    with open(tmp, 'w', newline='\n') as fo:
        fo.write('\n\n'.join('\n'.join(' '.join(l) for l in p) for p in best) + '\n')
    os.replace(tmp, mine)
    print('%s: total %d (published %d), wins %d, gain %d' % (case, tot, sum(pubLen), wins, gain))

if __name__ == '__main__':
    main()
