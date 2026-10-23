# Collect RES blocks from work/sparse/res/*.txt (case inferred from file name: contains 'medium','large','enormous'),
# verify each program exactly, keep the shortest per pattern in work/sparse/mine_<case>.txt, and write the candidate
# file work/cand/sparse/out_level-6_<case>.txt = per pattern min(mine, published). Validates the result.
import sys, os, glob
L6 = 'C:/Users/domok/OneDrive/Desktop/ccc 2026/level6'
sys.path.insert(0, L6)
from sim6 import parse_input, parse_output, geometry, simulate

CASES = {'medium': '2-medium', 'large': '3-large', 'enormous': '4-enormous'}
W = L6 + '/work/sparse'

def read_res(path):
    out = []
    lines = open(path).read().split('\n')
    i = 0
    while i < len(lines):
        t = lines[i].split()
        if len(t) == 3 and t[0] == 'RES':
            idx, n = int(t[1]), int(t[2])
            prog = [lines[i + 1 + j].split() for j in range(n)]
            if all(len(x) == 6 for x in prog): out.append((idx, prog))
            i += 1 + n
        else:
            i += 1
    return out

def load_mine(case):
    p = W + '/mine_%s.txt' % case
    d = {}
    if os.path.exists(p):
        for idx, prog in read_res(p): d[idx] = prog
    return d

def save_mine(case, d):
    with open(W + '/mine_%s.txt' % case, 'w', newline='\n') as f:
        for idx in sorted(d):
            f.write('RES %d %d\n' % (idx, len(d[idx])))
            for l in d[idx]: f.write(' '.join(l) + '\n')

def main():
    for short, case in CASES.items():
        cases = parse_input(L6 + '/in_level-6_%s.txt' % case)
        geos = {}
        def ok(i, prog):
            r, g, cnt = cases[i]
            if r not in geos: geos[r] = geometry(r)
            return simulate(r, g, prog, geos[r]) >= cnt
        mine = load_mine(case)
        changed = False
        for f in glob.glob(W + '/res/*%s*.txt' % short):
            for idx, prog in read_res(f):
                if idx in mine and len(mine[idx]) <= len(prog): continue
                if ok(idx, prog):
                    mine[idx] = prog; changed = True
                else:
                    print('INVALID', f, idx)
        if changed: save_mine(case, mine)
        pub = parse_output(L6 + '/out_level-6_%s.txt' % case, len(cases))
        final = []; nb = 0; gain = 0
        for i in range(len(cases)):
            if i in mine and len(mine[i]) < len(pub[i]):
                final.append(mine[i]); nb += 1; gain += len(pub[i]) - len(mine[i])
            else:
                final.append(pub[i])
        os.makedirs(L6 + '/work/cand/sparse', exist_ok=True)
        tmp = L6 + '/work/cand/sparse/tmp_%s.txt' % case
        with open(tmp, 'w', newline='\n') as f:
            f.write('\n\n'.join('\n'.join(' '.join(l) for l in p) for p in final) + '\n')
        # validate whole file
        progs = parse_output(tmp, len(cases))
        bad = sum(1 for i in range(len(cases)) if not ok(i, progs[i]))
        tot = sum(len(p) for p in progs)
        if bad == 0 and len(progs) == len(cases):
            os.replace(tmp, L6 + '/work/cand/sparse/out_level-6_%s.txt' % case)
        print('%s: mine %d, better than published %d (gain %d), total %d, invalid %d' % (case, len(mine), nb, gain, tot, bad))

if __name__ == '__main__':
    main()
