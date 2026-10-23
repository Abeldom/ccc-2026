# Merge program files: shortest valid program per pattern; always includes the published best.
# Usage: python merge.py <case> <out file> <file> [<file> ...]
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

def valid_prog(p):
    return all(len(l) == 6 and all(c in 'FBT-' for c in l) for l in p)

def main():
    case, outp = sys.argv[1], sys.argv[2]
    files = [os.path.join(L6, 'out_level-6_%s.txt' % case)] + sys.argv[3:]
    cases = parse_input(os.path.join(L6, 'in_level-6_%s.txt' % case))
    P = []
    for f in files:
        if not os.path.exists(f): continue
        try: pr = read_programs(f)
        except Exception: continue
        if len(pr) == len(cases): P.append((f, pr))
        else: print('skip', f, len(pr))
    geos = {}; best = []; tot = 0; wins = {}
    for i, (r, g, cnt) in enumerate(cases):
        if r not in geos: geos[r] = geometry(r)
        opts = sorted(((len(pr[i]), k) for k, (f, pr) in enumerate(P)), key=lambda x: x)
        ch = None
        for n, k in opts:
            p = P[k][1][i]
            if valid_prog(p) and simulate(r, g, p, geos[r]) >= cnt: ch = (n, k); break
        assert ch is not None
        best.append(P[ch[1]][1][i]); tot += ch[0]
        if ch[1] != 0 and ch[0] < len(P[0][1][i]): wins[P[ch[1]][0]] = wins.get(P[ch[1]][0], 0) + 1
    tmp = outp + '.tmp'
    open(tmp, 'w', newline='\n').write('\n'.join('\n'.join(' '.join(l) for l in p) + '\n' for p in best))
    os.replace(tmp, outp)
    nw = sum(1 for i in range(len(cases)) if len(best[i]) < len(P[0][1][i]))
    print('%s: total %d, beats published on %d patterns %s' % (case, tot, nw, wins))

main()
