import sys
sys.path.insert(0, r'C:\Users\domok\OneDrive\Desktop\ccc 2026\level6')
from sim6 import parse_input, parse_output

def cells(r, g):
    P = set(); I = set()
    for y, row in enumerate(g):
        for x, ch in enumerate(row):
            if ch in 'MO':
                u = (x - y - r) // 2; v = (x + y - 3 * r) // 2
                I.add((u, v))
                if ch == 'M': P.add((u, v))
    return P, I

NB = [(1, 1), (-1, -1), (0, 1), (0, -1), (1, 0), (-1, 0)]

def main(case):
    cases = parse_input(r'C:\Users\domok\OneDrive\Desktop\ccc 2026\level6\in_level-6_%s.txt' % case)
    progs = parse_output(r'C:\Users\domok\OneDrive\Desktop\ccc 2026\level6\out_level-6_%s.txt' % case, len(cases))
    for i, ((r, g, cnt), p) in enumerate(zip(cases, progs)):
        P, I = cells(r, g)
        K = {(u + a, v + b) for (u, v) in P for a, b in NB} & I - P
        out = []
        for ax, f in (('H', lambda c: c[1] - c[0]), ('U', lambda c: c[0]), ('V', lambda c: c[1])):
            lp = {}; lk = {}
            for c in P: lp[f(c)] = lp.get(f(c), 0) + 1
            for c in K: lk[f(c)] = lk.get(f(c), 0) + 1
            ls = sorted(lp)
            ks = [lk.get(l, 0) for l in ls]
            span = max(ls) - min(ls) + 1
            out.append('%s:%2d/%2d k%3d mx%d' % (ax, len(ls), span, sum(ks), max(ks)))
        print('%3d r=%2d c=%d |P|=%3d |K|=%3d best=%3d  %s' % (i, r, cnt, len(P), len(K), len(p), '  '.join(out)))

if __name__ == "__main__": main(sys.argv[1])
