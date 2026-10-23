# potential.py <case>: per pattern, ideal sweep cost (startup + extent*count) vs current best
import sys
sys.path.insert(0, r'C:\Users\domok\OneDrive\Desktop\ccc 2026\level6')
from sim6 import parse_input
case = sys.argv[1]
L6 = r'C:\Users\domok\OneDrive\Desktop\ccc 2026\level6'
cases = parse_input(L6 + '\\in_level-6_%s.txt' % case)
best = {}
for line in open(L6 + '\\work\\best_lengths_%s.txt' % case):
    if line.startswith('#'): continue
    t = line.split()
    best[int(t[0])] = int(t[4])
tot_b = tot_i = 0
rows = []
for i, (r, g, cnt) in enumerate(cases):
    P = []
    for y, row in enumerate(g):
        for x, ch in enumerate(row):
            if ch == 'M':
                u, v, h = (x - y - r) // 2, (x + y - 3 * r) // 2, y - r
                P.append((h, u, v))
    bestc = 10**9
    for a in range(3):
        lo = min(p[a] for p in P); hi = max(p[a] for p in P)
        for s, f in ((lo, hi), (hi, lo)):
            c = abs(s) + 1 + abs(f - s)
            # subsequent builds: sweep back and forth: extent+1 each
            c += (cnt - 1) * (abs(f - s) + 1)
            bestc = min(bestc, c)
    rows.append((i, r, len(P), cnt, best[i], bestc))
    tot_b += best[i]; tot_i += min(bestc, best[i])
for x in sorted(rows, key=lambda x: x[5] - x[4]):
    print('%3d r=%2d P=%3d cnt=%d best=%3d ideal_sweep=%3d gain=%d' % (x + (x[4] - x[5],)))
print('total best', tot_b, 'if ideal sweep where better', tot_i)
