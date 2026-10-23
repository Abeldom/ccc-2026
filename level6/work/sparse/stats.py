import sys
L6 = 'C:/Users/domok/OneDrive/Desktop/ccc 2026/level6'
sys.path.insert(0, L6)
from sim6 import parse_input, parse_output, geometry

def acc_ticks(r, pattern, prog, geo):
    interior, border, per, ray = geo
    n = len(per)
    pat = {(y, x) for y, row in enumerate(pattern) for x, ch in enumerate(row) if ch == 'M'}
    nb = [(0, 2), (0, -2), (1, 1), (1, -1), (-1, 1), (-1, -1)]
    skirt = ({(y + dy, x + dx) for (y, x) in pat for dy, dx in nb} - pat) & interior
    pos = [e * r for e in range(6)]; act = [False] * 6; matter = set(); res = []
    for t, line in enumerate(prog):
        for e, ins in enumerate(line):
            if ins == 'F': pos[e] = (pos[e] + 1) % n
            elif ins == 'B': pos[e] = (pos[e] - 1) % n
            elif ins == 'T': act[e] = not act[e]
        hit = {}
        for e in range(6):
            if act[e]:
                for c in ray[pos[e]][e]: hit[c] = hit.get(c, 0) + 1
        for c, k in hit.items():
            if k == 2: matter.add(c)
            elif k >= 3: matter.discard(c)
        if pat <= matter and not (skirt & matter):
            matter -= pat; res.append(t + 1)
    return res, pat, skirt

maxP = int(sys.argv[1]) if len(sys.argv) > 1 else 30
for case in ['2-medium', '3-large', '4-enormous']:
    cases = parse_input(L6 + '/in_level-6_%s.txt' % case)
    progs = parse_output(L6 + '/out_level-6_%s.txt' % case, len(cases))
    geos = {}
    tot = 0; totlb = 0; n = 0; sumc = 0
    rows = []
    for i, ((r, g, cnt), p) in enumerate(zip(cases, progs)):
        if r not in geos: geos[r] = geometry(r)
        a, pat, sk = acc_ticks(r, g, p, geos[r])
        if len(pat) > maxP: continue
        lb1 = 0
        for (y, x) in pat:
            u = (x - y - r) // 2; v = (x + y - 3 * r) // 2; h = y - r
            lb1 = max(lb1, 1 + min(abs(u), abs(v), abs(h)))
        a = a[:cnt]
        d = [a[0]] + [a[j] - a[j - 1] for j in range(1, len(a))]
        lb = lb1 + (cnt - 1)
        tot += len(p); totlb += lb; n += 1; sumc += cnt
        rows.append((len(p) - lb, i, r, len(pat), len(sk), cnt, len(p), lb1, d))
    rows.sort(reverse=True)
    print('==', case, 'n', n, 'sumcount', sumc, 'lines', tot, 'weakLB', totlb)
    for x in rows[:int(sys.argv[2]) if len(sys.argv) > 2 else 40]:
        print('gap %3d  i %3d r %d P %3d K %3d cnt %2d best %3d lb1 %d builds %s' % (x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8][:8]))
