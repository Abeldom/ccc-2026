# trace.py <input> <output> <case> [tick_to_dump ...]
import sys
sys.path.insert(0, '../..')
from sim6 import parse_input, parse_output, geometry
inp, outp, ci = sys.argv[1], sys.argv[2], int(sys.argv[3])
dumps = set(int(x) for x in sys.argv[4:])
cases = parse_input(inp)
progs = parse_output(outp, len(cases))
r, g, cnt = cases[ci]
prog = progs[ci]
interior, border, per, ray = geometry(r)
n = len(per)
pat = {(y, x) for y, row in enumerate(g) for x, ch in enumerate(row) if ch == 'M'}
nb = [(0, 2), (0, -2), (1, 1), (1, -1), (-1, 1), (-1, -1)]
skirt = ({(y + dy, x + dx) for (y, x) in pat for dy, dx in nb} - pat) & interior
SG = [-1, 1, 1, 1, -1, -1]; AXN = 'HUVHUV'
pos = [e * r for e in range(6)]; act = [False] * 6; matter = set(); acc = 0
def off(e):
    k = (pos[e] - e * r) % n
    if k > n // 2: k -= n
    return k
for t, line in enumerate(prog, 1):
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
    miss = len(pat - matter); dirty = len(skirt & matter)
    em = ' '.join('%s%+d%s' % (AXN[e], SG[e] * off(e), '*' if act[e] else ' ') for e in range(6))
    accd = ''
    if pat <= matter and not (skirt & matter):
        matter -= pat; acc += 1; accd = ' ACCEPT'
    print('%3d %s  %s  miss %3d dirty %3d%s' % (t, ' '.join(line), em, miss, dirty, accd))
    if t in dumps:
        grid = [list(row) for row in g]
        for (y, x) in interior:
            P = (y, x) in pat; K = (y, x) in skirt; m = (y, x) in matter
            grid[y][x] = ('M' if m else 'm') if P else (('k' if m else 'o') if K else ('#' if m else '.'))
        for row in grid: print('      ' + ''.join(row))
