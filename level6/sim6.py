# Level 6 validator: exact port of the simulation in sandbox-level-6.html.
# Usage: python sim6.py <input file> <output file>
#   Simulates every program against its pattern, checks accepted >= required count,
#   and prints the score (total number of instruction lines).
#
# Geometry (same as levels 4/5): board radius r, character grid of 2r+1 rows x 4r+1 cols.
# Border = first/last cell of each row + whole top/bottom row; the rest are interior cells.
# Perimeter = 6r border cells, clockwise from the left corner going up; emitter e (0-based)
# starts at index e*r, inactive. Emitter e fires in direction DIRS[e].
# Per tick: (1) apply instructions, (2) count beam hits on interior cells: exactly 2 -> matter,
# 3+ -> destroyed, (3) if every pattern 'M' is matter and no pattern 'O' that touches a pattern
# 'M' is matter, the pattern is accepted and its matter removed (at most once per tick).
import sys

DIRS = [(0, 2), (1, 1), (1, -1), (0, -2), (-1, -1), (-1, 1)]   # E, SE, SW, W, NW, NE as (row, col)


def geometry(r):
    """Return (interior set, border set, clockwise perimeter list, ray table).
    ray[p][e] = list of interior cells hit by emitter e standing on perimeter index p."""
    interior, border = set(), set()
    for i in range(2 * r + 1):
        d = abs(i - r)
        m = 2 * r - d
        for j in range(m + 1):
            (border if d == r or j == 0 or j == m else interior).add((i, d + 2 * j))
    per = []
    for k in range(r): per.append((r - k, k))               # upper-left side, going up
    for k in range(r): per.append((0, r + 2 * k))           # top row, going right
    for k in range(r): per.append((k, 3 * r + k))           # upper-right side, going down
    for k in range(r): per.append((r + k, 4 * r - k))       # lower-right side, going down
    for k in range(r): per.append((2 * r, 3 * r - 2 * k))   # bottom row, going left
    for k in range(r): per.append((2 * r - k, r - k))       # lower-left side, going up
    assert len(per) == 6 * r == len(set(per)) and set(per) == border
    allc = interior | border
    ray = []
    for (y, x) in per:
        row = []
        for dy, dx in DIRS:
            cells, cy, cx = [], y + dy, x + dx
            while (cy, cx) in allc:                          # like the sandbox: walk while inside the hexagon
                if (cy, cx) in interior: cells.append((cy, cx))
                cy += dy; cx += dx
            row.append(cells)
        ray.append(row)
    return interior, border, per, ray


def parse_input(path):
    """List of (r, pattern_rows(list of str), count)."""
    t = open(path).read().split()
    n = int(t[0]); p = 1; cases = []
    for _ in range(n):
        r, rows, cols, cnt = map(int, t[p:p + 4]); p += 4
        g = t[p:p + rows]; p += rows
        assert all(len(s) == cols for s in g)
        cases.append((r, g, cnt))
    assert p == len(t)
    return cases


def parse_output(path, n):
    """Split output into n programs (lists of 6-token lists), separated by empty lines."""
    progs, cur = [], []
    for line in open(path).read().split('\n'):
        if line.strip() == '':
            if cur: progs.append(cur); cur = []
            continue
        toks = line.split()
        assert len(toks) == 6 and all(c in 'FBT-' for c in toks), 'bad line: %r' % line
        cur.append(toks)
    if cur: progs.append(cur)
    return progs


def simulate(r, pattern, prog, geo=None):
    """Return number of accepted patterns after running prog."""
    interior, border, per, ray = geo or geometry(r)
    n = len(per)
    pat = {(y, x) for y, row in enumerate(pattern) for x, ch in enumerate(row) if ch == 'M'}
    nb = [(0, 2), (0, -2), (1, 1), (1, -1), (-1, 1), (-1, -1)]
    skirt = {(y + dy, x + dx) for (y, x) in pat for dy, dx in nb} - pat
    skirt &= interior                          # pattern 'O' cells touching 'M' (border cells never hold matter)
    pos = [e * r for e in range(6)]
    act = [False] * 6
    matter = set()
    acc = 0
    for line in prog:
        for e, ins in enumerate(line):
            if ins == 'F': pos[e] = (pos[e] + 1) % n
            elif ins == 'B': pos[e] = (pos[e] - 1) % n
            elif ins == 'T': act[e] = not act[e]
        hit = {}
        for e in range(6):
            if act[e]:
                for c in ray[pos[e]][e]:
                    hit[c] = hit.get(c, 0) + 1
        for c, k in hit.items():
            if k == 2: matter.add(c)
            elif k >= 3: matter.discard(c)
        if pat <= matter and not (skirt & matter):
            matter -= pat
            acc += 1
    return acc


def main(inp, outp):
    cases = parse_input(inp)
    progs = parse_output(outp, len(cases))
    assert len(progs) == len(cases), 'expected %d programs, got %d' % (len(cases), len(progs))
    total, bad = 0, 0
    geos = {}
    for i, ((r, g, cnt), prog) in enumerate(zip(cases, progs)):
        if r not in geos: geos[r] = geometry(r)
        a = simulate(r, g, prog, geos[r])
        total += len(prog)
        if a < cnt:
            bad += 1
            print('case %d: accepted %d < %d' % (i, a, cnt))
    print('%s: %d programs, %d invalid, score (total lines) = %d' % (outp, len(cases), bad, total))
    return bad == 0


if __name__ == '__main__':
    sys.exit(0 if main(sys.argv[1], sys.argv[2]) else 1)
