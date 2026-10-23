# Level 4: simulate 6 emitters on the border of a hex board, driven by a program.
# Character grid: 2r+1 rows, 4r+1 cols; row i (d=|i-r|) has cells at cols d, d+2, ..., d+2(2r-d).
# Border cells: first/last cell of each row, plus the whole top and bottom row. The rest are interior.
# Each tick: (1) every emitter executes its instruction (F=clockwise step, B=counterclockwise
# step, T=toggle on/off, -=nothing); (2) every ACTIVE emitter fires a beam through interior
# cells along its direction; (3) cells hit by exactly 2 beams become matter 'M',
# cells hit by 3+ beams become empty 'O', cells hit by 1 beam stay unchanged (as in level 3).
import sys

# direction of each emitter 1..6 as (row, col) step: right, down-right, down-left, left, up-left, up-right
DIRS = [(0, 2), (1, 1), (1, -1), (0, -2), (-1, -1), (-1, 1)]

def border_cycle(r):
    """Border cells in clockwise order, starting at the left corner (emitter 1) going up."""
    path = []
    for k in range(r):            # left corner -> up the upper-left side (excluding top-left corner)
        path.append((r - k, k))
    for k in range(r):            # top-left corner -> along the top row (excluding top-right corner)
        path.append((0, r + 2 * k))
    for k in range(r):            # top-right corner -> down the upper-right side
        path.append((k, 3 * r + k))
    for k in range(r):            # right corner -> down the lower-right side
        path.append((r + k, 4 * r - k))
    for k in range(r):            # bottom-right corner -> along the bottom row leftwards
        path.append((2 * r, 3 * r - 2 * k))
    for k in range(r):            # bottom-left corner -> up the lower-left side
        path.append((2 * r - k, r - k))
    return path

def simulate(r, prog):
    rows, cols = 2 * r + 1, 4 * r + 1
    interior, border = set(), set()
    for i in range(rows):
        d = abs(i - r)
        m = 2 * r - d             # last cell index in this row
        for j in range(m + 1):
            c = d + 2 * j
            (border if d == r or j == 0 or j == m else interior).add((i, c))
    path = border_cycle(r)
    assert len(path) == 6 * r == len(set(path)) and set(path) == border   # sanity check
    n = len(path)
    # emitter e starts at corner index: e1=0, e2=r, e3=2r, e4=3r, e5=4r, e6=5r
    pos = [e * r for e in range(6)]
    active = [False] * 6
    matter = set()
    for line in prog:
        for e, ins in enumerate(line):
            if ins == 'F': pos[e] = (pos[e] + 1) % n
            elif ins == 'B': pos[e] = (pos[e] - 1) % n
            elif ins == 'T': active[e] = not active[e]
        cnt = {}
        for e in range(6):
            if not active[e]: continue
            dr, dc = DIRS[e]
            y, x = path[pos[e]]
            y += dr; x += dc
            while (y, x) in interior:
                cnt[(y, x)] = cnt.get((y, x), 0) + 1
                y += dr; x += dc
        for cell, k in cnt.items():
            if k == 2: matter.add(cell)
            elif k >= 3: matter.discard(cell)
    # render
    g = [['-'] * cols for _ in range(rows)]
    for (y, x) in border: g[y][x] = '*'
    for (y, x) in interior: g[y][x] = 'M' if (y, x) in matter else 'O'
    return '\n'.join(''.join(row) for row in g) + '\n'

def solve(name):
    t = open('in_level-4_' + name + '.txt').read().split()
    n = int(t[0]); p = 1; out = []
    for _ in range(n):
        r, l = int(t[p]), int(t[p + 1]); p += 2
        prog = [t[p + 6 * i: p + 6 * i + 6] for i in range(l)]; p += 6 * l
        out.append(simulate(r, prog))
    assert p == len(t)            # all input consumed
    open('out_level-4_' + name + '.txt', 'w', newline='\n').write('\n'.join(out) + '\n')

for a in sys.argv[1:]: solve(a)
