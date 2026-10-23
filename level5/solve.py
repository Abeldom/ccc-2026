# Level 5: same emitter simulation as level 4, plus a pick-up arm that submits a pattern.
# Per tick: (1) run instruction line (F/B move along border, T toggle), (2) active beams
# create matter (exactly 2 beams) / destroy matter (3+ beams), (3) pattern check:
#   - every matter cell of the pattern must be matter on the board (exact position), and
#   - no OTHER matter on the board may be hex-adjacent to a pattern cell.
# If valid, the pattern's matter is removed and the counter increases. Output: count per program.
import sys

# beam direction of emitter 1..6 as (row, col) step on the character grid
DIRS = [(0, 2), (1, 1), (1, -1), (0, -2), (-1, -1), (-1, 1)]
# the 6 hex neighbours of a cell on the character grid (same steps, all directions)
NEIGH = [(0, 2), (0, -2), (1, 1), (1, -1), (-1, 1), (-1, -1)]

def border_cycle(r):
    """Border cells in clockwise order, starting at the left corner (emitter 1) going up."""
    path = []
    for k in range(r): path.append((r - k, k))              # left corner up to before top-left corner
    for k in range(r): path.append((0, r + 2 * k))          # top row, left to right
    for k in range(r): path.append((k, 3 * r + k))          # upper-right side down
    for k in range(r): path.append((r + k, 4 * r - k))      # lower-right side down
    for k in range(r): path.append((2 * r, 3 * r - 2 * k))  # bottom row, right to left
    for k in range(r): path.append((2 * r - k, r - k))      # lower-left side up
    return path

def simulate(r, pattern, prog):
    rows = 2 * r + 1
    interior, border = set(), set()
    for i in range(rows):
        d = abs(i - r)
        m = 2 * r - d
        for j in range(m + 1):
            (border if d == r or j == 0 or j == m else interior).add((i, d + 2 * j))
    path = border_cycle(r)
    assert len(path) == 6 * r == len(set(path)) and set(path) == border   # sanity check
    n = len(path)
    # cells forming the pattern; also the set of cells around it that must be free of other matter
    pat = {(y, x) for y, row in enumerate(pattern) for x, ch in enumerate(row) if ch == 'M'}
    assert pat <= interior
    ring = {(y + dy, x + dx) for (y, x) in pat for dy, dx in NEIGH} - pat
    pos = [e * r for e in range(6)]      # emitter e starts at its corner on the border cycle
    active = [False] * 6
    matter = set()
    count = 0
    for line in prog:
        for e, ins in enumerate(line):                      # 1. execute instructions
            if ins == 'F': pos[e] = (pos[e] + 1) % n
            elif ins == 'B': pos[e] = (pos[e] - 1) % n
            elif ins == 'T': active[e] = not active[e]
        cnt = {}
        for e in range(6):                                  # 2. beams create / destroy matter
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
        if pat and pat <= matter and not (ring & matter):   # 3. pattern detection
            matter -= pat                                   # submit: remove its matter
            count += 1
    return count

def solve(name):
    t = open('in_level-5_' + name + '.txt').read().split()
    n = int(t[0]); p = 1; out = []
    for _ in range(n):
        r, prow, _pcol = int(t[p]), int(t[p + 1]), int(t[p + 2]); p += 3
        pattern = t[p:p + prow]; p += prow                  # pattern lines contain no spaces
        l = int(t[p]); p += 1
        prog = [t[p + 6 * i: p + 6 * i + 6] for i in range(l)]; p += 6 * l
        out.append(str(simulate(r, pattern, prog)))
    assert p == len(t)                                      # all input consumed
    open('out_level-5_' + name + '.txt', 'w', newline='\n').write('\n'.join(out) + '\n')

for a in sys.argv[1:]: solve(a)
