# Level 6 baseline ("walker"): always valid, not optimal.
# Emitters 1 (horizontal line) and 2 (SE diagonal) are both active; their single crossing cell is
# the only cell with 2 hits, so each tick creates exactly one cell. Moving each emitter by at most
# one step per tick moves the crossing to any of the 6 neighbours (or 2 "jump" cells), so we walk a
# DFS tour over the pattern cells, never touching a non-pattern cell. Each finished tour submits the
# pattern once; the next build walks the same tour backwards.
# Usage: python baseline.py <input> <output>
import sys
from sim6 import parse_input

# cell (y, x) on the character grid -> line offsets: u = SE-diagonal, v = SW-diagonal, h = row
def coords(y, x, r):
    return (x - y - r) // 2, (x + y - 3 * r) // 2, y - r

# emitter 1 sits on row h at signed offset k1 = -h, emitter 2 on diagonal u at k2 = +u
def walker_tour(cells):
    """DFS over pattern cells where one tick can move (dh, du) by at most 1 each.
    Returns a walk (list of (h, u)) visiting every cell; the deepest branch is visited last."""
    cs = set(cells)
    start = min(cells, key=lambda c: abs(c[0]) + abs(c[1]))
    steps = [(a, b) for a in (-1, 0, 1) for b in (-1, 0, 1) if (a, b) != (0, 0)]
    nbrs = {c: [(c[0] + a, c[1] + b) for a, b in steps if (c[0] + a, c[1] + b) in cs] for c in cs}
    # build DFS tree
    parent, order, seen, stack = {start: None}, [], {start}, [start]
    children = {c: [] for c in cs}
    while stack:
        c = stack.pop()
        order.append(c)
        for n in nbrs[c]:
            if n not in seen:
                seen.add(n); parent[n] = c; children[c].append(n); stack.append(n)
    assert len(seen) == len(cs)
    depth = {}
    for c in reversed(order):                      # subtree height, to visit the deepest child last
        depth[c] = 1 + max((depth[k] for k in children[c]), default=0)
    walk = []
    def visit(c, last):
        walk.append(c)
        kids = sorted(children[c], key=lambda k: depth[k])
        for i, k in enumerate(kids):
            is_last = last and i == len(kids) - 1
            visit(k, is_last)
            if not is_last:
                walk.append(c)                     # walk back to c (re-creating existing matter is harmless)
    sys.setrecursionlimit(10000)
    visit(start, True)
    return walk

def solve_case(r, grid, count):
    cells = [(coords(y, x, r)[2], coords(y, x, r)[0])
             for y, row in enumerate(grid) for x, ch in enumerate(row) if ch == 'M']
    walk = walker_tour(cells)
    prog = []
    k1, k2 = 0, 0
    t1, t2 = -walk[0][0], walk[0][1]
    while (k1, k2) != (t1, t2):                    # move both emitters (still inactive) to the first cell
        line = ['-'] * 6
        if k1 != t1: d = 1 if t1 > k1 else -1; k1 += d; line[0] = 'F' if d > 0 else 'B'
        if k2 != t2: d = 1 if t2 > k2 else -1; k2 += d; line[1] = 'F' if d > 0 else 'B'
        prog.append(line)
    prog.append(['T', 'T', '-', '-', '-', '-'])    # switch both on: creates the first cell
    def step_to(c):
        nonlocal k1, k2
        line = ['-'] * 6
        d1, d2 = -c[0] - k1, c[1] - k2
        assert abs(d1) <= 1 and abs(d2) <= 1
        if d1: line[0] = 'F' if d1 > 0 else 'B'
        if d2: line[1] = 'F' if d2 > 0 else 'B'
        k1, k2 = -c[0], c[1]
        prog.append(line)
    for c in walk[1:]: step_to(c)                  # first build -> first submission
    for b in range(1, count):                      # later builds: walk the tour back and forth
        w = walk[::-1] if b % 2 == 1 else walk
        prog.append(['-'] * 6)                     # stay one tick to re-create the cell we stand on
        for c in w[1:]: step_to(c)
    return prog

def main(inp, outp):
    cases = parse_input(inp)
    out = []
    for r, g, cnt in cases:
        out.append('\n'.join(' '.join(l) for l in solve_case(r, g, cnt)) + '\n')
    open(outp, 'w', newline='\n').write('\n'.join(out))

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
