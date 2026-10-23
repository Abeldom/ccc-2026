# Level 3: several emitters; a beam passes over 'O' and 'M' cells and stops at anything else.
# Count how many beams cross each cell (same direction steps as level 2), then:
#   1 beam -> unchanged, exactly 2 beams -> 'M' (matter created), 3+ beams -> 'O' (destroyed).
import sys
D = {'1': (0, 2), '2': (1, 1), '3': (1, -1), '4': (0, -2), '5': (-1, -1), '6': (-1, 1)}
def solve(name):
    L = open('in_level-3_' + name + '.txt').read().split('\n')
    n = int(L[0]); p = 1; out = []
    for _ in range(n):
        while not L[p].strip(): p += 1          # skip blank separator lines
        rows = int(L[p].split()[1]); p += 1     # header: radius rows cols
        g = [list(x.rstrip('\r')) for x in L[p:p + rows]]; p += rows
        cnt = {}                                # (row, col) -> number of beams hitting it
        for r in range(rows):
            for c in range(len(g[r])):
                if g[r][c] in D:
                    dr, dc = D[g[r][c]]; y, x = r + dr, c + dc
                    while 0 <= y < rows and 0 <= x < len(g[y]) and g[y][x] in 'OM':
                        cnt[(y, x)] = cnt.get((y, x), 0) + 1; y += dr; x += dc
        for (y, x), k in cnt.items():           # apply counts after all beams are traced
            if k == 2: g[y][x] = 'M'
            elif k >= 3: g[y][x] = 'O'
        out.append('\n'.join(''.join(r) for r in g) + '\n')
    open('out_level-3_' + name + '.txt', 'w', newline='\n').write('\n'.join(out) + '\n')
for a in sys.argv[1:]: solve(a)
