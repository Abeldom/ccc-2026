# Level 2: shoot a laser from the emitter ('1'-'6') and turn every 'O' on its path into 'M'.
# On the character grid each hex neighbour is a step of (row, col):
# 1 = right (0,+2), 4 = left (0,-2), 2 = down-right, 3 = down-left, 5 = up-left, 6 = up-right.
# The beam stops at the first non-'O' character (border '*', edge of grid).
import sys
D = {'1': (0, 2), '2': (1, 1), '3': (1, -1), '4': (0, -2), '5': (-1, -1), '6': (-1, 1)}
def solve(name):
    L = open('in_level-2_' + name + '.txt').read().split('\n')
    n = int(L[0]); p = 1; out = []
    for _ in range(n):
        while not L[p].strip(): p += 1          # skip blank separator lines
        rows = int(L[p].split()[1]); p += 1     # header: radius rows cols
        g = [list(x.rstrip('\r')) for x in L[p:p + rows]]; p += rows
        for r in range(rows):
            for c in range(len(g[r])):
                if g[r][c] in D:                # found an emitter: walk along its direction
                    dr, dc = D[g[r][c]]; y, x = r + dr, c + dc
                    while 0 <= y < rows and 0 <= x < len(g[y]) and g[y][x] == 'O':
                        g[y][x] = 'M'; y += dr; x += dc
        out.append('\n'.join(''.join(r) for r in g) + '\n')
    open('out_level-2_' + name + '.txt', 'w', newline='\n').write('\n'.join(out) + '\n')
for a in sys.argv[1:]: solve(a)
