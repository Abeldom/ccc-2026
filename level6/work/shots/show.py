# show a program: per tick the active lines and effect on P/K
import sys
sys.path.insert(0, r'C:\Users\domok\OneDrive\Desktop\ccc 2026\level6')
from sim6 import parse_input, parse_output, geometry

SGN = [-1, 1, 1, 1, -1, -1]
AX = 'HUVHUV'

def main(inp, outp, idx):
    cases = parse_input(inp)
    progs = parse_output(outp, len(cases))
    r, g, cnt = cases[idx]
    prog = progs[idx]
    interior, border, per, ray = geometry(r)
    pat = {(y, x) for y, row in enumerate(g) for x, ch in enumerate(row) if ch == 'M'}
    nb = [(0, 2), (0, -2), (1, 1), (1, -1), (-1, 1), (-1, -1)]
    skirt = ({(y + dy, x + dx) for (y, x) in pat for dy, dx in nb} - pat) & interior
    print('r=%d cnt=%d |P|=%d |K|=%d len=%d' % (r, cnt, len(pat), len(skirt), len(prog)))
    print('\n'.join(g))
    n = len(per)
    k = [0] * 6; act = [False] * 6; matter = set(); acc = 0
    for t, line in enumerate(prog):
        for e, ins in enumerate(line):
            if ins == 'F': k[e] += 1
            elif ins == 'B': k[e] -= 1
            elif ins == 'T': act[e] = not act[e]
        hit = {}
        for e in range(6):
            if act[e]:
                for c in ray[(e * r + k[e]) % n][e]: hit[c] = hit.get(c, 0) + 1
        S = {c for c, v in hit.items() if v == 2}; D = {c for c, v in hit.items() if v >= 3}
        for c in S: matter.add(c)
        for c in D: matter.discard(c)
        lines = []
        for e in range(6):
            if act[e] and abs(k[e]) < r: lines.append('%s%+d' % (AX[e], SGN[e] * k[e]))
            elif act[e]: lines.append('%s(dead)' % AX[e])
        okp = len(pat & matter); badk = len(skirt & matter)
        accd = pat <= matter and not (skirt & matter)
        if accd: matter -= pat; acc += 1
        print('%3d %s  %-34s S:P%2d K%2d D:P%2d K%2d | have P %d/%d dirtyK %d %s' % (t + 1, ' '.join(line), ' '.join(lines), len(S & pat), len(S & skirt), len(D & pat), len(D & skirt), okp, len(pat), badk, 'ACCEPT %d' % acc if accd else ''))

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]))
