# Visualize a program for one pattern: python viz.py <input> <output> <pattern index>
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from sim6 import parse_input, parse_output, geometry

def main(inp, outp, idx):
    cases = parse_input(inp)
    progs = parse_output(outp, len(cases))
    r, g, cnt = cases[idx]
    prog = progs[idx]
    interior, border, per, ray = geometry(r)
    n = len(per)
    pat = {(y, x) for y, row in enumerate(g) for x, ch in enumerate(row) if ch == 'M'}
    nb = [(0, 2), (0, -2), (1, 1), (1, -1), (-1, 1), (-1, -1)]
    skirt = {(y + dy, x + dx) for (y, x) in pat for dy, dx in nb} - pat
    skirt &= interior
    pos = [e * r for e in range(6)]
    act = [False] * 6
    matter = set()
    acc = 0
    print('r=%d |P|=%d cnt=%d len=%d' % (r, len(pat), cnt, len(prog)))
    for t, line in enumerate(prog):
        for e, ins in enumerate(line):
            if ins == 'F': pos[e] = (pos[e] + 1) % n
            elif ins == 'B': pos[e] = (pos[e] - 1) % n
            elif ins == 'T': act[e] = not act[e]
        hit = {}
        for e in range(6):
            if act[e]:
                for c in ray[pos[e]][e]:
                    hit[c] = hit.get(c, 0) + 1
        created = set(); destroyed = set()
        for c, k in hit.items():
            if k == 2:
                if c not in matter: created.add(c)
                matter.add(c)
            elif k >= 3:
                if c in matter: destroyed.add(c)
                matter.discard(c)
        ok = pat <= matter and not (skirt & matter)
        offs = []
        for e in range(6):
            k = (pos[e] - e * r) % n
            if k > n // 2: k -= n
            offs.append('%+d%s' % (k, '*' if act[e] else ' '))
        print('tick %d: %s   %s  newP=%d newK=%d destroyed=%d%s' % (t + 1, ' '.join(line), ' '.join(offs),
              len(created & pat), len(created & skirt), len(destroyed & (pat | skirt)), '  ACCEPT' if ok else ''))
        rows = []
        for y in range(2 * r + 1):
            s = ''
            for x in range(4 * r + 1):
                c = (y, x)
                if c in pat: ch = ('#' if c in created else 'M') if c in matter else ('x' if c in destroyed else 'm')
                elif c in skirt: ch = ('!' if c in created else 'K') if c in matter else ('o' if c in destroyed else '.')
                elif c in interior: ch = ('+' if c in matter else ' ')
                elif c in border: ch = '*'
                else: ch = ' '
                s += ch
            rows.append(s)
        print('\n'.join(rows))
        if ok:
            matter -= pat; acc += 1
    print('accepted', acc)

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]))
