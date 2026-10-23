# Build my candidate file: per pattern the shortest VALID program among
#   published best (L6/out_level-6_<case>.txt), my existing candidate file, and my sol files ("# idx" format).
# Usage: python mkcand.py <case> [sol files...]   (writes L6/work/cand/cycles/out_level-6_<case>.txt)
import os, sys, glob
L6 = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, L6)
from sim6 import parse_input, geometry, simulate


def read_out(path):
    progs, cur = [], []
    for line in open(path).read().split('\n'):
        if not line.strip():
            if cur: progs.append(cur); cur = []
            continue
        cur.append(line.split())
    if cur: progs.append(cur)
    return progs


def read_sol(path):
    res, cur = {}, None
    for line in open(path):
        line = line.rstrip('\n')
        if line.startswith('#'):
            cur = int(line.split()[1]); res.setdefault(cur, []).append([])
        elif line.strip() and cur is not None:
            res[cur][-1].append(line.split())
    return res


def wellformed(p):
    return len(p) > 0 and all(len(l) == 6 and all(c in ('F', 'B', 'T', '-') for c in l) for l in p)


def main(case, sols):
    cases = parse_input(os.path.join(L6, 'in_level-6_%s.txt' % case))
    n = len(cases)
    pub = read_out(os.path.join(L6, 'out_level-6_%s.txt' % case))
    assert len(pub) == n
    mine_path = os.path.join(L6, 'work', 'cand', 'cycles', 'out_level-6_%s.txt' % case)
    cands = [[(len(pub[i]), 'pub', pub[i])] for i in range(n)]
    if os.path.exists(mine_path):
        mine = read_out(mine_path)
        if len(mine) == n:
            for i in range(n): cands[i].append((len(mine[i]), 'mine', mine[i]))
    for sp in sols:
        for i, pl in read_sol(sp).items():
            for p in pl: cands[i].append((len(p), os.path.basename(sp), p))
    geos = {}
    out, wins, total, pubtot = [], [], 0, 0
    for i, (r, g, cnt) in enumerate(cases):
        if r not in geos: geos[r] = geometry(r)
        pubtot += len(pub[i])
        chosen = None
        for ln, src, p in sorted(cands[i], key=lambda x: (x[0], x[1] != 'pub')):
            if wellformed(p) and simulate(r, g, p, geos[r]) >= cnt:
                chosen = (ln, src, p); break
        assert chosen is not None
        out.append(chosen[2]); total += chosen[0]
        if chosen[0] < len(pub[i]): wins.append((i, len(pub[i]), chosen[0], chosen[1]))
    os.makedirs(os.path.dirname(mine_path), exist_ok=True)
    tmp = mine_path + '.tmp'
    open(tmp, 'w', newline='\n').write('\n'.join('\n'.join(' '.join(l) for l in p) + '\n' for p in out))
    os.replace(tmp, mine_path)
    print('%s: total %d (published %d), wins %d, saved %d' % (case, total, pubtot, len(wins), sum(a - b for _, a, b, _ in wins)))
    for w in wins: print('  win', *w)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2:])
