# write cycles/pub_lengths_<case>.txt (format of best_lengths) from the min of published file and my candidate
import os, sys
L6 = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
def read_out(path):
    progs, cur = [], []
    for line in open(path).read().split('\n'):
        if not line.strip():
            if cur: progs.append(cur); cur = []
            continue
        cur.append(line)
    if cur: progs.append(cur)
    return progs
for case in sys.argv[1:]:
    toks = open(os.path.join(L6, 'in_level-6_%s.txt' % case)).read().split()
    n = int(toks[0]); pos = 1; info = []
    for i in range(n):
        r, rows, cols, c = map(int, toks[pos:pos + 4]); g = toks[pos + 4:pos + 4 + rows]; pos += 4 + rows
        info.append((r, sum(row.count('M') for row in g), c))
    pub = read_out(os.path.join(L6, 'out_level-6_%s.txt' % case))
    mine_p = os.path.join(L6, 'work', 'cand', 'cycles', 'out_level-6_%s.txt' % case)
    mine = read_out(mine_p) if os.path.exists(mine_p) else pub
    with open(os.path.join(os.path.dirname(__file__), 'pub_lengths_%s.txt' % case), 'w') as f:
        f.write('# pattern_index r |P| count best_lines\n')
        for i in range(n):
            f.write('%d %d %d %d %d\n' % (i, info[i][0], info[i][1], info[i][2], min(len(pub[i]), len(mine[i]))))
    print(case, sum(min(len(pub[i]), len(mine[i])) for i in range(n)))
