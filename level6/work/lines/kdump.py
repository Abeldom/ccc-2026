import sys
sys.path.insert(0, r'C:\Users\domok\OneDrive\Desktop\ccc 2026\level6')
from sim6 import parse_input
from stats import cells, NB

def main(case, idx, axis):
    cases = parse_input(r'C:\Users\domok\OneDrive\Desktop\ccc 2026\level6\in_level-6_%s.txt' % case)
    r, g, cnt = cases[idx]
    P, I = cells(r, g)
    K = {(u + a, v + b) for (u, v) in P for a, b in NB} & I - P
    print('\n'.join(g))
    # coords (h,u,v)
    def co(c): return (c[1] - c[0], c[0], c[1])
    A = 'HUV'.index(axis)
    O = [i for i in range(3) if i != A]
    lines = sorted({co(c)[A] for c in P})
    for l in lines:
        ks = sorted(co(c) for c in K if co(c)[A] == l)
        ps = sorted(co(c) for c in P if co(c)[A] == l)
        print('%s=%+d  P%2d  K%d: %s' % (axis, l, len(ps), len(ks), ' '.join('(%s%+d,%s%+d)' % ('HUV'[O[0]], k[O[0]], 'HUV'[O[1]], k[O[1]]) for k in ks)))

if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]), sys.argv[3])
