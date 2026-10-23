# per-pattern line counts of output files vs best_lengths
import sys, os, glob
L6 = r"C:\Users\domok\OneDrive\Desktop\ccc 2026\level6"
def lens(path):
    txt = open(path).read().strip('\n')
    progs = txt.split('\n\n')
    return [len([l for l in p.split('\n') if l.strip()]) for p in progs]
def best(case):
    rows = []
    for l in open(os.path.join(L6, 'work', 'best_lengths_%s.txt' % case)):
        if l.startswith('#') or not l.strip(): continue
        rows.append(list(map(int, l.split())))
    return rows
for case in ['2-medium', '3-large', '4-enormous']:
    b = best(case)
    print(case, 'best total', sum(x[4] for x in b))
    for d in sorted(glob.glob(os.path.join(L6, 'work', 'cand', '*'))):
        f = os.path.join(d, 'out_level-6_%s.txt' % case)
        if not os.path.exists(f): continue
        try: l = lens(f)
        except Exception as e: print(' ', d, 'err', e); continue
        if len(l) != len(b): print('  ', os.path.basename(d), 'count mismatch', len(l)); continue
        wins = sum(1 for x, y in zip(l, b) if x == y[4]); strict = 0
        print('   %-15s total %6d  ties/wins %3d' % (os.path.basename(d), sum(l), wins))
