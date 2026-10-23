# Harvest every complete output file that the solver agents left anywhere under work/ (intermediate
# results included) into work/cand/_harvest/<id>/out_level-6_<case>.txt, so the merge sees all tactics.
# A file qualifies for a case if it splits into exactly N programs of valid instruction lines.
import os, sys, glob, shutil, hashlib
L6 = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, L6)
from sim6 import parse_input

CASES = ['2-medium', '3-large', '4-enormous']
N = {c: len(parse_input(os.path.join(L6, 'in_level-6_%s.txt' % c))) for c in CASES}

def count_programs(path):
    n, cur = 0, 0
    try:
        for line in open(path, errors='replace'):
            s = line.strip()
            if not s:
                if cur: n += 1; cur = 0
                continue
            t = s.split()
            if len(t) != 6 or any(x not in ('F', 'B', 'T', '-') for x in t): return -1
            cur += 1
    except OSError: return -1
    return n + (1 if cur else 0)

def main():
    root = os.path.join(L6, 'work')
    out = os.path.join(root, 'cand', '_harvest')
    shutil.rmtree(out, ignore_errors=True)
    seen, k = set(), 0
    for path in glob.glob(os.path.join(root, '**', '*.txt'), recursive=True):
        rel = os.path.relpath(path, root)
        if rel.startswith(os.path.join('cand', '_harvest')) or rel.startswith('tmp_best'): continue
        if os.path.getsize(path) < 1000: continue
        n = count_programs(path)
        for c in CASES:
            if n == N[c]:
                h = hashlib.md5(open(path, 'rb').read()).hexdigest()
                if (c, h) in seen: continue
                seen.add((c, h)); k += 1
                d = os.path.join(out, 'h%03d' % k); os.makedirs(d, exist_ok=True)
                shutil.copy(path, os.path.join(d, 'out_level-6_%s.txt' % c))
                print('%s <- %s' % (c, rel))
    print('harvested %d files' % k)

if __name__ == '__main__':
    main()
