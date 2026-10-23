# Keep level6/out_level-6_<case>.txt = best valid solution found so far.
# For each case (example, small, medium, large, enormous): take the shortest valid program per pattern
# across all candidate dirs in work/cand/ and the current best file, run the post-optimizer on it,
# re-validate with sim6.py and replace the published file only if the score got lower.
# Also writes level6/SCORES.txt. Usage (from level6/): python work/update_best.py
import os, subprocess, sys, glob

L6 = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CASES = ['0-example', '1-small', '2-medium', '3-large', '4-enormous']
sys.path.insert(0, L6)
from sim6 import parse_input, parse_output, geometry, simulate

def score_of(case, path):
    """Total lines if every program is valid, else None."""
    if not os.path.exists(path): return None
    cases = parse_input(os.path.join(L6, 'in_level-6_%s.txt' % case))
    try: progs = parse_output(path, len(cases))
    except AssertionError: return None
    if len(progs) != len(cases): return None
    geos = {}
    for (r, g, cnt), p in zip(cases, progs):
        if r not in geos: geos[r] = geometry(r)
        if simulate(r, g, p, geos[r]) < cnt: return None
    return sum(len(p) for p in progs)

def main():
    os.chdir(L6)
    import time
    lock = os.path.join('work', '.update.lock')
    for _ in range(600):                                       # one merge at a time (background loop + manual runs)
        try: fd = os.open(lock, os.O_CREAT | os.O_EXCL); os.close(fd); break
        except FileExistsError:
            if time.time() - os.path.getmtime(lock) > 1800: os.remove(lock)   # stale lock
            time.sleep(1)
    try: run()
    finally: os.remove(lock)

def run():
    dirs = sorted(d for d in glob.glob(os.path.join('work', 'cand', '*')) + glob.glob(os.path.join('work', 'cand', '_harvest', '*'))
                  if os.path.isdir(d))
    lines = []
    for case in CASES:
        pub = 'out_level-6_%s.txt' % case
        tmpdir = os.path.join('work', 'tmp_best'); os.makedirs(tmpdir, exist_ok=True)
        merged = os.path.join(tmpdir, 'merged_%s.txt' % case)
        opt = os.path.join(tmpdir, 'out_level-6_%s.txt' % case)
        srcs = [d for d in dirs if os.path.exists(os.path.join(d, pub))]
        cur = score_of(case, pub)
        if cur is not None: srcs.append('.')                     # current published best competes too
        subprocess.run([sys.executable, 'combine.py', case, merged] + srcs, check=True, capture_output=True)
        subprocess.run([os.path.join('work', 'postopt', 'postopt.exe'), 'in_level-6_%s.txt' % case, merged, opt],
                       check=True, capture_output=True)
        new = score_of(case, opt)
        if new is not None and (cur is None or new < cur):
            open(pub, 'w', newline='\n').write(open(opt).read())
            cur = new
        lines.append('%-12s %8s lines' % (case, cur))
        print(lines[-1], flush=True)
    open('SCORES.txt', 'w').write('Current best (total instruction lines, lower is better):\n' + '\n'.join(lines) + '\n')

if __name__ == '__main__':
    main()
