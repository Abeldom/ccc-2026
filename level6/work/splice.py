# Build splicing: after every acceptance the relevant board (pattern + skirt) is empty, so the rest of a
# program only depends on the emitter configuration. Every candidate program is cut at its acceptance ticks
# into "builds" (start config -> end config, length). A shortest path over these builds (from the initial
# config, exactly `count` builds) can mix builds from different solvers and repeat good cyclic ones.
# Usage (from level6/): python work/splice.py <case>   -> writes work/cand/splice/out_level-6_<case>.txt
import os, sys, glob, heapq
L6 = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, L6); os.chdir(L6)
from sim6 import parse_input, geometry, simulate
from combine import read_programs

def builds(r, pattern, prog, geo):
    """Split prog into builds: list of (start_config, end_config, lines)."""
    interior, border, per, ray = geo
    n = len(per)
    pat = {(y, x) for y, row in enumerate(pattern) for x, ch in enumerate(row) if ch == 'M'}
    nb = [(0, 2), (0, -2), (1, 1), (1, -1), (-1, 1), (-1, -1)]
    skirt = ({(y + dy, x + dx) for (y, x) in pat for dy, dx in nb} - pat) & interior
    pos, act, matter = [e * r for e in range(6)], [False] * 6, set()
    out, start, cur = [], (tuple(pos), tuple(act)), []
    for line in prog:
        for e, ins in enumerate(line):
            if ins == 'F': pos[e] = (pos[e] + 1) % n
            elif ins == 'B': pos[e] = (pos[e] - 1) % n
            elif ins == 'T': act[e] = not act[e]
        hit = {}
        for e in range(6):
            if act[e]:
                for c in ray[pos[e]][e]: hit[c] = hit.get(c, 0) + 1
        for c, k in hit.items():
            if k == 2: matter.add(c)
            elif k >= 3: matter.discard(c)
        cur.append(line)
        if pat <= matter and not (skirt & matter):
            matter -= pat
            end = (tuple(pos), tuple(act))
            out.append((start, end, cur)); start, cur = end, []
            matter &= ~0 if False else matter      # free cells are irrelevant; nothing else carries over
    return out

def main(case):
    cases = parse_input('in_level-6_%s.txt' % case)
    files = [p for p in glob.glob('work/**/*.txt', recursive=True) + ['out_level-6_%s.txt' % case]
             if 'tmp_best' not in p and os.path.getsize(p) >= 1000]
    progsets = []
    for p in files:
        try: pr = read_programs(p)
        except Exception: continue
        if len(pr) == len(cases): progsets.append(pr)
    geos, result, before, after = {}, [], 0, 0
    best_pub = read_programs('out_level-6_%s.txt' % case)
    for i, (r, g, cnt) in enumerate(cases):
        if r not in geos: geos[r] = geometry(r)
        edges, seen = {}, set()
        for pr in progsets:
            p = pr[i]
            if p is None: continue
            key = tuple(map(tuple, p))
            if key in seen: continue
            seen.add(key)
            for s, e, lines in builds(r, g, p, geos[r]):
                k = (s, e)
                if k not in edges or len(lines) < len(edges[k]): edges[k] = lines
        adj = {}
        for (s, e), lines in edges.items(): adj.setdefault(s, []).append((len(lines), e, lines))
        init = (tuple(e * r for e in range(6)), (False,) * 6)
        # DP over number of builds: dist[config] after k builds
        layer = {init: (0, None)}
        hist = []
        for k in range(cnt):
            nxt = {}
            for c, (d, _) in layer.items():
                for L, e, lines in adj.get(c, []):
                    if e not in nxt or d + L < nxt[e][0]: nxt[e] = (d + L, (c, lines))
            hist.append(nxt); layer = nxt
            if not layer: break
        cur_best = best_pub[i]
        before += len(cur_best)
        if layer:
            endc = min(layer, key=lambda c: layer[c][0])
            parts, c = [], endc
            for k in range(cnt - 1, -1, -1):
                d, (prev, lines) = hist[k][c]; parts.append(lines); c = prev
            prog = [l for part in reversed(parts) for l in part]
            if len(prog) < len(cur_best) and simulate(r, g, prog, geos[r]) >= cnt: cur_best = prog
        after += len(cur_best)
        result.append(cur_best)
    os.makedirs('work/cand/splice', exist_ok=True)
    open('work/cand/splice/out_level-6_%s.txt' % case, 'w', newline='\n').write(
        '\n'.join('\n'.join(' '.join(l) for l in p) + '\n' for p in result))
    print('%s: %d -> %d' % (case, before, after))

if __name__ == '__main__':
    main(sys.argv[1])
