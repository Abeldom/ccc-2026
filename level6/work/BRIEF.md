# Level 6 brief (Cloudflight Coding Contest "Nanofabricator") — shared by all solver agents

Folder: `C:\Users\domok\OneDrive\Desktop\ccc 2026\level6` (call it L6). Python 3.14 (no numpy) and
g++ 16 (MSYS2, `g++ -O2 -std=c++20`, OpenMP available via `-fopenmp`) are installed; 16 CPU cores.

## Task
Each input file `L6/in_level-6_<case>.txt` (cases: `0-example 1-small 2-medium 3-large 4-enormous`) holds N
patterns. For each pattern: header `r rows cols count`, then the board drawing (`-`, `*` border, `O` empty,
`M` pattern matter). Output: for each pattern one program (lines of 6 space-separated instructions from
`F B T -`), programs separated by ONE empty line, in input order. A program is valid iff simulating it
accepts the pattern at least `count` times. **Score = total number of instruction lines over all programs;
lower is better.** Every line counts (including all-`-` lines).

Sizes: small 10 patterns (r=3); medium 210 (r=3..5, |P| 4..59, count 1..19, sum count 2045);
large 190 (r=6..8, |P| 6..131, count 1..6, sum 689); enormous 140 (r=9..12, |P| 18..297, count 1..2).

## Exact rules (verified against the official sandbox JS; `L6/sim6.py` is the reference validator)
- Board radius r: char grid 2r+1 rows x 4r+1 cols. Row i (d=|i-r|) has cells at cols d, d+2, ..., 4r-d.
  Border = first/last cell of each row + whole top & bottom row. Others are interior (3r(r-1)+1 cells).
- 6 emitters ride the 6r border cells (clockwise list starting at the left corner going up). Emitter e
  (0-based) starts at index e*r, inactive. Instructions per tick, one per emitter, all simultaneous:
  `F` = +1 clockwise, `B` = -1, `T` = toggle active, `-` = nothing. (An emitter cannot move and toggle in
  the same tick.)
- Each tick: (1) apply the instruction line; (2) every ACTIVE emitter fires a straight ray in its fixed
  direction (e0 E, e1 SE, e2 SW, e3 W, e4 NW, e5 NE) through interior cells; count hits per cell:
  exactly 2 hits -> cell becomes matter; >=3 hits -> cell becomes empty (destroyed); 0/1 hit -> unchanged.
  Matter persists between ticks. (3) Detection: if every pattern `M` cell is matter AND no "skirt" cell
  is matter, the pattern is accepted (+1) and exactly the pattern cells' matter is removed. At most one
  acceptance per tick. Skirt K = interior cells not in P that are hex-adjacent to a P cell.
  All other interior cells ("free" cells) are never checked: matter there is harmless forever.
- `python L6/sim6.py <input> <output>` validates a whole file and prints the score. `L6/combine.py`
  picks the shortest valid program per pattern across candidate dirs.

## Geometry facts (verified for r=2..12) — use these, they make everything simple
Cell at char (y,x): u = (x-y-r)/2 (SE-diagonal index), v = (x+y-3r)/2 (SW-diagonal index), h = y-r = v-u
(row index). Interior <=> |u|,|v|,|h| <= r-1. Hex neighbours in (u,v): (+1,+1) E, (-1,-1) W, (0,+1) SE,
(0,-1) NW, (+1,0) NE, (-1,0) SW. Three line families ("axes"): H lines {h=c}, U lines {u=c}, V lines {v=c},
each with 2r-1 interior lines c in [-(r-1), r-1]. Each cell lies on exactly one line of each axis; two
lines of different axes cross in exactly one point (maybe outside the board).
Emitter e at signed offset k from its start (index = e*r + k mod 6r) covers exactly ONE full line:
axis = [H,U,V,H,U,V][e], line = sign[e]*k with sign = [-1,+1,+1,+1,-1,-1], valid when |k| <= r-1.
For |k| >= r the ray contains no interior cell ("dead zone" = effectively off even if active).
So: every emitter starts on the central line of its axis (offset 0, the longest chord); reaching line c
costs |c| moves; the two emitters of an axis (0&3 on H, 1&4 on U, 2&5 on V) are interchangeable, and when
both sit on the same line ("duplicate") that WHOLE line gets 2 hits -> a full line of matter in one tick
(minus cells also hit by another active line, which get 3+ hits and are destroyed). To move a duplicate
pair together emitter 0 and 3 use opposite letters (e.g. e0 `B` + e3 `F` = both move to h+1); same for
1/4 and 2/5.
Per tick with lines A (H), B (U), C (V) (multisets of size <= 2): created set S = cells with exactly 2
hits = pairwise crossings not on a third line, plus duplicate lines minus crossings; destroyed D = 3+ hits.

## Key structural insights
1. Only P and K cells matter. Free cells can absorb any junk.
2. Canonical state: right after an acceptance (and at the start), all P cells are empty and all K cells
   are empty. So the only state carried from one acceptance to the next is the emitter configuration
   (offset+active of the 6 emitters). The total cost = sum of "builds"; build(C) = shortest tick sequence
   from config C with an empty board to the next detection; it ends in some config C'. For count>1 the
   steady state is a (min-mean) cycle of builds; e.g. a periodic schedule with period L costs ~ L per
   acceptance after the first.
3. Static configs: if a single config C has S(C) ⊇ P and S(C) ∩ K = ∅, then move emitters (inactive) to
   their lines, toggle them all in one tick, and every following all-`-` tick re-accepts:
   cost = max|offset| + count. Only ~9/210 medium patterns allow this (probe: `L6/work/static_probe.cpp`).
4. The two emitters of an axis can be split (two distinct lines) or doubled (full line). Moving an active
   emitter off the board edge (|k| = r) switches it off for the same 1 instruction a toggle costs.
5. Walker: two active emitters of different axes create exactly their crossing cell; moving each by
   <=1 per tick moves the crossing to any neighbour (or a distance-2 jump). `L6/baseline.py` paints with
   one walker along a DFS tour: always valid, ~1-2 ticks per cell per acceptance.
6. Lower bound per cell: a cell at (u,v,h) cannot exist before tick 1 + min(|u|,|v|,|h|) (needs a
   duplicate on its nearest line) — e.g. tick 1 with e0+e3 toggled creates the whole middle row.

## Current reference scores (lower is better)
baseline (walker): example 100, small 473, medium 38575, large 32946, enormous 30445.
Official example output: 446 for the example file.

## Deliverables for every solver agent
- Code in `L6/work/<your-name>/` (any language). Candidate outputs in `L6/work/cand/<your-name>/` named
  `out_level-6_<case>.txt` for all 5 cases, each with exactly N programs (where your method fails or is
  worse, copy the baseline program for that pattern from `L6/work/cand/baseline/`).
- Validate every file with `python L6/sim6.py`; never hand in an invalid program.
- Do not modify sim6.py, baseline.py, combine.py, the inputs, or other agents' folders.
- Keep any single long computation in the background or under ~30 minutes; use multiple cores.
