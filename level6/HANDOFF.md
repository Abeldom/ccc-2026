# Level 6 handoff

## Task
For each pattern, output a program (lines of 6 instructions `F B T -`) that makes the nanofabricator accept the
pattern `count` times. Score = total number of lines; lower is better. Only medium, large, enormous count.

## Current best (all valid, in this folder)
| file | lines |
|---|---|
| out_level-6_2-medium.txt | 5941 |
| out_level-6_3-large.txt | 5334 |
| out_level-6_4-enormous.txt | 5054 |
Baseline (simple walker) was 38575 / 32946 / 30445.

## Tools (in level6/)
- `sim6.py` — exact validator, verified tick-by-tick against `sandbox-level-6.html`.
  `python sim6.py in_level-6_3-large.txt out_level-6_3-large.txt` prints invalid count + score.
- `baseline.py` — always-valid walker solver.
- `combine.py` — per pattern, keeps the shortest valid program across candidate folders.
- `work/update_best.py` — merges all `work/cand/*` folders, runs the post-optimizer, publishes only improvements.
- `work/harvest.py` — collects complete output files from solver folders into `work/cand/_harvest/`.
- `work/postopt/postopt.cpp` — shortens valid programs (delete or merge lines, re-validated).
- `work/watch_gain.py` — loop: harvest + merge every 3 minutes, logs to `UPDATES.txt`.
- `work/BRIEF.md` — full rules, geometry facts and insights (read this first).
- `work/best_lengths_<case>.txt` — per-pattern best line counts.

## Key facts (verified)
- Cell coords: u=(x-y-r)/2, v=(x+y-3r)/2, h=y-r=v-u; interior iff |u|,|v|,|h| <= r-1.
- Emitter e at signed offset k covers one full line: axis [H,U,V,H,U,V][e], line sign[e]*k with
  sign=[-1,+1,+1,+1,-1,-1]; |k| >= r gives no beam. All start on the central line, inactive.
- 2 hits = matter, 3+ hits = destroyed. Two facing emitters (0&3, 1&4, 2&5) on the same line create a whole
  line in one tick; a crossing third beam removes a cell ("cut").
- Only pattern cells P and their empty neighbours K matter; other cells are irrelevant.
- After each acceptance P and K are empty, so only the emitter configuration carries over. Repeats should be a
  short cycle; a fixed setup (repeat every tick) is possible for only 9 medium patterns, 1 large, 0 enormous.
- Proven lower bound for medium is about 4019 (non-static patterns need at least 2 ticks per repeat).

## Solver approaches tried (code in work/<name>/)
- shots — macro planner over beam setups (produced most of the best).
- beam — tick-level beam search with exact simulation, plus re-solving the tails of the best programs.
- cycles — short periodic cycles for repeated patterns (medium).
- sweep — full-line sweeps with cutting beams (large/enormous).
- lines — user idea: cover the pattern with whole lines (opposite emitter pairs) plus cuts, ordered to save travel.
- sparse — exact or wide search for small patterns.

## Where the room is
- Medium: about 2.9 lines per repeat; target is 2 (cycles of length 2). The big medium patterns with high counts
  cost the most.
- Large/enormous: big patterns (130-230 cells) cost 110-190 lines. They lie on about 15 lines of one direction,
  so line-based building (whole lines plus cuts, about 4-5 skirt cells per line) should cut them a lot.
