# Merge loop that stops (so the main session gets notified) as soon as the counted total
# (medium + large + enormous) has dropped by >= 50 lines since the last report.
# Every 3 minutes: harvest all solver outputs, merge + post-optimize + publish best, log to level6/UPDATES.txt.
# Usage (from level6/): python work/watch_gain.py <last reported total>
import os, sys, time, subprocess, datetime
L6 = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(L6)
COUNTED = ['2-medium', '3-large', '4-enormous']

def scores():
    out = {}
    for line in open('SCORES.txt'):
        p = line.split()
        if len(p) >= 2 and p[0] in COUNTED: out[p[0]] = int(p[1])
    return out

base = int(sys.argv[1])
step = int(sys.argv[2]) if len(sys.argv) > 2 else 50
while True:
    subprocess.run([sys.executable, 'work/harvest.py'], capture_output=True)
    subprocess.run([sys.executable, 'work/update_best.py'], capture_output=True)
    s = scores(); total = sum(s.values())
    stamp = datetime.datetime.now().strftime('%H:%M')
    msg = '%s  medium %d  large %d  enormous %d  total %d' % (stamp, s['2-medium'], s['3-large'], s['4-enormous'], total)
    open('UPDATES.txt', 'a').write(msg + '\n')
    if total <= base - step:
        print('GAIN %d lines: %s' % (base - total, msg)); break
    time.sleep(180)
