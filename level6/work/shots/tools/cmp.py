import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from perpat import lens, best, L6
case = sys.argv[1]
files = sys.argv[2:]
b = best(case)
L = [lens(f) for f in files]
for i, row in enumerate(b):
    print(i, row[1], row[2], row[3], 'best', row[4], ' '.join('%5d' % l[i] for l in L))
