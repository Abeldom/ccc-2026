# per-pattern lower bounds for medium (from exhaustive 2-cycle search + static hints) vs my current candidate
import re,sys
sys.path.insert(0,'C:/Users/domok/OneDrive/Desktop/ccc 2026/level6')
W='C:/Users/domok/OneDrive/Desktop/ccc 2026/level6/work/'
from sim6 import parse_output, parse_input
cases=parse_input(W+'../in_level-6_2-medium.txt')
cur=parse_output(W+'cand/cycles/out_level-6_2-medium.txt',len(cases))
hint={}
for l in open(W+'cycles/two_medium.txt'):
    m=re.match(r'(\d+) .*static=(-?\d+) two=(-?\d+) bestP=(-?\d+)',l)
    if m: hint[int(m.group(1))]=tuple(map(int,m.groups()[1:]))
res=[]
for i,(r,g,c) in enumerate(cases):
    st,two,bp=hint.get(i,(-1,-1,-1))
    if st>=0: lb=c+st if False else c+1
    elif two>0: lb=2*c if bp==0 else 2*c+1
    elif two==0: lb=3*c
    else: lb=2*c
    res.append((len(cur[i])-lb,i,c,len(cur[i]),lb))
res.sort(reverse=True)
imp=[x for x in res if x[0]>0 and x[2]>1]
print('improvable count>1:',len(imp),'gap sum',sum(x[0] for x in imp))
print(' '.join(str(x[1]) for x in imp))
open(W+'cycles/imp.idx','w').write(' '.join(str(x[1]) for x in imp))
for x in imp[:60]: print(x)
