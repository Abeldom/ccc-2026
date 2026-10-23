import re,sys
W='C:/Users/domok/OneDrive/Desktop/ccc 2026/level6/work/'
best={}
for l in open(W+'best_lengths_2-medium.txt'):
    if l.startswith('#'): continue
    i,r,P,c,b=map(int,l.split()); best[i]=(r,P,c,b)
def exh(fn):
    d={}
    for l in open(W+'cycles/'+fn):
        m=re.match(r'(\d+) .*found=(\d) bestP=(-?\d+)',l)
        if m: d[int(m.group(1))]=(int(m.group(2)),int(m.group(3)))
        elif 'skip' in l: d[int(l.split()[0])]=('skip',)
    return d
e2=exh('exh2_medium.txt'); e3=exh('exh3_medium.txt')
tot_gain=0
rows=[]
for i in sorted(best):
    r,P,c,b=best[i]
    s2=e2.get(i); s3=e3.get(i)
    if s2 and s2[0]==1: lb=2*c; kind='L2'
    elif s2 and s2[0]==0:
        if s3 and s3[0]==1: lb=2+3*(c-1); kind='L3'
        elif s3 and s3[0]==0: lb=2+4*(c-1); kind='noL3'  # rough
        else: lb=2+3*(c-1); kind='?3'
    else: lb=0; kind='skip'
    g=b-lb
    rows.append((g,i,r,P,c,b,kind,s2,s3))
rows.sort(reverse=True)
for x in rows: print(*x)
from collections import Counter
print(Counter(x[6] for x in rows))
print('sum gap',sum(x[0] for x in rows if x[6]!='skip'))
