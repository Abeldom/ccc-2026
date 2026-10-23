import sys
sys.path.insert(0, '../..')
from sim6 import parse_input
def coords(y, x, r):
    return (x - y - r) // 2, (x + y - 3 * r) // 2, y - r
def read_progs(path):
    progs, cur = [], []
    for line in open(path).read().split('\n'):
        if line.strip()=='':
            if cur: progs.append(cur); cur=[]
            continue
        cur.append(line.split())
    if cur: progs.append(cur)
    return progs
case = sys.argv[1]
cases = parse_input('../../in_level-6_%s.txt' % case)
bl = read_progs('../cand/baseline/out_level-6_%s.txt' % case)
tot=0
for i,(r,g,cnt) in enumerate(cases):
    P=set()
    for y,row in enumerate(g):
        for x,ch in enumerate(row):
            if ch=='M': u,v,h=coords(y,x,r); P.add((u,v))
    K=set()
    for (u,v) in P:
        for du,dv in [(1,1),(-1,-1),(0,1),(0,-1),(1,0),(-1,0)]:
            c=(u+du,v+dv)
            if c in P: continue
            if max(abs(c[0]),abs(c[1]),abs(c[1]-c[0]))<=r-1: K.add(c)
    tot+=len(bl[i])
    print(i, 'r',r,'cnt',cnt,'P',len(P),'K',len(K),'base',len(bl[i]), 'ratio %.2f'%(len(bl[i])/len(P)/cnt))
print('total', tot)
