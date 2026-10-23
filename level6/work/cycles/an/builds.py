# print accept ticks (build lengths) of published programs
import sys
L6='C:/Users/domok/OneDrive/Desktop/ccc 2026/level6'
sys.path.insert(0,L6)
from sim6 import parse_input, parse_output, geometry
def acc_ticks(r, pattern, prog, geo):
    interior, border, per, ray = geo
    n = len(per)
    pat = {(y, x) for y, row in enumerate(pattern) for x, ch in enumerate(row) if ch == 'M'}
    nb = [(0, 2), (0, -2), (1, 1), (1, -1), (-1, 1), (-1, -1)]
    skirt = ({(y + dy, x + dx) for (y, x) in pat for dy, dx in nb} - pat) & interior
    pos = [e * r for e in range(6)]; act = [False] * 6; matter = set(); res=[]
    for t,line in enumerate(prog):
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
        if pat <= matter and not (skirt & matter):
            matter -= pat; res.append(t+1)
    return res
case=sys.argv[1]; f=sys.argv[2] if len(sys.argv)>2 else L6+'/out_level-6_%s.txt'%case
cases=parse_input(L6+'/in_level-6_%s.txt'%case); progs=parse_output(f,len(cases))
geos={}
for i,((r,g,cnt),p) in enumerate(zip(cases,progs)):
    if cnt<2: continue
    if r not in geos: geos[r]=geometry(r)
    a=acc_ticks(r,g,p,geos[r])[:cnt]
    d=[a[0]]+[a[j]-a[j-1] for j in range(1,len(a))]
    print(i,'r',r,'cnt',cnt,'len',len(p),'builds',d)
