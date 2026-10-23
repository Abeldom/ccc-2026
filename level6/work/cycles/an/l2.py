exec(open('analyze.py').read().split('tot_gain=0')[0])
tot=0; idx=[]
for i in sorted(best):
    r,P,c,b=best[i]
    s2=e2.get(i)
    if s2 and s2[0]==1 and c>1:
        simple=s2[1]+2*c
        if b>2*c: idx.append(i)
        if b>simple: print(i,'cnt',c,'best',b,'simple',simple); tot+=b-simple
print('total simple gain',tot)
print(' '.join(map(str,idx)))
open('../l2.idx','w').write(' '.join(map(str,idx)))
