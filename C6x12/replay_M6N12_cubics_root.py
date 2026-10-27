#!/usr/bin/env python3
"""Independent positive coverage and generic group-law checks for all cubics."""
from pathlib import Path
import ast,json
from collections import Counter
import sys;sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'common'))
from M3N18_finite_arithmetic_root import FF
P=Path(__file__).resolve().parent
a=ast.literal_eval((P/'data'/'M6N12_cubic_control_arithmetic_root.txt').read_text().replace('[*','[').replace('*]',']'))
fps,f,keys,recs=a;assert len(recs)==48 and len(set(map(tuple,keys)))==96
for p,fp in fps:
    tr=[];cl=[]
    for n in range(1,5):
        tr.append(-n*fp[18-n]-sum(fp[18-j]*tr[n-j-1] for j in range(1,n)))
        v=p**n+1-tr[-1]-sum(d*cl[d-1] for d in range(1,n) if n%d==0)
        assert v%n==0;cl.append(v//n)
    assert cl=={7:[16,8,96,468],13:[32,144,576,7116]}[p]
F=FF(7,f)
class EC6:
    def __init__(self,aa):self.a1,self.a2,self.a3,self.a4,self.a6=aa
    def has(self,P):
        x,y=P
        return F.add(F.add(F.mul(y,y),F.mul(self.a1,F.mul(x,y))),F.mul(self.a3,y))==F.add(F.add(F.add(F.pow(x,3),F.mul(self.a2,F.mul(x,x))),F.mul(self.a4,x)),self.a6)
    def slope(self,P):
        x,y=P
        return F.div(F.sub(F.add(F.add(F.sc(3,F.mul(x,x)),F.sc(2,F.mul(self.a2,x))),self.a4),F.mul(self.a1,y)),F.add(F.sc(2,y),F.add(F.mul(self.a1,x),self.a3)))
    def add(self,P,Q):
        if P is None:return Q
        if Q is None:return P
        x,y=P;X,Y=Q
        if x==X:
            if F.add(F.add(y,Y),F.add(F.mul(self.a1,x),self.a3))==F.z:return None
            sl=self.slope(P)
        else:sl=F.div(F.sub(Y,y),F.sub(X,x))
        nu=F.sub(y,F.mul(sl,x))
        xx=F.sub(F.sub(F.add(F.mul(sl,sl),F.mul(self.a1,sl)),self.a2),F.add(x,X))
        yy=F.sub(F.neg(F.mul(F.add(sl,self.a1),xx)),F.add(nu,self.a3))
        return xx,yy
    def mul(self,n,P):
        R=None
        while n:
            if n&1:R=self.add(R,P)
            P=self.add(P,P);n//=2
        return R
    def uv(self,P,Q):
        x,y=Q;a3=F.add(F.sc(2,y),F.add(self.a3,F.mul(self.a1,x)))
        a2=F.add(F.sc(3,x),self.a2)
        a4=F.sub(F.add(F.add(F.sc(3,F.mul(x,x)),F.sc(2,F.mul(self.a2,x))),self.a4),F.mul(self.a1,y))
        sl=F.div(a4,a3);aa1=F.add(F.sc(2,sl),self.a1)
        aa2=F.sub(F.sub(a2,F.mul(sl,sl)),F.mul(self.a1,sl));scale=F.div(aa2,a3)
        b=F.neg(F.div(F.pow(aa2,3),F.pow(a3,2)));c=F.sub(F.o,F.mul(aa1,scale))
        r=F.div(b,c);s=F.div(F.mul(c,c),F.sub(b,c));u=F.div(F.sub(r,F.o),F.sub(r,s))
        q=F.mul(F.mul(scale,scale),F.sub(P[0],x))
        v=F.sub(F.sub(F.div(F.mul(F.mul(q,F.mul(u,u)),F.pow(F.sub(u,F.o),2)),F.add(u,F.o)),F.mul(u,u)),F.o)
        return u,v
covered=set();wpvals=set()
for i,record in recs:
    vals=list(map(tuple,record));assert len(vals)==11
    E=EC6(vals[:5]);PP,QQ=tuple(vals[5:7]),tuple(vals[7:9]);uv=tuple(vals[9:11])
    assert E.has(PP) and E.has(QQ)
    assert E.mul(6,PP) is None and E.mul(3,PP) is not None and E.mul(2,PP) is not None
    assert E.mul(12,QQ) is None and E.mul(6,QQ) is not None and E.mul(4,QQ) is not None
    assert E.mul(3,PP)!=E.mul(6,QQ)
    R,S=E.mul(2,PP),E.mul(4,QQ);assert R not in [S,E.mul(2,S)]
    x,y=R;X,Y=S
    lp=F.sub(F.sub(Y,y),F.mul(E.slope(R),F.sub(X,x)))
    lq=F.sub(F.sub(y,Y),F.mul(E.slope(S),F.sub(x,X)))
    wp=F.neg(F.div(lp,lq));assert wp!=F.o and F.pow(wp,3)==F.o;wpvals.add(wp)
    assert E.uv(PP,QQ)==uv
    ims={F.key(E.uv(E.mul(t,PP),E.mul(t,QQ))) for t in [1,5]}
    assert len(ims)==2 and not covered.intersection(ims);covered.update(ims)
assert covered==set(map(tuple,keys)) and len(wpvals)==1
log=(P/'logs'/'control_6_12_D4_cubics_root.log').read_text()
assert 'M6N12_CUBIC_GENERIC_REFERENCE 42' in log and 'M6N12_CUBIC_CONTROL_COMPLETE' in log
assert log.count('M6N12_CUBIC_EXCLUDED')==48 and not any(x in log for x in ['Runtime error','Assertion failed'])
out={'degree':8,'group':[6,12],'cubic_places':96,'diamond_representatives':48,
 'independent_level_orders':[6,12],'subgroup_order':72,'all_positive_coordinates_and_images_checked':True,
 'nonprincipal_A13_records':48,'generic_reference_degree':42,'Newton_counts_checked':True}
(P/'data'/'replay_M6N12_cubics_root.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2));print('M6N12_CUBIC_INDEPENDENT_REPLAY_DONE')
