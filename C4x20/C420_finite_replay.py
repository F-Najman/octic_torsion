#!/usr/bin/env python3
"""Independent finite arithmetic, orbit, Newton, stabilizer and coverage replay.
Only Python integer arithmetic is used. No Magma or imported proof helper.
--report writes a new summary; the default is strictly read-only.
"""
import argparse,ast,collections,hashlib,itertools,json
from pathlib import Path
P=Path(__file__).resolve().parents[1]/'C4x20'/'data';NC=P
def read(p):return ast.literal_eval(p.read_text().replace('[*','[').replace('*]',']'))
def freeze(a):return tuple(map(freeze,a)) if isinstance(a,list) else a

class Field:
    def __init__(self,poly):
        assert poly[-1]==1;self.m=len(poly)-1;self.poly=tuple(poly)
        self.zero=(0,)*self.m;self.one=(1,)+(0,)*(self.m-1)
    def c(self,n):return (n%3,)+(0,)*(self.m-1)
    def add(self,a,b):return tuple((x+y)%3 for x,y in zip(a,b))
    def neg(self,a):return tuple(-x%3 for x in a)
    def sub(self,a,b):return self.add(a,self.neg(b))
    def mul(self,a,b):
        v=[0]*(2*self.m-1)
        for i,x in enumerate(a):
            if x:
                for j,y in enumerate(b):v[i+j]+=x*y
        for k in range(2*self.m-2,self.m-1,-1):
            z=v[k]%3
            if z:
                for j in range(self.m):v[k-self.m+j]-=z*self.poly[j]
        return tuple(x%3 for x in v[:self.m])
    def pow(self,a,n):
        ans=self.one
        while n:
            if n&1:ans=self.mul(ans,a)
            a=self.mul(a,a);n//=2
        return ans
    def inv(self,a):
        assert a!=self.zero
        z=self.pow(a,3**self.m-2);assert self.mul(a,z)==self.one;return z
    def div(self,a,b):return self.mul(a,self.inv(b))
    def key(self,q,t):
        ans=[]
        for _ in range(self.m//2):
            ans.append(q+t);q=self.pow(q,9);t=self.pow(t,9)
        return min(ans)

class Curve:
    def __init__(self,F,B):self.f=F;self.B=B
    def valid(self,P):
        if P is None:return True
        f=self.f;x,y=P
        return f.add(f.add(f.mul(y,y),f.mul(x,y)),f.mul(self.B,y))==f.add(f.pow(x,3),f.mul(self.B,f.mul(x,x)))
    def neg(self,P):
        if P is None:return None
        f=self.f;x,y=P;return x,f.neg(f.add(f.add(y,x),self.B))
    def add(self,P,Q):
        if P is None:return Q
        if Q is None:return P
        f=self.f;x,y=P;u,v=Q
        if x==u:
            if Q==self.neg(P):return None
            assert P==Q
            lam=f.div(f.sub(f.mul(f.c(2),f.mul(self.B,x)),y),f.add(f.add(f.mul(f.c(2),y),x),self.B))
        else:lam=f.div(f.sub(v,y),f.sub(u,x))
        nu=f.sub(y,f.mul(lam,x))
        xx=f.sub(f.sub(f.sub(f.add(f.mul(lam,lam),lam),self.B),x),u)
        yy=f.sub(f.sub(f.neg(f.mul(f.add(lam,f.one),xx)),nu),self.B)
        R=(xx,yy);assert self.valid(R);return R
    def mul(self,n,P):
        if n<0:return self.mul(-n,self.neg(P))
        R=None
        while n:
            if n&1:R=self.add(R,P)
            P=self.add(P,P);n//=2
        return R
    def normalized(self,P,Q):
        f=self.f;r,t0=P
        s=f.div(f.sub(f.mul(f.c(2),f.mul(self.B,r)),t0),f.add(f.add(f.mul(f.c(2),t0),r),self.B))
        a1=f.add(f.one,f.mul(f.c(2),s))
        a2=f.sub(f.sub(self.B,s),f.mul(s,s))
        scale=f.inv(a1);b=f.mul(f.mul(scale,scale),a2)
        xx=f.mul(f.pow(scale,2),f.sub(Q[0],r))
        yy=f.mul(f.pow(scale,3),f.sub(f.sub(Q[1],t0),f.mul(s,f.sub(Q[0],r))))
        q=f.add(f.one,f.div(f.mul(f.c(2),b),xx))
        t=f.add(f.one,f.div(f.mul(f.c(4),yy),f.mul(f.add(q,f.one),xx)))
        assert b==f.div(f.mul(f.sub(f.mul(q,q),f.one),f.sub(f.mul(t,t),f.one)),f.c(16))
        return q,t

def params():
    return [(a,b,c,d,u) for u in [1,2] for a,b,c,d in itertools.product(range(4),repeat=4) if (a*d-b*c)%4==1]
def groupmul(g,h):
    a,b,c,d,u=g;e,f,k,l,v=h
    M=((a*e+b*k)%4,(a*f+b*l)%4,(c*e+d*k)%4,(c*f+d*l)%4)
    uv=u*v%5
    if uv in [3,4]:M=tuple(-x%4 for x in M);uv=-uv%5
    return M+(uv,)

def cusps():
    # Tate labels: Q20=zeta20^a q^(b/20), P4=zeta4^c q^(d/4).
    from math import gcd
    labels={(a,b,c,d) for a,b in itertools.product(range(20),repeat=2)
            if gcd(gcd(a,b),20)==1 for c,d in itertools.product(range(4),repeat=2) if (c*b-d*a)%4==1}
    def orbit(x):
        a,b,c,d=x
        return frozenset((s*(a+k*b)%20,s*b%20,s*(c+k*d)%4,s*d%4) for s in [-1,1] for k in range(20))
    cs=sorted(set(map(orbit,labels)),key=lambda s:sorted(s));ix={x:i for i,c in enumerate(cs) for x in c}
    def perm(t):
        ans=[]
        for c in cs:
            imgs={ix[(t*a%20,b,t*e%4,d)] for a,b,e,d in c};assert len(imgs)==1;ans.append(imgs.pop())
        assert len(set(ans))==48;return ans
    gal={t:perm(t) for t in [1,9,13,17]};left=set(range(48));cycles=[]
    while left:
        x=min(left);cy=[]
        while x not in cy:cy.append(x);left.remove(x);x=gal[9][x]
        assert x==cy[0];cycles.append(cy)
    sig=collections.Counter(map(len,cycles));assert sig=={1:24,2:12}
    assert all(all(pi[cy[0]]==cy[0] for pi in gal.values()) for cy in cycles if len(cy)==1)
    assert all(gal[9][cy[0]] in cy for cy in cycles if len(cy)<=2)
    return {'geometric':48,'closed_degrees':dict(sig),'rational_cusps_global_over_Qi':24,
            'admissible_q_for_complement_1':[13,17,29],'admissible_q_for_complement_2':[29]}

def main():
    a=read(P/'C420_positive_fresh_arithmetic.txt');old=read(NC/'preflight_4_20_D4_Aut96_arithmetic.txt')
    pars=params();assert len(pars)==96
    for g in pars:
        assert groupmul(g,(1,0,0,1,1))==g
        assert {groupmul(g,h) for h in pars}==set(pars)
    counts=[];orbit_sizes=[]
    for d,poly,keys,labels,orbs,indices in a:
        assert [d,keys,labels,orbs,indices]==old[d-2]
        F=Field(poly);seen=set();counts.append(len(keys));orbit_sizes.append([])
        assert len(set(map(tuple,keys)))==len(keys)
        for key in keys:
            q,t=tuple(key[:2*d]),tuple(key[2*d:]);assert F.key(q,t)==tuple(key)
        for i,(images,index) in enumerate(zip(orbs,indices),1):
            key=keys[index-1];q,t=tuple(key[:2*d]),tuple(key[2*d:])
            B=F.div(F.mul(F.sub(F.mul(q,q),F.one),F.sub(F.mul(t,t),F.one)),F.c(16));E=Curve(F,B)
            PP=(F.zero,F.zero)
            xx=F.div(F.mul(F.add(q,F.one),F.sub(F.mul(t,t),F.one)),F.c(8))
            yy=F.div(F.mul(F.mul(F.pow(F.add(q,F.one),2),F.pow(F.sub(t,F.one),2)),F.add(t,F.one)),F.c(32));QQ=(xx,yy)
            assert E.valid(PP) and E.valid(QQ)
            assert E.mul(4,PP) is None and E.mul(2,PP) is not None
            assert E.mul(20,QQ) is None and E.mul(10,QQ) is not None and E.mul(4,QQ) is not None
            assert E.mul(2,PP)!=E.mul(10,QQ)
            Q4=E.mul(5,QQ);Q5=E.mul(4,QQ)
            actual=[]
            for aa,bb,cc,dd,u in pars:
                R=E.add(E.mul(aa,PP),E.mul(bb,Q4));S=E.add(E.add(E.mul(cc,PP),E.mul(dd,Q4)),E.mul(-u,Q5))
                x,y=E.normalized(R,S);actual.append(F.key(x,y))
            assert actual==list(map(tuple,images)),(d,i,'independent group-action mismatch')
            orb=set(actual);expected={tuple(k) for k,l in zip(keys,labels) if l==i}
            assert orb==expected and not(seen&orb);seen|=orb;orbit_sizes[-1].append(len(orb))
        assert seen==set(map(tuple,keys))
        print('C420_FINITE_FIELD_ORBITS_PASS',d,len(keys),len(orbs),orbit_sizes[-1],flush=True)
    assert counts==[96,320,888]
    g,rows=read(NC/'preflight_4_20_D4_Frob_arithmetic.txt');assert g==25
    for p,step,poly,expected,order in rows:
        assert len(poly)==51 and poly[-1]==1
        # Full coefficient functional equation, Newton identities, closed places.
        assert all(poly[k]==p**(25-k)*poly[50-k] for k in range(26))
        tr=[0]
        for n in range(1,4*step+1):tr.append(-n*poly[50-n]-sum(poly[50-j]*tr[n-j] for j in range(1,n)))
        cl=[]
        for n in range(1,5):
            num=p**(n*step)+1-tr[n*step]-sum(e*cl[e-1] for e in range(1,n) if n%e==0)
            assert num%n==0;cl.append(num//n)
        assert cl==expected
        assert order==sum(poly)*(sum(x*(-1)**k for k,x in enumerate(poly)) if step==2 else 1)
        if p==3:assert cl==[24,108,320,888] and [cl[j]-[24,12,0,0][j] for j in range(4)]==[0,96,320,888]
    print('C420_NEWTON_COUNTS_PASS',flush=True)
    cusp=cusps();print('C420_CUSP_GALOIS_PASS',cusp,flush=True)
    tests,rem=read(NC/'sieve_4_20_D4_Aut96_arithmetic.txt');tail,trem=read(NC/'M4N20_pair_tail_root_arithmetic.txt')
    assert len(tests)==103 and len(tail)==111
    complete=collections.defaultdict(dict)
    for typ,d,i,j,q,pr in tests:
        key=(typ,d,i,j);assert q not in complete[key];complete[key][q]=pr
    for i,j,q,pr in tail:
        key=(2,2,i,j);assert q not in complete[key];complete[key][q]=pr
    expected={(1,d,i,0) for d,num in [(2,2),(3,4),(4,12)] for i in range(1,num+1)}|{(2,2,i,j) for i in [1,2] for j in range(1,97)}
    assert set(complete)==expected and len(expected)==210
    survivors=sorted(k for k,v in complete.items() if all(v.values()))
    assert survivors==[(1,4,4,0),(1,4,12,0),(2,2,1,82),(2,2,2,84)]
    assert rem==[[1,4,4,0],[1,4,12,0],[2,2,1,82]] and trem==[[2,84]]
    residual=read(NC/'M4N20_q29_residuals_root_arithmetic.txt')
    assert residual==[[2,1,0,29,0],[2,2,0,29,0]]+[[3,i,0,29,0] for i in range(1,5)]+[[4,4,0,29,1],[4,12,0,29,1],[2,1,82,29,1]]
    for k,v in complete.items():
        typ,d,i,j=k
        assert set(v)<=({29} if typ==1 and d==2 else {13,17})
        assert all(x in [0,1] for x in v.values())
    print('C420_SUPPORT_COVERAGE_PASS singles=18 pairs=192 prefix82+tail110=192 remaining=4',flush=True)
    spar,records=read(NC/'M4N20_stabilizers_root_arithmetic.txt');assert list(map(tuple,spar))==pars
    assert {tuple(r[:3])+(0 if r[0]==1 else r[3],) for r in records}==set(survivors)
    Hs=[]
    for typ,d,i,j,key,images,st in records:
        fix=[k+1 for k,x in enumerate(images) if x==key];assert fix==st
        H={pars[k-1] for k in st};assert all(groupmul(x,y) in H for x in H for y in H)
        if typ==1:
            assert len(H)==8 and len(set(map(tuple,images)))==12
            assert images==a[d-2][4][i-1] and key==a[d-2][2][j-1]
        else:
            assert len(H)==4
            assert all(groupmul(x,x)==(1,0,0,1,1) for x in H)
            keys=a[0][2];poly=a[0][1];F=Field(poly)
            # Each pair image must be obtained from its two independently checked input-place orbits.
            ownindex=a[0][5][i-1];assert key==sorted([keys[ownindex-1],keys[j-1]])
            repimages=a[0][4][i-1];partnerimages=a[0][4][a[0][3][j-1]-1]
            h=pars[partnerimages.index(keys[j-1])]
            expectedims=[sorted([repimages[k],partnerimages[pars.index(groupmul(g,h))]]) for k,g in enumerate(pars)]
            assert images==expectedims
            # Every partner action is recovered by the exact group law from
            # independently checked representative actions, rather than trusted.
            assert len(set(map(freeze,images)))==24
            assert all(len(im)==2 and all(x in keys for x in im) for im in images)
            Hs.append((i,H))
    H1=dict(Hs)[1];H2=dict(Hs)[2]
    assert H1=={(1,0,0,1,1),(3,2,2,3,1),(0,3,1,0,2),(2,1,3,2,2)}
    T=(1,0,1,1,1);Ti=(1,0,3,1,1)
    assert {groupmul(groupmul(T,h),Ti) for h in H1}==H2
    print('C420_STABILIZERS_PASS quartic=8,8 pair=4,4 conjugacy=SL2',flush=True)
    report={'complete':True,'NC_counts':counts,'Aut96_orbit_sizes':orbit_sizes,'cusp_audit':cusp,
            'single_supports':18,'quadratic_pairs_including_repetition':192,'witness_records':214,
            'remaining_before_global_completion':[list(k) for k in survivors],'global_completion':'C420_compact_manifest.json',
            'stabilizer_orders':[8,8,4,4],'status':'All4 remaining cases excluded by non-CM stabilizers plus compact Qi quotient; CM excluded by CCRS.'}
    print('C420_FULL_FINITE_REPLAY_DONE complete=true',flush=True)
    return report
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--report',type=Path);args=ap.parse_args();report=main()
    if args.report:args.report.write_text(json.dumps(report,indent=2)+'\n')
