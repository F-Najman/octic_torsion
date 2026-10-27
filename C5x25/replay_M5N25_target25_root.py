"""Independent binary fields, complete elliptic scan, and cusp-space ranks."""
from pathlib import Path
import ast, hashlib, json, re
P=Path(__file__).resolve().parent/'data'

def read_magma(name):
    s=(P/name).read_text().replace('[*','[').replace('*]',']')
    return ast.literal_eval(s)

class BinaryField:
    def __init__(self, polynomial):
        self.d=len(polynomial)-1;self.q=1<<self.d
        self.pol=sum(c<<i for i,c in enumerate(polynomial))
        def raw(a,b):
            out=0
            while b:
                if b&1:out^=a
                b>>=1;a<<=1
                if a&self.q:a^=self.pol
            return out
        self.table=[[raw(a,b) for b in range(self.q)] for a in range(self.q)]
        self.inv=[0]+[next(b for b in range(1,self.q) if raw(a,b)==1) for a in range(1,self.q)]
        self.square=[raw(a,a) for a in range(self.q)]
        self.tr=[]
        for a in range(self.q):
            t=0;b=a
            for _ in range(self.d):t^=b;b=self.square[b]
            assert t in (0,1);self.tr.append(t)
        self.as_roots={}
        for a in range(self.q):self.as_roots.setdefault(self.square[a]^a,[]).append(a)
        assert len(self.as_roots)==self.q//2
    def mul(self,a,b):return self.table[a][b]
    def div(self,a,b):assert b;return self.mul(a,self.inv[b])
    def pow(self,a,n):
        out=1
        while n:
            if n&1:out=self.mul(out,a)
            a=self.square[a];n>>=1
        return out
    def bits(self,a):return tuple((a>>i)&1 for i in range(self.d))
    def key16(self,xy):
        return min(tuple(v for a in xy for v in self.bits(a)),
                   tuple(v for a in xy for v in self.bits(self.pow(a,16))))

class Elliptic:
    def __init__(self,F,a):self.F=F;self.a1,self.a2,self.a3,self.a4,self.a6=a
    def on(self,P):
        if P is None:return True
        F=self.F;m=F.mul;x,y=P
        return F.square[y]^m(self.a1,m(x,y))^m(self.a3,y)==m(x,F.square[x])^m(self.a2,F.square[x])^m(self.a4,x)^self.a6
    def add(self,P,Q):
        if P is None:return Q
        if Q is None:return P
        F=self.F;m=F.mul;x,y=P;u,v=Q
        if x==u:
            if y^v^m(self.a1,x)^self.a3==0:return None
            assert P==Q
            lam=F.div(F.square[x]^self.a4^m(self.a1,y),m(self.a1,x)^self.a3)
        else:lam=F.div(y^v,x^u)
        nu=y^m(lam,x);xx=F.square[lam]^m(self.a1,lam)^self.a2^x^u
        yy=m(lam^self.a1,xx)^nu^self.a3
        out=(xx,yy);assert self.on(out);return out
    def multiple(self,n,P):
        out=None
        while n:
            if n&1:out=self.add(out,P)
            P=self.add(P,P);n>>=1
        return out
    def tate(self,P):
        F=self.F;m=F.mul;x,y=P
        aa3=self.a3^m(self.a1,x);aa2=x^self.a2
        aa4=F.square[x]^self.a4^m(self.a1,y)
        c0=F.div(aa4,aa3);a2=aa2^F.square[c0]^m(self.a1,c0)
        b=F.div(m(a2,F.square[a2]),F.square[aa3]);c=F.div(m(self.a1,a2)^aa3,aa3)
        r=F.div(b,c);s=F.div(F.square[c],b^c);t=m(r,s)^1
        xy=(F.div(s^r,t),F.div(t,F.square[s]^s^r^1))
        return (b,c),xy

def rank(F,rows):
    piv={}
    for row in rows:
        row=list(row)
        for j in sorted(piv):
            a=row[j]
            if a:row=[x^F.mul(a,y) for x,y in zip(row,piv[j])]
        if any(row):
            j=next(i for i,v in enumerate(row) if v)
            row=[F.div(v,row[j]) for v in row];piv[j]=row
    return len(piv)

f,records,reps,tests=read_magma('M5N25_target25_good_arithmetic_root.txt')
F=BinaryField(f);assert F.d==8 and len(records)==100 and len(reps)==12
native_keys=set()
for row in records:
    assert len(row)==11
    E=Elliptic(F,row[:5]);Q=tuple(row[5:7]);assert E.on(Q)
    assert E.multiple(25,Q) is None and E.multiple(5,Q) is not None
    bc,xy=E.tate(Q);assert bc==tuple(row[7:9]) and xy==tuple(row[9:11])
    assert F.pow(xy[0],16)!=xy[0] or F.pow(xy[1],16)!=xy[1]
    key=F.key16(xy);assert key not in native_keys;native_keys.add(key)

# An independent exhaustive scan of both twists of every ordinary j.
# E: y^2+xy=x^3+a2*x^2+a6, j=1/a6. The two trace classes
# of a2 are the two F256 twists. Trace tests count points exactly.
twist=next(a for a in range(F.q) if F.tr[a])
scan_keys=set();curves=0;qualifying_j=[]
for j in range(1,F.q):
    a6=F.inv[j]
    for a2 in (0,twist):
        count=2+2*sum(F.tr[x^a2^F.div(a6,F.square[x])]==0 for x in range(1,F.q))
        if count!=250:continue
        E=Elliptic(F,[1,a2,0,0,a6]);pts=[]
        y0=F.pow(a6,128);assert F.square[y0]==a6;pts.append((0,y0))
        for x in range(1,F.q):
            rhs=x^a2^F.div(a6,F.square[x])
            for z in F.as_roots.get(rhs,[]):pts.append((x,F.mul(x,z)))
        assert len(pts)+1==count and all(E.on(Q) for Q in pts)
        five=[Q for Q in pts if E.multiple(5,Q) is None]
        if len(five)!=24:continue
        curves+=1;qualifying_j.append(j)
        twentyfive=[Q for Q in pts if E.multiple(25,Q) is None and E.multiple(5,Q) is not None]
        assert len(twentyfive)==100
        for Q in twentyfive:scan_keys.add(F.key16(E.tate(Q)[1]))
assert curves==4 and scan_keys==native_keys
seen=set();orbit_sizes=[]
for row in reps:
    E=Elliptic(F,row[:5]);Q=tuple(row[5:7])
    orbit={F.key16(E.tate(E.multiple(d,Q))[1]) for d in range(1,13) if d%5}
    assert orbit<=native_keys and not orbit&seen
    seen|=orbit;orbit_sizes.append(len(orbit))
assert seen==native_keys
assert len(tests)==12 and [v[0] for v in tests]==list(range(1,13))
assert all(d==2 and hd==24 and sum(ds)>=1 and principal==0 for _,d,hd,ds,principal in tests)
print('M5N25_INDEPENDENT_ELLIPTIC_SCAN_PASS',curves,len(scan_keys),orbit_sizes,flush=True)

f16,spaces=read_magma('M5N25_target25_cusp_spaces_root.txt')
G=BinaryField(f16);assert G.d==4 and len(spaces)==18
assert all(rank(G,v[0])==11 and rank(G,v[1])==10 for v in spaces)
assert all(rank(G,spaces[i][0]+spaces[j][0])==12 for i in range(18) for j in range(i+1,18))
print('M5N25_INDEPENDENT_CUSP_SPACE_PASS degree2_divisors=171 h0=1',flush=True)

log=(P.parent/'logs'/'rank_J25_Qzeta5_root.log').read_text()
rows=re.findall(r'level=(\d+) order=(\d+) conductor=(\d+) dimension=(\d+) winding=(\d+) coefficient_degree=(\d+) eigen6_nontrivial=(true|false)',log)
assert len(rows)==8
for N,o,c,d,w,e,nontriv in rows:
    N,o,c,d,w,e=map(int,(N,o,c,d,w,e))
    if N==25:assert d==({1:0,2:0,5:1,10:2}[o])
    elif o in (5,10):assert d==w==({5:7,10:8}[o])
    if o in (5,10):assert c==25 and e==4 and nontriv=='true'
assert 'J25_CYCLO5_SPLIT11_ODD_ORDER 5014032203125' in log
assert 'J25_CYCLO5_FINITE16_TOTAL_POINTS 18' in log
assert 'J25_CYCLO5_RANK_SUPPORT_COMPLETE' in log
names=['M5N25_target25_good_arithmetic_root.txt','M5N25_target25_cusp_spaces_root.txt','rank_J25_Qzeta5_root.log']
report={'complete':True,'independent_full_ordinary_j_twist_scan':True,'qualifying_twists':curves,
        'all_quadratic_places':100,'diamond_orbits':12,'orbit_sizes':orbit_sizes,
        'native_nonprincipal_hecke_tests':12,'cusp_divisors_with_h0_one':171,
        'rank_twist_support_checked':True,'split11_order_odd':True,
        'inputs':{n:hashlib.sha256((P.parent/'logs'/n if n.endswith('.log') else P/n).read_bytes()).hexdigest() for n in names}}
(P/'replay_M5N25_target25_root.json').write_text(json.dumps(report,indent=2)+'\n')
print('M5N25_TARGET25_INDEPENDENT_REPLAY_COMPLETE')
