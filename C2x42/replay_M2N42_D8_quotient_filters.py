#!/usr/bin/env python3
"""Replay all D8 coverage, actual-relation arithmetic, and residue multipliers.

Magma supplies proved finite Jacobian presentations and exact divisor classes.
Python independently reconstructs all support cases, all integer-lattice
membership decisions, and every moduli-map residue degree from finite-field
elliptic-curve arithmetic. This is a necessary-condition filter, not exclusion.
"""
import ast, hashlib, json, time
from collections import Counter
from pathlib import Path
from sympy import Matrix, ilcm
from sympy.matrices.normalforms import hermite_normal_form
import sys;sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'common'))
from finite_field5_checked import FF, EC, pack

P=Path(__file__).resolve().parent; start=time.monotonic(); bound={}
def read(name):
    raw=(P/'data'/name).read_bytes(); bound[name]=hashlib.sha256(raw).hexdigest()
    return ast.literal_eval(raw.decode().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))

class FiniteGroup:
    def __init__(self, inv, relations, hg=()):
        assert inv[-1]==0 and all(v[-1]==0 for v in list(relations)+list(hg))
        self.inv=inv; self.n=len(inv)-1
        self.R=hermite_normal_form(Matrix([r[:-1] for r in relations]).T)
        self.B=hermite_normal_form(Matrix([r[:-1] for r in list(relations)+list(hg)]).T)
        self.index=abs(int(self.B.det()))
        self.inverse=self.B.inv()
        self.den=int(ilcm(*[a.q for a in self.inverse])) if self.n>1 else int(self.inverse[0].q)
        self.M=(self.den*self.inverse).applyfunc(int)
    def key(self,v):
        assert len(v)==self.n+1 and v[-1]==0
        return tuple(int(a)%self.den for a in self.M*Matrix(v[:-1]))
    def zero(self,v): return not any(self.key(v))
    def equal(self,a,b):
        v=self.R.inv()*Matrix([a[i]-b[i] for i in range(self.n)])
        return a[-1]==b[-1]==0 and all(x.q==1 for x in v)

def add(*vs): return [sum(x) for x in zip(*vs)]
def mul(E,P,n):
    A=None
    while n:
        if n&1:A=E.add(A,P)
        P=E.add(P,P);n//=2
    return A
def tatexy(F,c,d,Q):
    x,y=Q; a3=F.sc(2,y); a2=F.add(F.sc(3,x),c)
    a4=F.add(F.add(F.sc(3,F.mul(x,x)),F.sc(2,F.mul(c,x))),d)
    slope=F.div(a4,a3); a1=F.sc(2,slope); a2=F.sub(a2,F.mul(slope,slope))
    b=F.neg(F.div(F.pow(a2,3),F.mul(a3,a3)))
    cc=F.neg(F.div(F.sub(F.mul(a1,a2),a3),a3))
    r=F.div(b,cc);s=F.div(F.mul(cc,cc),F.sub(b,cc))
    tt=F.add(F.sub(F.mul(r,s),F.sc(2,r)),1)
    return F.div(F.sub(s,r),tt),F.div(tt,F.add(F.sub(F.sub(F.mul(s,s),s),r),1))
def isog2(F,E,Q,r):
    A=F.add(E.c,F.sc(3,r));B=F.add(F.add(F.sc(3,F.mul(r,r)),F.sc(2,F.mul(E.c,r))),E.d)
    x=F.sub(Q[0],r);xx=F.mul(x,x)
    return F.neg(F.sc(2,A)),F.sub(F.mul(A,A),F.sc(4,B)),(
        F.div(F.mul(Q[1],Q[1]),xx),F.div(F.mul(Q[1],F.sub(B,xx)),xx))

shards=[read(f'filter_2_42_D8_small_quotients_shard{s}_arithmetic.txt') for s in range(4)]
gd,fields=shards[0][:2]
assert all(a[:2]==[gd,fields] for a in shards)
gy,gz=gd; Y=FiniteGroup(gy[0],gy[1],gy[3]); Z=FiniteGroup(gz[0],gz[1],gz[3])
assert gy[0]==[2,2,12,36,0] and gy[4]==[2,2,6,18] and Y.index==4
assert gz[0]==[2184,0] and gz[4]==[364] and Z.index==6
meta=read('M2N42_small_targets_full_group_shard3.txt');assert meta[0]==gd
csdegrees,csclasses,qinv,qrelations,qimages=meta[1:]
QZ=FiniteGroup(qinv,qrelations)
assert QZ.index==6
def qmap(v):return [sum(v[i]*qimages[i][j] for i in range(len(v))) for j in range(len(qinv))]
assert all(QZ.zero(qmap(v)) for v in gz[1]+gz[3])
assert QZ.key(qmap([1,0]))!=(0,)
assert len({QZ.key(qmap([i,0])) for i in range(2184)})==6
compG=[{(0,0)}]
for d in range(1,5):
    ss=set()
    for e,v in zip(csdegrees,csclasses):
        if e<=d:ss.update(tuple(add(v,list(w))) for w in compG[d-e])
    compG.append(ss)
compkeys=[{Z.key(v) for v in ss} for ss in compG]
assert list(map(len,compkeys))==[1,1,1,2,2]
for d in range(5):
    assert {QZ.key(qmap(v)) for v in compG[d]}=={QZ.key(v) for v in gz[5][d]}

masters={}; expected={}; Fdict={d:FF(f) for d,f in fields}
for D in range(4,9):
    a=read(f'M2N42_D8_direct_full{D}_arithmetic.txt');assert a[0]==D and a[2] and a[1]==dict(fields)[D]
    masters[D]=a[8];expected[D]={i:r[1] for i,r in enumerate(a[8],1) if r[0]==D}
assert [len(expected[D]) for D in range(4,9)]==[5,18,64,316,1372]
source_cache={}
def source(D,key):
    ck=(D,tuple(key))
    if ck not in source_cache:
        F=Fdict[D];b,c=pack(key[:D]),pack(key[D:]);E=EC(F,b,c);Q=(b,b)
        source_cache[ck]=(F,E,Q)
    return source_cache[ck]
def check_first_weight(D,key,ed):
    F,E,Q=source(D,key);ey,wy,ez,wz=ed
    assert F.degree(E.bc((0,0),mul(E,Q,3)))==ey and D==ey*wy
    assert F.degree(tatexy(F,E.c,E.d,mul(E,Q,2)))==ez and D==ez*wz

singles=[];small=[];pairs=[];inputs=[];seen={D:set() for D in range(4,9)}
own={D:[] for D in range(4,9)};passed=Counter()
for sh,a in enumerate(shards):
    _,_,ss,sm,pp,survs=a;want=[]
    for row in ss:
        D,i,key,ed,cy,cz,flag=row
        assert expected[D][i]==key and i not in seen[D]
        assert (D<8 and sh==0) or (D==8 and (i-1)%4==sh)
        seen[D].add(i);check_first_weight(D,key,ed)
        yes=Y.zero(cy) and tuple(-x%Z.den for x in Z.key(cz)) in compkeys[8-D]
        assert bool(flag)==yes
        own[D].append(row);singles.append(row)
        if yes:want.append([1,D,i]);passed[D]+=1
    if sh==0:
        small=sm;pairs=pp
        all4=read('M2N42_D8_direct_small_all_points.txt');assert len(all4)==1 and all4[0][0]==4
        assert len(sm)==126
        for i,row in enumerate(sm,1):
            D,j,key,ed,cy,cz=row;assert D==4 and j==i and key==sum(all4[0][2][i-1],[])
            check_first_weight(D,key,ed)
        assert len(pp)==630
        for idx,row in enumerate(pp):
            i,j,cy,cz,flag=row;assert (i,j)==(idx//126+1,idx%126+1)
            assert Y.equal(cy,add(own[4][i-1][4],sm[j-1][4]))
            assert Z.equal(cz,add(own[4][i-1][5],sm[j-1][5]))
            yes=Y.zero(cy) and Z.zero(cz);assert bool(flag)==yes
            if yes:want.append([2,4,i,4,j])
    else:assert sm==[] and pp==[]
    assert survs==want;inputs+=survs
    log=(P/'logs'/f'filter_2_42_D8_small_quotients_shard{sh}.log').read_text()
    assert f'M2N42_SMALL_QUOTIENT_FILTER_DONE {sh}' in log and 'complete=true' in log
for D in range(4,9):assert seen[D]==set(expected[D])
assert len(singles)==1775 and len(pairs)==630 and len(inputs)==288
assert [passed[D] for D in range(4,9)]==[3,0,6,15,78]
assert sum(r[-1] for r in pairs)==186
def supports(desc):
    if desc[0]==1:
        D,i=desc[1:];return [(D,masters[D][i-1][1])]
    assert desc[0]==2 and desc[1]==desc[3]==4
    return [(4,own[4][desc[2]-1][2]),(4,small[desc[4]-1][2])]

second=[read(f'filter_2_42_D8_three_2isogenies_v2_shard{s}_arithmetic.txt') for s in range(4)]
gd2=second[0][0];assert all(a[0]==gd2 and a[1]==fields for a in second)
Z2=FiniteGroup(gd2[0],gd2[1],gd2[3]);assert Z2.index==6 and gd2[4]==[364]
assert gd2[0:2]==gz[0:2] and gd2[4]==gz[4]
assert Counter(gd2[2])==Counter(gz[2]) and Counter(gd2[5])==Counter(csdegrees)
# Fresh class-group calls can use a different generator or rational base cusp.
# Recompute each call's own complement classes, without identifying its basis.
comps2=[{(0,0)}]
for dd in range(1,5):
    ss=set()
    for ee,v in zip(gd2[5],gd2[6]):
        if ee<=dd:ss.update(tuple(add(v,list(w))) for w in comps2[dd-ee])
    comps2.append(ss)
compkeys2=[{Z2.key(v) for v in ss} for ss in comps2]
assert list(map(len,compkeys2))==[1,1,1,2,2]
assert [{Z2.key(v) for v in ss} for ss in gd2[7][3]]==compkeys2
images={};tests={};remaining=[];types=Counter()
for sh,a in enumerate(second):
    for mapno,D,key,e,w,cl in a[2]:
        F,E,Q=source(D,key);T=mul(E,Q,21)
        r=0 if mapno==1 else T[0] if mapno==2 else F.sub(F.neg(E.c),T[0])
        cp,dp,A=isog2(F,E,mul(E,Q,2),r)
        assert F.degree(tatexy(F,cp,dp,A))==e and D==e*w
        ck=(mapno,D,tuple(key))
        if ck in images:assert Z2.equal(images[ck],cl)
        else:images[ck]=cl
    for row in a[3]:
        k,desc,checks,flag=row;assert k not in tests and (k-1)%4==sh and desc==inputs[k-1]
        tests[k]=row;sp=supports(desc);cd=8-sum(d for d,key in sp);rr=[]
        assert cd in [0,1,2,4]
        for mapno,cl,hit in checks:
            assert mapno==len(rr)+1
            computed=add(*[images[mapno,D,tuple(key)] for D,key in sp])
            assert Z2.equal(cl,computed)
            yes=tuple(-x%Z2.den for x in Z2.key(cl)) in compkeys2[cd]
            assert bool(hit)==yes;rr.append(yes)
            if not yes:assert len(rr)==len(checks)
        assert bool(flag)==(len(rr)==3 and all(rr))
        if flag:remaining.append((k,desc));types[tuple(d for d,key in sp)]+=1
    assert a[4]==[r[1] for r in a[3] if r[-1]]
    log=(P/'logs'/f'filter_2_42_D8_three_2isogenies_v2_shard{sh}.log').read_text()
    assert f'M2N42_ISOG21_FILTER_DONE {sh} 72' in log and 'complete=true' in log
assert set(tests)==set(range(1,289)) and len(remaining)==41
remaining.sort()
out={'status':'complete necessary-condition filters; NOT an exclusion',
 'source_group':[2,42],'source_genus':61,'prime':5,
 'source_single_tests':1775,'source_quartic_pairs_with_repetition':630,
 'target_genera':[4,5],'first_target_cusp_image_indices':[4,6],
 'finite_cusp_complement_cosets_degree0_to4':[1,1,1,2,2],
 'first_joint_filter_survivors':288,'three_two_isogeny_filter_survivors':41,
 'remaining_support_types':{str(k):v for k,v in types.items()},
 'remaining':remaining,'all_five_map_residue_degree_weights_independently_checked':True,
 'all_actual_relation_class_additions_and_memberships_checked':True,
 'verification_boundary':__doc__,'sha256':bound,'seconds':round(time.monotonic()-start,3)}
(P/'data'/'replay_M2N42_D8_quotient_filters.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2));print('M2N42_D8_ALL_QUOTIENT_FILTERS_REPLAY_DONE')
