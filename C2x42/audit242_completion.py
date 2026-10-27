#!/usr/bin/env python3
"""Read-only independent final C242 coverage/multiplicity/deck arithmetic audit."""
import argparse,ast,collections,hashlib,json,sys
from pathlib import Path
P=Path(__file__).resolve().parent;N=P/'data'
sys.path.insert(0,str(P.parent/'common'))
from finite_field5_checked import FF,EC,pack
bound={}
def rd(name):
    p=N/name;bound[name]=hashlib.sha256(p.read_bytes()).hexdigest()
    return ast.literal_eval(p.read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
def mul(E,Q,n):
    R=None
    while n:
        if n&1:R=E.add(R,Q)
        Q=E.add(Q,Q);n//=2
    return R
fields,rows=rd('M2N42_full_Aut_representative_supports.m');F={d:FF(f) for d,f in fields}
reps={i:(desc,sp) for i,desc,sp in rows};assert len(reps)==30
expectedtypes=collections.Counter(tuple(d for d,ky in sp) for desc,sp in reps.values())
assert expectedtypes=={(4,):2,(7,):1,(4,4):23,(8,):4}
sets=[];datasets=[rd('C242_g31_bounded_Hecke_v3_arithmetic.txt')]+[rd(f'C242_g31_remaining_shard{k}_arithmetic.txt') for k in [0,1]]
results={};allweights=[]
for sh,ds in enumerate(datasets):
    metadata={}
    for ell,d,key,e,w,ed,wd,outputs in ds[0]:
        assert d==e*w==ed*wd
        assert sum(k for k,ee,ww in outputs)==ell+1
        assert all(d*k==ee*ww and k>0 and ee>0 and ww>0 for k,ee,ww in outputs)
        assert sum(ee*ww for k,ee,ww in outputs)==d*(ell+1)
        ck=(ell,d,tuple(key));assert ck not in metadata;metadata[ck]=(e,w,ed,wd)
        allweights+=outputs
    got=set()
    for i,desc,ell,pr in ds[1]:
        assert i not in results and desc==reps[i][0];results[i]=(ell,pr);got.add(i)
        sp=reps[i][1];nc=sum(d for d,ky in sp);assert ell==(17 if nc==4 else 11)
        assert all((ell,d,tuple(ky)) in metadata for d,ky in sp)
    assert len(got)==(4 if sh==0 else 13);sets.append(got)
assert set(results)==set(reps) and all(not(a&b) for k,a in enumerate(sets) for b in sets[k+1:])
assert {i for i,(q,pr) in results.items() if pr}=={118,193}
assert sum(not pr for q,pr in results.values())==28
for sh in [0,1]:
    text=(P/'logs'/f'replay_C242_genus31_serialized_RR_shard{sh}.log').read_text()
    assert f'C242_G31_SERIALIZED_RR_REPLAY_DONE {sh} 13 complete=true' in text
ref=rd('C242_g31_pair118_refinement_arithmetic.txt');assert len(ref)==3
assert ref[1]==[[118,reps[118][0],13,True],[118,reps[118][0],17,False]]
metas={}
for q,d,key,e,w,ed,wd,outputs in ref[0]:
    assert d==e*w==ed*wd
    assert sum(k for k,ee,ww in outputs)==q+1 and all(d*k==ee*ww for k,ee,ww in outputs)
    assert sum(ee*ww for k,ee,ww in outputs)==d*(q+1)
    metas[q,d,tuple(key)]=True
for q in [13,17]:assert all((q,d,tuple(key)) in metas for d,key in reps[118][1])
assert len(ref[2])==1
log=(P/'logs'/'refine_C242_genus31_pair118_v2.log').read_text()
assert 'C242_G31_RESIDUAL_NECESSARY_RESULT 118 [* 2, 4, 3, 4, 85 *] 17 false' in log
assert 'C242_G31_PAIR118_REFINEMENT_DONE 2 complete=true' in log
desc,sp=reps[193];assert desc==[2,4,5,4,38] and [d for d,key in sp]==[4,4]
def original(sp):return collections.Counter((d,pack(key[:d]),pack(key[d:])) for d,key in sp)
def act(sp,p,r,u):
    ans=[]
    for d,key in sp:
        f=F[d];b,c=pack(key[:d]),pack(key[d:]);E=EC(f,b,c);Q=(b,b)
        pts=[(0,0),mul(E,Q,21)];pts.append(E.add(*pts));QQ=E.add(pts[r],mul(E,Q,2*u))
        ans.append((d,*f.key(E.bc(pts[p],QQ))))
    return collections.Counter(ans)
assert original(sp)==act(sp,2,1,10)
for d,key in sp:
    f=F[d];b,c=pack(key[:d]),pack(key[d:]);E=EC(f,b,c);Q=(b,b);T=mul(E,Q,21)
    assert T is not None and E.add(T,T) is None and T!=(0,0)
    # Parameter(2,1,10) gives Q'=-Q; simultaneous negation is the
    # literal deck transformation(P2,Q)->(P2+21Q,Q).
    assert E.add(T,mul(E,Q,20))==(Q[0],f.neg(Q[1]))
# Independent exact Tate-label check of the genus31 quotient and all cusp allowances.
from math import gcd,ceil
n=42;labels={(a,b,c,d) for a,b in __import__('itertools').product(range(n),repeat=2)
    if gcd(gcd(a,b),n)==1 for c,d in __import__('itertools').product(range(2),repeat=2) if (c*b-d*a)%2==1}
def orbit(v):
    a,b,c,d=v
    return frozenset((s*(a+k*b)%n,s*b%n,s*(c+k*d)%2,s*d%2) for s in [-1,1] for k in range(n))
cs=sorted(set(map(orbit,labels)),key=lambda s:sorted(s));ix={x:i for i,c in enumerate(cs) for x in c};assert len(cs)==72
def perm(fn):
    out=[]
    for c in cs:
        im={ix[fn(x)] for x in c};assert len(im)==1;out.append(im.pop())
    assert len(set(out))==72;return out
eps=perm(lambda x:(29*x[0]%42,29*x[1]%42,x[2],x[3]))
frob=perm(lambda x:(5*x[0]%42,x[1],x[2],x[3]))
assert all(eps[eps[i]]==i and eps[i]!=i for i in range(72))
assert (2*61-2)//2==2*31-2
left=set(range(72));cycles=[]
while left:
    i=min(left);cy=[]
    while i not in cy:cy.append(i);left.remove(i);i=frob[i]
    assert i==cy[0];cycles.append(cy)
assert collections.Counter(map(len,cycles))=={1:18,2:9,3:6,6:3}
for q,cd in [(11,1),(17,4)]:
    sig=perm(lambda x:(q*x[0]%42,x[1],x[2],x[3]))
    assert all(sig[c[0]] in c for c in cycles if len(c)<=cd)
assert ceil(325*1152/32768)==12
report={'complete':True,'representatives':30,'initial_nonprincipal':28,'pair118_killed_by':17,
        'pair193_deck_exact':True,'remaining':0,'initial_partition_sizes':list(map(len,sets)),
        'source_support_types':{str(k):v for k,v in expectedtypes.items()},
        'kernel_weight_records_checked':len(allweights),'quotient_genus':31,
        'source_geometric_gonality_lower_bound':12,'source_cusp_signature':dict(collections.Counter(map(len,cycles))),
        'boundary':'Native exact RR decisions are supplied by completed logs/serialized replays; global quartic42 absence is DNTheorem1.1.',
        'sha256':bound}
print(json.dumps(report,indent=2));print('C242_INDEPENDENT_COMPLETION_AUDIT_PASS remaining=0')
ap=argparse.ArgumentParser();ap.add_argument('--report',type=Path);args=ap.parse_args()
if args.report:args.report.write_text(json.dumps(report,indent=2)+'\n')
