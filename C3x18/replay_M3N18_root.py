#!/usr/bin/env python3
"""Independent finite, lattice, coverage, and proof-assembly checks for (3,18)."""
import ast, hashlib, importlib.util, itertools, json, re
from collections import Counter
from pathlib import Path
from sympy import Matrix, Poly, Symbol, ZZ
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form
import sys;sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'common'))
from M3N18_finite_arithmetic_root import FF, EC

P=Path(__file__).resolve().parent
hashes={}
def read(n):
    b=(P/'data'/n).read_bytes();hashes[n]=hashlib.sha256(b).hexdigest()
    return ast.literal_eval(b.decode().replace('[*','[').replace('*]',']'))

def load_group(n,p,expected):
    g=read(n);assert g[0]==p and g[3][-1]==0 and len(g[8])==21
    fp=g[8];tr=[];cl=[]
    for d in range(1,5):
        tr.append(-d*fp[20-d]-sum(fp[20-j]*tr[d-j-1] for j in range(1,d)))
        a=p**d+1-tr[-1]-sum(e*cl[e-1] for e in range(1,d) if d%e==0)
        assert a%d==0;cl.append(a//d)
    assert cl==expected
    assert sum(fp)==__import__('math').prod(x for x in g[3] if x)
    return g

g7=load_group('M3N18_9multiple_fresh_p7_group_data.m',7,[18,27,150,513])
g13=load_group('M3N18_9multiple_aut27_p13_group_data.m',13,[18,54,852,7506])
assert g7[1:3]==g13[1:3] and len(g7[1])==len(set(map(tuple,g7[1])))==24
assert Counter(g7[2])=={1:18,3:6} and all(len(s)==65 for s in g7[1])
H=hermite_normal_form(Matrix([v[:-1] for v in g7[4]+g7[5]]).T)
inv=smith_normal_form(H,domain=ZZ)
assert sorted(abs(inv[i,i]) for i in range(inv.rows) if abs(inv[i,i])!=1)==[3,3,9]
spec=importlib.util.spec_from_file_location('tate',P.parent/'common'/'cusp_conditions_relative.py')
tate=importlib.util.module_from_spec(spec);spec.loader.exec_module(tate)
for p in [7,13]:
    c=tate.compute(3,18,p,4)
    assert Counter(c['closed_cusp_degrees'])=={1:18,3:6}
    assert c['rational_residue_cusps_not_F_rational']==0
    assert c['allowed_residues_by_complement_degree']==[[1,7,13]]*3

pars=[(pow(d,-1,3),b,c,d) for b,c in itertools.product(range(3),repeat=2) for d in [1,5,7]]
assert len(set(pars))==27
for a,b,c,d in pars:
    for A,B,C,D in pars:
        v=((A*a+6*B*c)%3,(A*b+B*d)%3,(C*a+D*c)%3,(6*C*b+D*d)%18)
        if v[3] not in [1,5,7]:v=((-v[0])%3,(-v[1])%3,(-v[2])%3,(-v[3])%18)
        assert v in pars

pairings={7:set(),13:set()};checked=0
def images(p,d,f,qt):
    global checked
    assert Poly.from_list(f[::-1],Symbol('x'),modulus=p).is_irreducible
    F=FF(p,f);qt=tuple(map(tuple,qt));q,t=qt;E=EC(F,q,t)
    O=(F.z,F.z);Q=(F.neg(t),F.mul(t,t));R=E.mul(6,Q)
    assert E.mul(3,O) is None and E.mul(18,Q) is None
    assert E.mul(9,Q) is not None and E.mul(6,Q) is not None
    assert O not in [R,E.mul(2,R)]
    xx,yy=R
    sl=F.div(F.sub(F.sc(3,F.mul(xx,xx)),F.mul(E.a,yy)),F.add(F.sc(2,yy),F.add(F.mul(E.a,xx),E.b)))
    wp=F.div(F.neg(yy),F.sub(F.mul(sl,xx),yy))
    assert wp!=F.o and F.pow(wp,3)==F.o and not any(wp[1:])
    pairings[p].add(wp[0])
    conj=[];v=qt
    for _ in range(d):conj.append(v);v=tuple(F.pow(x,p) for x in v)
    assert len(set(conj))==d and v==qt
    out=[]
    for a,b,c,dd in pars:
        PP=E.add(E.mul(a,O),E.mul(6*b,Q));QQ=E.add(E.mul(c,O),E.mul(dd,Q))
        assert E.mul(9,QQ)==E.mul(9,Q)
        out.append(F.key(E.normal(PP,QQ)))
    checked+=1
    return set(out)

exports,tests,rem=read('sieve_3_18_D4_fullaut_root_v4_arithmetic.txt')
for row in exports:
    d,n,idx,labels,keys=row
    data=read(f'M3N18_9multiple_fresh_p7_degree{d}_data.m')
    assert n=={2:27,3:144,4:513}[d]==len(keys)==len(data[3])
    assert len(set(map(tuple,keys)))==n and len(labels)==n
    fs=data[2];F=FF(7,fs)
    assert [F.key(tuple(map(tuple,v))) for v in data[4]]==list(map(tuple,keys))
    seen={}
    for j,i in enumerate(idx,1):
        ims=images(7,d,fs,data[4][i-1])
        assert ims=={tuple(k) for k,l in zip(keys,labels) if l==j}
        for k in ims:assert k not in seen;seen[k]=j
    assert set(seen)==set(map(tuple,keys))
    assert [seen[tuple(k)] for k in keys]==labels
assert len(exports[0][2])==1 and len(exports[1][2])==8 and len(exports[2][2])==19
assert rem==[[3,4,0],[3,5,0],[3,7,0],[3,8,0]]+[[4,j,0] for j in [1,5,6,7,11,12,17,19]]
for j in [4,5,7,8]:assert Counter(exports[1][3])[j]==9
assert len([r for r in tests if r[0]==22])==27
assert all(r[-1]==0 for r in tests if r[0]==22)
assert not any(r[-1] for r in tests if r[0]==2)
by={}
for d,i,j,q,flag in tests:
    assert q in [13,19] and flag in [0,1]
    by.setdefault((d,i,j),[]).append((q,flag))
for k,v in by.items():
    assert v in [[(13,0)],[(13,1),(19,0)],[(13,1),(19,1)]]
    assert (v[-1][1]==1)==(list(k) in rem)
log=(P/'logs'/'sieve_3_18_D4_fullaut_root_v4.log').read_text()
assert 'complete true' in log and 'M3N18_ROOT_HECKE_GENERIC_REFERENCE 2 13 28' in log
assert not any(x in log for x in ['Runtime error','Assertion failed'])

data4=read('M3N18_9multiple_aut27_p13_degree4_data.m')
_,d,f,classes,coords,idx,keys,labels,orbs=data4
assert d==4 and len(keys)==7506 and len(idx)==len(coords)==len(orbs)==278
assert len(set(map(tuple,keys)))==7506
seen={};F=FF(13,f)
for j,(i,qt,orb) in enumerate(zip(idx,coords,orbs),1):
    assert F.key(tuple(map(tuple,qt)))==tuple(keys[i-1])
    ims=images(13,4,f,qt)
    assert len(ims)==27 and ims==set(map(tuple,orb))
    for k in ims:assert k not in seen;seen[k]=j
assert set(seen)==set(map(tuple,keys)) and [seen[tuple(k)] for k in keys]==labels

diag=json.loads((P/'data'/'M3N18_9multiple_fresh_diagonal_match_data_root.json').read_text())
matched=[diag['target'][i-1] for i,s in diag['matches']]
assert Counter(x[0] for x in matched)=={3:270,4:9}
sel,out,left=read('M3N18_p13_postdiagonal_arithmetic_root.txt');assert not left
assert sel==[sorted({r[1] for r in matched if r[0]==t}) for t in [3,4]]
assert len(sel[0])==54 and len(sel[1])==9
data3=read('M3N18_9multiple_aut27_p13_degree3_data.m')
F3=FF(13,data3[2]);seen3=set();checked4=set()
for d,i,j,orb,q,flag in out:
    assert flag==0
    if d==3:
        assert q==0 and orb==9 and j==i
        ims=images(13,3,data3[2],data3[4][i-1]);assert len(ims)==9
        assert not seen3.intersection(ims);seen3.update(ims)
    else:
        assert d==4 and q==7 and orb==27 and idx[i-1]==j
        checked4.add(i)
assert seen3=={F3.key(tuple(map(tuple,data3[4][i-1]))) for i in sel[0]}
assert checked4==set(sel[1])
postlog=(P/'logs'/'sieve_3_18_p13_postdiagonal_root.log').read_text()
assert 'M3N18_P13_POSTDIAG_GENERIC_REFERENCE 7 4 32' in postlog
assert 'M3N18_P13_POSTDIAG_COMPLETE []' in postlog
assert not any(x in postlog for x in ['Runtime error','Assertion failed'])
assert all(len(v)==1 for v in pairings.values())
report={'degree':8,'group':[3,18],'genus':10,'relative_degree':4,
 'source_NC_counts':[0,27,144,513],'target_NC_counts':[0,54,846,7506],
 'source_Hecke_orbits':[1,8,19],'source_quartics_after_Hecke':216,
 'target_quartic_orbits':278,'independent_elliptic_representatives':checked,
 'Aut_order':27,'cusp_subgroup_order':20253807,'quotient_exponent_bound':9,
 'matched_cubic_cases':270,'matched_cubic_supports':54,'cubic_orbits_stabilizer_excluded':6,
 'matched_quartic_orbits_Hecke_excluded':9,'remaining':0,
 'Weil_pairings':{str(p):next(iter(v)) for p,v in pairings.items()},
 'boundary':'Independent finite-field and elliptic arithmetic, all automorphism images/partitions, cusp Galois action, Newton counts, and proof assembly. The separate HNF comparator reconstructs all26225 states. Native proved class groups and exact IsPrincipal are not reimplemented in Python.',
 'input_sha256':hashes}
(P/'data'/'replay_M3N18_root.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2));print('M3N18_COMPLETE_EXCLUSION_REPLAY_DONE')
