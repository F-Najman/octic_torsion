#!/usr/bin/env python3
"""Independent integer-lattice replay of the full-level10 two-prime sieve."""
import ast, hashlib, json
from math import prod
from pathlib import Path
from sympy import Matrix, ilcm
from sympy.matrices.normalforms import hermite_normal_form
P=Path(__file__).resolve().parent
f=P/'data'/'diagonal_full10_D2_root_arithmetic.txt'
sigs,invs,rels,cusps,source,target,matches=ast.literal_eval(f.read_text().replace('[*','[').replace('*]',']').replace('false','False').replace('true','True'))
assert len(sigs)==len(set(map(tuple,sigs)))==36
assert len(source)==150 and len(target)==810 and matches==[]
def hnf(rows,n):
 h=hermite_normal_form(Matrix(rows).T)
 assert h.shape==(n,n) and h.det()!=0
 return h
def keyfun(h):
 hi=h.inv(); den=int(ilcm(*[a.q for a in hi])); m=[[int(den*a) for a in row] for row in hi.tolist()]
 return lambda v:tuple(sum(a*b for a,b in zip(row,v))%den for row in m)
r=[];n=[];orders=[];c=[]
for inv,rel,cs in zip(invs,rels,cusps):
 assert inv[-1]==0
 nn=len(inv)-1;n.append(nn)
 assert all(len(v)==nn+1 and v[-1]==0 for v in rel+cs)
 rr=[v[:-1] for v in rel];cc=[v[:-1] for v in cs];r.append(rr);c.append(cc)
 order=abs(int(hnf(rr,nn).det()));orders.append(order)
 assert order==prod(inv[:-1])
 assert order//abs(int(hnf(rr+cc,nn).det()))==prod([2,6,30,30,30,30,30,30,30,30])
rows=[a+[0]*n[1] for a in r[0]]+[[0]*n[0]+a for a in r[1]]+[a+b for a,b in zip(*c)]
h=hnf(rows,sum(n));diag=keyfun(h)
horder=orders[0]*orders[1]//abs(int(h.det()))
assert horder==prod([2,6,30,30,30,30,30,30,30,30])
assert orders[0]//horder==81
# The diagonal subgroup has the same order as both cusp images. Hence both
# projections are isomorphisms, independently certifying equal cusp kernels.
q11=keyfun(hnf(r[0]+c[0],n[0]));k31=keyfun(hnf(r[1],n[1]))
lookup={}
for i,cl,pred in source:
 assert 1<=i<=150 and len(cl)==n[0]+1 and len(pred)==n[1]+1 and cl[-1]==pred[-1]==0
 assert any(q11(cl[:-1])) and not any(q11([3*a for a in cl[:-1]]))
 assert not any(diag([3*a for a in cl[:-1]]+pred[:-1]))
 lookup.setdefault(diag([3*a for a in cl[:-1]]+[0]*n[1]),[]).append(i)
found=[]
for j,cl,hit in target:
 assert 1<=j<=810 and len(cl)==n[1]+1 and cl[-1]==0
 ii=lookup.get(diag([0]*n[0]+[-3*a for a in cl[:-1]]),[])
 assert bool(ii)==hit
 found.extend((i,j) for i in ii)
assert found==[]
report={'degree':8,'group':[10,10],'relative_degree':2,'source_places':150,'target_places':810,'cusp_group_order':horder,'source_quotient_order':81,'cusp_transport_kernels_equal':True,'multiplier':3,'matches':0,'arithmetic_sha256':hashlib.sha256(f.read_bytes()).hexdigest()}
(P/'data'/'replay_full10_D2_root.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2));print('FULL10_D2_ROOT_INDEPENDENT_HNF_REPLAY_DONE')
