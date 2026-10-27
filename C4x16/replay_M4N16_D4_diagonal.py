#!/usr/bin/env python3
"""Actual-relation integer replay of the completed F5/F9 M4 diagonal.

This verifies the exact surviving necessary conditions. It makes no
exclusion claim: the completed calculation has 17 surviving pairs.
"""
import ast,hashlib,json
from collections import Counter
from pathlib import Path
from sympy import Matrix,zeros
from sympy.matrices.normalforms import hermite_normal_form
P=Path(__file__).resolve().parent/'data'
name='diagonal_4_16_D4_p5_p3_arithmetic.txt';raw=(P/name).read_bytes()
a=ast.literal_eval(raw.decode().replace('[*','[').replace('*]',']'))
sigs,ds,inv5,r5,inv9,r9,cusp,source,target,matches,nc,reps,coverage,rat9,pc,counts,fps=a
assert Counter(ds)=={1:16,2:4,4:4} and len(sigs)==24 and len(set(map(tuple,sigs)))==24
assert nc==[0,0,320,1808] and reps==[5,33] and counts==[120,33]
assert len(source)==1184 and len(target)==153 and len(matches)==17
def tors(rows):
 assert all(r[-1]==0 for r in rows)
 return [r[:-1] for r in rows]
r5=tors(r5);r9=tors(r9);c5=tors(cusp[0]);c9=tors(cusp[1]);n5=len(inv5)-1;n9=len(inv9)-1
R5=hermite_normal_form(Matrix(r5).T);R9=hermite_normal_form(Matrix(r9).T)
H5=hermite_normal_form(Matrix(r5+c5).T);H9=hermite_normal_form(Matrix(r9+c9).T)
order5=int(abs(R5.det()));order9=int(abs(R9.det()));h5=order5//int(abs(H5.det()));h9=order9//int(abs(H9.det()))
assert order5==10485760000 and h5==h9==327680000 and abs(H5.det())==32
assert all(x.q==1 for x in (2*H5.inv())) # exponent(G5/C5)=2
def member(B,v):return all(x.q==1 for x in B.inv()*Matrix(v))
graphrows=[r+[0]*n9 for r in r5]+[[0]*n5+r for r in r9]+[u+v for u,v in zip(c5,c9)]
BG=hermite_normal_form(Matrix(graphrows).T);BGi=BG.inv()
assert order5*order9//int(abs(BG.det()))==h5 # both cusp projections inject
R5i=R5.inv();R9i=R9.inv()
def key(Bi,v):return tuple(x%1 for x in Bi*Matrix(v))
def sub(u,v):return [x-y for x,y in zip(u,v)]
sourcebase={};pred={};insideC5=0
rat=[i+1 for i,d in enumerate(ds) if d==1]
for k,(d,i,j,cl,pr) in enumerate(source,1):
 assert cl[-1]==pr[-1]==0
 cl=cl[:-1];pr=pr[:-1]
 assert d in [3,4] and (j in rat if d==3 else j==0)
 assert all(x.q==1 for x in BGi*Matrix([2*x for x in cl]+pr))
 base=cl if j==0 else sub(cl,c5[j-1]);bk=key(R5i,base)
 assert sourcebase.setdefault((d,i),bk)==bk
 insideC5+=member(H5,cl)
 pred.setdefault(key(R9i,pr),[]).append(k)
assert insideC5==0
assert {i for d,i in sourcebase if d==3}==set(range(1,65))
assert {i for d,i in sourcebase if d==4}==set(range(1,161))
assert Counter((d,i) for d,i,j,cl,pr in source)==Counter({**{(3,i):16 for i in range(1,65)},**{(4,i):1 for i in range(1,161)}})
got=[]
for k,(d,i,j,cl) in enumerate(target,1):
 assert cl[-1]==0
 if d==3:
  assert 1<=i<=5 and 1<=j<=24;want=[x+y for x,y in zip(pc[0][i-1],rat9[j-1])]
 else:
  assert d==4 and 1<=i<=33 and j==0;want=pc[1][i-1]
 assert key(R9i,cl[:-1])==key(R9i,want[:-1])
 ky=key(R9i,[2*x for x in cl[:-1]])
 if ky in pred:got.append([k,pred[ky]])
assert got==matches
assert Counter((d,i,j) for d,i,j,cl in target)==Counter({**{(3,i,j):1 for i in range(1,6) for j in range(1,25)},**{(4,i,0):1 for i in range(1,34)}})
for d,keys,labels,idx in coverage:
 assert len(keys)==nc[d-1] and len(keys)==len(labels) and len(idx)==reps[d-3]
 assert len(set(map(tuple,keys)))==len(keys)
 assert all(len(k)==4*d and set(k)<={0,1,2} for k in keys)
 assert set(labels)==set(range(1,len(idx)+1)) and all(labels[i-1]==j for j,i in enumerate(idx,1))
 assert all(v<=64 and 64%v==0 for v in Counter(labels).values())
report={'group':[4,16],'relative_degree':4,'fields':[5,9],'cusp_group_order':h5,
 'finite_J5_order':order5,'finite_J9_order':order9,'quotient_J5_C_exponent':2,
 'cusp_graph_has_injective_projections':True,'source_divisors':1184,'source_classes_in_C':insideC5,
 'target_patterns':{'cubic_plus_cusp':120,'quartic':33},'matching_pairs':matches,
 'exclusion_proved':False,'input_sha256':hashlib.sha256(raw).hexdigest()}
(P/'replay_M4N16_D4_diagonal.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2));print('M4N16_D4_DIAGONAL_INTEGER_REPLAY_DONE')
