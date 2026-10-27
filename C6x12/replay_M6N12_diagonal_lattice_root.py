#!/usr/bin/env python3
"""Independently rebuild the112020-case diagonal and global unit reduction."""
from pathlib import Path
import ast,json
from collections import Counter,defaultdict
from math import prod
from sympy import Matrix,ilcm
from sympy.matrices.normalforms import hermite_normal_form
P=Path(__file__).resolve().parent
def read(n):return ast.literal_eval((P/'data'/n).read_text().replace('[*','[').replace('*]',']'))
def basis(rows,n):
    H=hermite_normal_form(Matrix(rows).T);assert H.shape==(n,n) and H.det();return H
def keymap(H):
    I=H.inv();den=int(ilcm(*[x.q for x in I]));A=[[int(den*x) for x in row] for row in I.tolist()]
    return lambda v:tuple(sum(x*y for x,y in zip(row,v))%den for row in A)
sigs,ds,invs,rels,cusps,source,records,matches,counts=read('diagonal_6_12_D4_root_arithmetic.txt')
assert len(sigs)==24 and len(set(map(tuple,sigs)))==24 and Counter(ds)=={1:16,2:8}
assert len(sigs[0])==141 and all(v[-1]==0 for R in rels+cusps for v in R)
assert invs[0][-1]==invs[1][-1]==0
n,m=[len(v)-1 for v in invs];r7,r13=rels;c7,c13=cusps
R7=basis([v[:-1] for v in r7],n);R13=basis([v[:-1] for v in r13],m)
k7,k13=keymap(R7),keymap(R13)
assert abs(int(R7.det()))==95551488==prod(x for x in invs[0] if x)
H7=basis([v[:-1] for v in r7+c7],n);H13=basis([v[:-1] for v in r13+c13],m)
assert abs(int(R7.det()/H7.det()))==abs(int(R13.det()/H13.det()))==5971968
assert abs(int(H7.det()))==16
DG=basis([v[:-1]+[0]*m for v in r7]+[[0]*n+v[:-1] for v in r13]+[v[:-1]+w[:-1] for v,w in zip(c7,c13)],n+m)
assert abs(int(R7.det()*R13.det()/DG.det()))==5971968
kd=keymap(DG);qh=keymap(H7);zero=(0,)*n
pred=defaultdict(list)
for i,r,cl,a in source:
    assert i in range(1,469) and r in [2,4] and cl[-1]==a[-1]==0
    assert qh([r*x for x in cl[:-1]])==zero
    if r==4:assert qh([2*x for x in cl[:-1]])!=zero
    assert not any(kd([r*x for x in cl[:-1]]+a[:-1]))
    pred[(r,k13(a[:-1]))].append(i)
assert [r[0] for r in source]==list(range(1,469))
found=[];types=Counter();descrs=set()
for descr,cl,hits in records:
    typ,i,j=descr;assert cl[-1]==0;descrs.add(tuple(descr));types[typ]+=1
    if typ==2:assert 1<=i<=144 and 1<=j<=528
    elif typ==3:assert 1<=i<=576 and 1<=j<=32
    elif typ==4:assert 1<=i<=7116 and j==1
    else:assert typ==22 and 1<=i<=j<=144
    exact=[]
    for r in [2,4]:exact+=pred.get((r,k13([r*x for x in cl[:-1]])),[])
    assert exact==hits
    if hits:found.append([descr,hits,cl])
assert len(descrs)==len(records)==112020 and counts==[76032,18432,7116,10440]
assert types=={2:76032,3:18432,4:7116,22:10440} and found==matches and len(found)==24
ids=sorted({i for r in found for i in r[1]});assert ids==list(range(127,133))+list(range(343,349))
assert all(source[i-1][1]==2 for i in ids)

ss,dd,lifts=read('M6N12_global_L8_cusp_classes_root.m');assert ss==sigs and dd==ds
assert [r[0] for r in lifts]==ids
U=Matrix(sigs);assert all(sum(ds[i]*U[i,j] for i in range(24))==0 for j in range(U.cols))
H=hermite_normal_form(U);assert H.shape==(24,23)
base=ds.index(1);e=Matrix([int(i==base) for i in range(24)])
assert abs(int(hermite_normal_form(U.row_join(e)).det()))==5971968
short=read('M6N12_reduced_L8_cusp_coefficients_magma_root.m');assert [r[0] for r in short]==ids
for (i,a),(j,v) in zip(lifts,short):
    assert i==j and len(a)==len(v)==24
    cl=source[i-1][2][:-1]
    mapped=[sum(a[t]*c7[t][k] for t in range(24)) for k in range(n)]
    assert not any(k7([x-2*y for x,y in zip(mapped,cl)]))
    old=a[:];old[base]+=8-sum(x*d for x,d in zip(old,ds))
    assert all(x>=0 for x in v) and sum(x*d for x,d in zip(v,ds))==8
    sol=H.gauss_jordan_solve(Matrix([x-y for x,y in zip(old,v)]))[0]
    assert all(x.q==1 for x in sol)
assert len({tuple(v) for i,v in short})==6
assert all(v in [w for i,w in short[:6]] for j,v in short)
assert len({tuple(s[j-1] for j in [34,58,27]) for s in sigs})==24
rows=read('M6N12_L8_char0_six_pencils_arithmetic_root.txt')
assert [r[1] for r in rows]==list(range(127,133))
assert all(r==[0,r[1],2,8,[8,8,8,8],[1,1,1,1]] for r in rows)
log=(P/'logs'/'RR_6_12_L8_char0_six_pencils_root.log').read_text()
assert 'M6N12_L8_FINITE_RR_COMPLETE' in log and not any(x in log for x in ['Runtime error','Assertion failed'])
report={'group':[6,12],'degree':8,'source_quartics':468,'target_patterns':112020,
 'matching_target_pairs':24,'remaining_source_indices':ids,'all_quotient_orders':2,
 'global_cusp_order':5971968,'integral_unit_lattice_rank':23,'unit_relation_index':5971968,
 'unique_effective_global_L8_classes':6,'char0_pencil_dimensions':[2]*6,'pole_degrees':[8]*6,
 'ramification_degrees':[[8]*4]*6,'ramification_multiplicities':[[1]*4]*6,
 'remaining':0,'boundary':'Independent integer relation arithmetic, complete diagonal assembly, global principal-lattice reductions, six-class coverage and native exact characteristic-zero outputs. Native class-group and ramification-factorization algorithms are not reimplemented in Python.'}
(P/'data'/'replay_M6N12_diagonal_lattice_root.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2));print('M6N12_DIAGONAL_LATTICE_RAMIFICATION_ASSEMBLY_DONE')
