#!/usr/bin/env python3
"""Independent HNF replay of exact cusp containment and the 12-class orbit."""
import ast, hashlib, json
from math import gcd, prod
from pathlib import Path
from sympy import Matrix, ilcm
from sympy.matrices.normalforms import hermite_normal_form
P=Path(__file__).resolve().parent/'data'
shas={}
def read(name):
 t=(P/name).read_text();shas[name]=hashlib.sha256(t.encode()).hexdigest()
 return ast.literal_eval(t.replace('[*','[').replace('*]',']'))
def basis(rows,n):
 H=hermite_normal_form(Matrix(rows).T);assert H.shape==(n,n) and H.det();return H
def keyfun(H):
 V=H.inv();den=int(ilcm(*[x.q for x in V]));A=[[int(den*x) for x in row] for row in V.tolist()]
 return lambda v:tuple(sum(a*b for a,b in zip(row,v))%den for row in A)
def finite(inv,R,cs):
 assert inv[-1]==0 and all(len(v)==len(inv) and v[-1]==0 for v in R+cs)
 n=len(inv)-1;r=[v[:-1] for v in R];h=basis(r,n);order=abs(int(h.det()));assert order==prod(inv[:-1])
 assert order//abs(int(basis(r+[v[:-1] for v in cs],n).det()))==1990656
 return n,r,order,keyfun(h),keyfun(basis(r+[v[:-1] for v in cs],n))
def linear(v,rows):return [sum(c*row[j] for c,row in zip(v,rows)) for j in range(len(rows[0]))]
a=read('diagonal_3_15_D4_arithmetic.txt');sigs,ds,I4,R4,I7,R7,C4,C7,IH,RH,E4,E7,S4,S7,pairs,counts=a
assert len(pairs)==60 and len(sigs)==20 and ds.count(1)==16 and ds.count(4)==4
n4,r4,o4,k4,modc4=finite(I4,R4,C4);n7,r7,o7,k7,modc7=finite(I7,R7,C7)
assert o4==5971968 and o7==149299200
B=basis([r+[0]*n7 for r in r4]+[[0]*n4+r for r in r7]+[x[:-1]+y[:-1] for x,y in zip(C4,C7)],n4+n7)
diag=keyfun(B);assert o4*o7//abs(int(B.det()))==1990656
exact=[]
for i,j in pairs:
 x,y=S4[i-1][2],S7[j-1][3]
 assert not any(modc4(x[:-1])) and not any(modc7(y[:-1]))
 if not any(diag(x[:-1]+y[:-1])):exact.append([i,j])
assert len(exact)==12
lift=read('M3N15_D4_exact_cusp_lifts.txt')
assert lift[:3]==[3,15,4] and lift[3:5]==[sigs,ds] and sorted(lift[8])==sorted(exact)
assert ds[lift[5]-1]==1 and len(lift[6])==len(lift[7])==12
for v,t in zip(lift[6],lift[7]):
 assert len(v)==20 and max(map(abs,v))<=1
 assert not any(k7([x-y for x,y in zip(linear(v,C7),t)][:-1]))
assert {k7(v[:-1]) for v in lift[7]}=={k7(S7[j-1][3][:-1]) for i,j in exact}
p13=read('diagonal_3_15_D4_p13_arithmetic.txt');p19=read('diagonal_3_15_D4_p19_arithmetic.txt')
assert gcd(gcd(o7,prod(p13[4][:-1])),prod(p19[4][:-1]))==3*1990656
orb=read('M3N15_D4_exact_orbit_arithmetic.txt');Io,Ro,Co,vs,ts,images,pair=orb
assert vs==lift[6] and len(images)==48 and pair[0]<=pair[1]
no,ro,oo,ko,modco=finite(Io,Ro,Co)
Bo=basis([r+[0]*no for r in r7]+[[0]*n7+r for r in ro]+[x[:-1]+y[:-1] for x,y in zip(C7,Co)],n7+no)
do=keyfun(Bo);assert o7*oo//abs(int(Bo.det()))==1990656
for x,y in zip(lift[7],ts):assert not any(do(x[:-1]+y[:-1]))
assert len({ko(t[:-1]) for t in ts})==12
assert {ko(t[:-1]) for t in ts}=={ko(x[:-1]) for x in images}
report={'degree':8,'group':[3,15],'cusp_group_order':1990656,'global_J_over_C_order_divides':3,
 'coarse_pairs':60,'all_coarse_classes_in_C_at_F4':True,'exact_unmultiplied_pairs':12,
 'global_signed_cusp_lifts':12,'max_abs_lift_coefficient':1,'automorphisms_tested':48,
 'single_global_class_orbit_size':12,'independent_actual_relation_HNF_replay':True,'input_sha256':shas}
(P/'replay_exact_3_15_D4.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2));print('M3N15_EXACT_CUSP_CONTAINMENT_ORBIT_REPLAY_DONE')
