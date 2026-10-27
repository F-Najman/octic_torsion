#!/usr/bin/env python3
"""Exact diagnostic diagonal after complete Hecke and stabilizer pruning."""
from pathlib import Path
import ast,json
from math import prod
from collections import Counter,defaultdict
from sympy import Matrix,ilcm
from sympy.matrices.normalforms import hermite_normal_form
P=Path(__file__).resolve().parent
def read(n):return ast.literal_eval((P/'data'/n).read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
def basis(rs,n):
 H=hermite_normal_form(Matrix(rs).T);assert H.shape==(n,n) and H.det();return H
def keyfun(H):
 V=H.inv();den=int(ilcm(*[v.q for v in V]));A=[[int(den*x) for x in row] for row in V.tolist()]
 return lambda v:tuple(sum(x*y for x,y in zip(row,v))%den for row in A)
g7=read('M3N18_9multiple_fresh_p7_group_data.m');g13=read('M3N18_9multiple_aut27_p13_group_data.m')
assert g7[1:3]==g13[1:3] and len(g7[1])==24
_,sigs,ds,i7,r7,c7,h7,rc7,f7=g7;_,_,_,i13,r13,c13,h13,rc13,f13=g13
assert i7[-1]==i13[-1]==0 and all(v[-1]==0 for v in r7+c7+rc7+r13+c13+rc13)
n7,n13=len(i7)-1,len(i13)-1
R7=basis([v[:-1] for v in r7],n7);R13=basis([v[:-1] for v in r13],n13)
o7,o13=abs(int(R7.det())),abs(int(R13.det()));assert o7==1640558367 and o13==177180303636
HC7=basis([v[:-1] for v in r7+c7],n7);HC13=basis([v[:-1] for v in r13+c13],n13)
assert o7//abs(int(HC7.det()))==o13//abs(int(HC13.det()))==20253807
DG=basis([v[:-1]+[0]*n13 for v in r7]+[[0]*n7+v[:-1] for v in r13]+[v[:-1]+w[:-1] for v,w in zip(c7,c13)],n7+n13)
assert o7*o13//abs(int(DG.det()))==20253807
key=keyfun(DG);k7=lambda v:key(v+[0]*n13);k13=lambda v:key([0]*n7+v)
# Independent Frobenius canonical key, for an explicitly presented prime field extension.
a7=read('M3N18_9multiple_fresh_p7_degree4_data.m');_,_,fp,classes,coords=a7;d=4;prime=7
assert fp[-1]==1 and len(fp)==5
def mul(a,b):
 c=[0]*7
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]=(c[i+j]+x*y)%7
 for k in range(6,3,-1):
  for j in range(4):c[k-4+j]=(c[k-4+j]-c[k]*fp[j])%7
 return c[:4]
def frob(a):
 z=[1,0,0,0];b=a
 for _ in range(7):z=mul(z,b)
 return z
def coordkey(v):
 v=[list(w) for w in v];out=[]
 for _ in range(4):out.append(tuple(v[0]+v[1]));v=[frob(w) for w in v]
 return min(out)
keys=[coordkey(v) for v in coords];assert len(set(keys))==513
ex,tests,rem=read('sieve_3_18_D4_fullaut_root_v4_arithmetic.txt')
assert [row[:2] for row in ex]==[[2,27],[3,144],[4,513]]
assert len([t for t in tests if t[0]==22])==27 and all(t[-1]==0 for t in tests if t[0]==22)
assert not any(t[-1] for t in tests if t[0]==2)
rs3=ex[1];cubic=[r[1] for r in rem if r[0]==3];assert cubic==[4,5,7,8]
assert all(Counter(rs3[3])[i]==9 for i in cubic)
row=ex[2];assert set(map(tuple,row[4]))==set(keys)
allowed={r[1] for r in rem if r[0]==4};assert allowed=={1,5,6,7,11,12,17,19}
selected={tuple(k) for k,l in zip(row[4],row[3]) if l in allowed};assert len(selected)==216
source=[(j+1,cl) for j,(cl,k) in enumerate(zip(classes,keys)) if k in selected]
S=defaultdict(list)
for j,cl in source:S[k7([9*x for x in cl[:-1]])].append(j)
# All lower-degree p13 places are included; quartics use a complete Aut27 partition.
a2=read('M3N18_9multiple_aut27_p13_degree2_data.m');a3=read('M3N18_9multiple_aut27_p13_degree3_data.m');a4=read('M3N18_9multiple_aut27_p13_degree4_data.m')
assert len(a2[3])==54 and len(a3[3])==846 and len(a4[6])==7506 and len(a4[3])==278
assert set(a4[7])==set(range(1,279)) and len(a4[8])==278
for j,ims in enumerate(a4[8],1):assert len(ims)==27 and set(map(tuple,ims))=={tuple(k) for k,l in zip(a4[6],a4[7]) if l==j}
C2,C3,C4=[a[3] for a in [a2,a3,a4]]
assert all(v[-1]==0 for v in C2+C3+C4)
rc=[v[:-1] for v in rc13];target=[];matches=[]
def test(typ,i,j,cl):
 idx=len(target)+1;target.append([typ,i,j,cl]);s=S.get(k13([-9*x for x in cl]));
 if s:matches.append([idx,s])
for i,cl in enumerate(C2,1):
 for j in range(18):
  for k in range(j,18):test(2,i,[j+1,k+1],[x+y+z for x,y,z in zip(cl[:-1],rc[j],rc[k])])
for i,cl in enumerate(C3,1):
 for j,r in enumerate(rc,1):test(3,i,j,[x+y for x,y in zip(cl[:-1],r)])
for i,cl in enumerate(C4,1):test(4,i,a4[5][i-1],cl[:-1])
for i,x in enumerate(C2,1):
 for j in range(i-1,54):test(22,i,j+1,[u+v for u,v in zip(x[:-1],C2[j][:-1])])
assert len(target)==26225
out={'degree':8,'group':[3,18],'complete_group_exclusion':False,'source7_surviving_quartics':216,'source7_9multiple_classes':len(S),'target13_cases':len(target),'target13_quartic_Aut_representatives':278,'matching_target_cases':len(matches),'matching_pairs':sum(len(s) for i,s in matches),'target_matching_patterns':dict(Counter(target[i-1][0] for i,s in matches)),'source_indices':sorted({j for i,s in matches for j in s})}
(P/'data'/'M3N18_9multiple_fresh_diagonal_matches_root.json').write_text(json.dumps(out,indent=2)+'\n');(P/'data'/'M3N18_9multiple_fresh_diagonal_match_data_root.json').write_text(json.dumps({'source':source,'target':target,'matches':matches},separators=(',',':'))+'\n')
print(json.dumps(out,indent=2));print('M3N18_9MULTIPLE_DIAGONAL_COMPARE_COMPLETE')
