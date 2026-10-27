#!/usr/bin/env python3
"""Exact Tate primitive-column cusp audit for Y57=X1(57)/<7>."""
from collections import Counter
from math import gcd
from pathlib import Path
import json
P=Path(__file__).resolve().parent/'data';n=57
labels={(a,b) for a in range(n) for b in range(n) if gcd(gcd(a,b),n)==1}
def orbit(x):
 a,b=x
 return frozenset(((s*(a+k*b))%n,s*b%n) for s in (-1,1) for k in range(n))
cls=sorted(set(map(orbit,labels)),key=lambda x:sorted(x));ix={x:i for i,c in enumerate(cls) for x in c};reps=[min(c) for c in cls]
def cycles(p):
 left=set(range(len(p)));ans=[]
 while left:
  i=min(left);o=[]
  while i not in o:o.append(i);left.remove(i);i=p[i]
  assert i==o[0];ans.append(o)
 return ans
units=[u for u in range(n) if gcd(u,n)==1]
gal={u:[ix[(u*a%n,b)] for a,b in reps] for u in units}
def diamond(u):return [ix[(pow(u,-1,n)*a%n,u*b%n)] for a,b in reps]
ho=cycles(diamond(7));assert {len(o) for o in ho}=={3};qi={i:j for j,o in enumerate(ho) for i in o}
qgal={u:[qi[gal[u][o[0]]] for o in ho] for u in units}
qd5=[qi[diamond(5)[o[0]]] for o in ho];assert {len(o) for o in cycles(qd5)}=={6}
qf=cycles(qgal[2]);sf=cycles(gal[2]);admissible=set()
for o in sf:
 qo={qi[i] for i in o};assert any(qo==set(z) for z in qf)
 if len(o)<=8:
  assert len(o)==len(qo);admissible.update(qo)
left=set(range(len(ho)));go=[]
while left:
 i=min(left);o={qgal[u][i] for u in units};left-=o;go.append(o)
gd={i:len(o) for o in go for i in o};rows=[];frobs={}
for o in qf:
 width=n//gcd(reps[ho[o[0]][0]][1],n);assert all(n//gcd(reps[i][1],n)==width for j in o for i in ho[j])
 q=o[0];at=q
 for k in range(6):
  if at==qgal[2][q]:break
  at=qd5[at]
 else:raise AssertionError('No cusp Frobenius')
 if width in frobs:assert frobs[width]==k
 frobs[width]=k
 adm=set(o)<=admissible;assert not(set(o)&admissible) or adm
 if adm:assert gd[o[0]]==len(o)
 rows.append({'degree':len(o),'width':width,'global_degree':gd[o[0]],'admissible_from_source_D8':adm})
assert Counter(map(len,sf))=={1:18,2:9,9:2,18:1}
assert Counter(map(len,qf))=={1:6,2:3,3:2,6:1}
assert Counter(r['degree'] for r in rows if r['admissible_from_source_D8'])=={1:6,2:3}
a=[1]+[0]*8
for r in rows:
 if r['admissible_from_source_D8']:
  for i in range(r['degree'],9):a[i]+=a[i-r['degree']]
assert a[8]==3654
out={'level':57,'H_generator':7,'H_order':3,'H_free_on_cusps':True,
     'source_F2_signature':dict(Counter(map(len,sf))),
     'quotient_F2_signature':dict(Counter(map(len,qf))),
     'all_admissible_closed_cusp_divisors_lift_Q':True,
     'generator_diamond':5,'Frobenius_exponents':frobs,'cusps':rows,
     'admissible_divisors_D8':3654,'single_cusp_RR_kernels':60}
(P/'cusp57_quotient.json').write_text(json.dumps(out,indent=2)+'\n')
print('C57_CUSP_QUOTIENT_AUDIT_DONE',out['source_F2_signature'],out['quotient_F2_signature'],frobs,3654)
