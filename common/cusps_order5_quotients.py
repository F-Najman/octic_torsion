#!/usr/bin/env python3
"""Tate primitive-column audit for the order-five quotients at55and77."""
from collections import Counter
from math import gcd
from pathlib import Path
import json
H=Path(__file__).resolve().parent/'data'
def cycles(p):
 left=set(range(len(p)));ans=[]
 while left:
  i=min(left);o=[]
  while i not in o:o.append(i);left.remove(i);i=p[i]
  assert i==o[0];ans.append(o)
 return ans
def coeff(ws,d):
 a=[0]*(d+1);a[0]=1
 for w in ws:
  for i in range(w,d+1):a[i]+=a[i-w]
 return a[d]
def audit(n,h):
 labels={(a,b) for a in range(n) for b in range(n) if gcd(gcd(a,b),n)==1}
 def orbit(x):
  a,b=x
  return frozenset(((s*(a+k*b))%n,s*b%n) for s in (-1,1) for k in range(n))
 cls=sorted(set(map(orbit,labels)),key=lambda x:sorted(x))
 ix={x:i for i,c in enumerate(cls) for x in c};reps=[min(c) for c in cls]
 units=[u for u in range(n) if gcd(u,n)==1]
 gal={u:[ix[(u*a%n,b)] for a,b in reps] for u in units}
 def diamond(u):return [ix[(pow(u,-1,n)*a%n,u*b%n)] for a,b in reps]
 ho=cycles(diamond(h));assert {len(o) for o in ho}=={5}
 qi={i:j for j,o in enumerate(ho) for i in o}
 qgal={u:[qi[gal[u][o[0]]] for o in ho] for u in units}
 dp=diamond(2);qp=[qi[dp[o[0]]] for o in ho]
 cover_degree=(n-1 if False else len(cycles(qp)[0]))
 assert {len(o) for o in cycles(qp)}=={4 if n==55 else 6}
 sf=cycles(gal[2]);qf=cycles(qgal[2])
 admissible=set()
 for o in sf:
  qset={qi[i] for i in o}
  assert any(qset==set(z) for z in qf)
  if len(o)<=8:
   assert len(o)==len(qset)
   admissible.update(qset)
 left=set(range(len(ho)));go=[]
 while left:
  i=min(left);o={qgal[u][i] for u in units};left-=o;go.append(o)
 gd={i:len(o) for o in go for i in o}
 rows=[];frobs={}
 for o in qf:
  width=n//gcd(reps[ho[o[0]][0]][1],n)
  assert all(n//gcd(reps[i][1],n)==width for j in o for i in ho[j])
  q=o[0];at=q;k=0
  while at!=qgal[2][q] or (k==0 and qgal[2][q]!=q):at=qp[at];k+=1
  assert k<cover_degree
  if width in frobs:assert frobs[width]==k
  frobs[width]=k
  adm=set(o)<=admissible
  assert not(set(o)&admissible) or adm
  rows.append({'degree':len(o),'width':width,'global_degree':gd[o[0]],'source_degree8_admissible':adm})
 assert all(r['global_degree']==1 for r in rows if r['degree']==1)
 expect=({55:0,11:1,5:2,1:3} if n==55 else {77:0,11:2,7:3,1:5})
 assert frobs==expect
 full=Counter(r['degree'] for r in rows);ads=Counter(r['degree'] for r in rows if r['source_degree8_admissible'])
 assert full==({1:4,2:2,4:2} if n==55 else {1:6,2:3,3:2,6:1})
 assert ads==({1:4,4:1} if n==55 else {1:6,3:2})
 count=coeff([r['degree'] for r in rows if r['source_degree8_admissible']],8)
 assert count==(201 if n==55 else 1854)
 return {'level':n,'H_generator':h,'H_order':5,'H_free_on_cusps':True,
  'source_F2_signature':dict(Counter(map(len,sf))),
  'quotient_F2_signature':dict(full),'admissible_F2_signature':dict(ads),
  'Frobenius_exponents_relative_to_diamond2':frobs,
  'all_rational_F2_cusps_lift_toQ':True,'closed_quotient_cusps':rows,
  'admissible_degree8_divisors':count,'rational_only_degree8_divisors':coeff([1]*full[1],8)}
out=[audit(55,16),audit(77,36)]
(H/'cusps_order5_quotients.json').write_text(json.dumps(out,indent=2)+'\n')
for r in out:print('OCTIC_ORDER5_QUOTIENT_CUSP_AUDIT_DONE',r['level'],r['admissible_F2_signature'],r['admissible_degree8_divisors'],r['Frobenius_exponents_relative_to_diamond2'])
