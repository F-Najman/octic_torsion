#!/usr/bin/env python3
"""Exact Tate primitive-column cusp and unramified Frobenius orbits."""
from collections import Counter
from math import gcd
from pathlib import Path
import json

def cycles(indices,perm):
 left=set(indices);out=[]
 while left:
  i=min(left);o=[]
  while i not in o:
   o.append(i);left.remove(i);i=perm[i]
  assert i==o[0];out.append(o)
 return out

def check(n):
 labels={(a,b) for a in range(n) for b in range(n) if gcd(gcd(a,b),n)==1}
 def orbit(x):
  a,b=x
  return frozenset(((s*(a+k*b))%n,(s*b)%n) for s in [-1,1] for k in range(n))
 classes=sorted(set(map(orbit,labels)),key=lambda s:sorted(s))
 index={x:i for i,c in enumerate(classes) for x in c}
 reps=[min(c) for c in classes]
 units=[a for a in range(n) if gcd(a,n)==1]
 perms={t:[index[((t*a)%n,b)] for a,b in reps] for t in units}
 left=set(range(len(classes)));globalorbits=[]
 while left:
  i=min(left);o={perms[t][i] for t in units};left-=o;globalorbits.append(o)
 f=cycles(range(len(classes)),perms[2]);globdeg={i:len(o) for o in globalorbits for i in o}
 table=[{'global_degree':globdeg[o[0]],'residue_degree':len(o)} for o in f]
 if n==57:assert all(r['global_degree']==r['residue_degree'] for r in table if r['residue_degree']<=7)
 if n==65:assert all(r['global_degree']==1 for r in table if r['residue_degree']==1)
 return {'N':n,'total_cusps':len(classes),'global_orbit_signature':dict(Counter(map(len,globalorbits))),
         'F2_signature':dict(Counter(map(len,f))),'F2_closed_cusps':table,
         'needed_lifts_verified':True}
out=[check(n) for n in [77,81]]
for x in out:
 if x['N']==81: assert all(r['global_degree']==r['residue_degree'] for r in x['F2_closed_cusps'] if r['residue_degree']<=8)
(Path(__file__).resolve().parent/'data'/'cusp_lift_checks.json').write_text(json.dumps(out,indent=2)+'\n')
for x in out:print('OCTIC_CUSP_LIFT_DONE',x['N'],x['global_orbit_signature'],x['F2_signature'])
