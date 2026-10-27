#!/usr/bin/env python3
"""Exact cusp freeness and Hecke complement conditions for the31-quotient."""
from collections import Counter
from math import gcd
from pathlib import Path
import json
n=42
labels={(a,b,c,d) for a in range(n) for b in range(n) if gcd(gcd(a,b),n)==1 for c in range(2) for d in range(2) if (c*b-d*a)%2==1}
def orbit(x):
 a,b,c,d=x
 return frozenset(((s*(a+k*b))%n,s*b%n,s*(c+k*d)%2,s*d%2) for s in [-1,1] for k in range(n))
cs=sorted(set(map(orbit,labels)),key=lambda s:sorted(s));lookup={x:i for i,c in enumerate(cs) for x in c};assert len(cs)==72
pp=[];fp=[]
for c in cs:
 a={lookup[(29*a%n,29*b%n,c,d)] for a,b,c,d in c};assert len(a)==1;pp.append(a.pop())
 a={lookup[(5*a%n,b,c,d)] for a,b,c,d in c};assert len(a)==1;fp.append(a.pop())
assert all(pp[i]!=i and pp[pp[i]]==i for i in range(72))
pairs=sorted({tuple(sorted([i,pp[i]])) for i in range(72)});pi={j:i for i,c in enumerate(pairs) for j in c}
qf=[]
for c in pairs:
 a={pi[fp[i]] for i in c};assert len(a)==1;qf.append(a.pop())
def cycles(perm):
 left=set(range(len(perm)));out=[]
 while left:
  i=min(left);c=[]
  while i not in c:c.append(i);left.remove(i);i=perm[i]
  assert i==c[0];out.append(c)
 return out
source=cycles(fp);target=cycles(qf)
assert Counter(map(len,source))=={1:18,2:9,3:6,6:3}
assert Counter(map(len,target))=={1:18,3:6}
allowed={}
for degree in [1,2,3,4]:
 aa=[]
 for q in range(1,n):
  if gcd(q,n)>1:continue
  perm=[]
  for c in cs:
   a={lookup[(q*a%n,b,c,d)] for a,b,c,d in c};assert len(a)==1;perm.append(a.pop())
  if all(perm[c[0]] in c for c in source if len(c)<=degree):aa.append(q)
 allowed[degree]=aa
assert 11 in allowed[1] and 17 in allowed[4]
out={'source_genus':61,'quotient_genus':31,'source_geometric_cusps':72,'free_delta29_on_cusps':True,'source_f5_cusp_signature':dict(Counter(map(len,source))),'quotient_f5_cusp_signature':dict(Counter(map(len,target))),'allowed_Hecke_residues_by_complement_degree':allowed,'noncuspidal_freeness_reason':'delta29 fixes the14-point3Q, hence a7-point; a nontrivial elliptic automorphism cannot fix it','rational_base_cusp_label':[0,1,1,0]}
p=Path(__file__).resolve().parent;(p/'data'/'C242_forget3sign_Tate.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2));print('C242_FORGET3SIGN_TATE_DONE')
