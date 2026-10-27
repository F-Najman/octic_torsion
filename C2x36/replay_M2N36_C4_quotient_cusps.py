#!/usr/bin/env python3
"""Exact geometric Tate cusp permutation and RH for the survivor C4 stabilizer."""
from collections import Counter
from math import gcd
from pathlib import Path
import json
m,n=2,36
labels={(a,b,c,d) for a in range(n) for b in range(n) if gcd(gcd(a,b),n)==1 for c in range(2) for d in range(2) if (c*b-d*a)%2==1}
def orbit(x):
 a,b,c,d=x
 return frozenset(((s*(a+k*b))%n,s*b%n,s*(c+k*d)%2,s*d%2) for s in [-1,1] for k in range(n))
classes=sorted(set(map(orbit,labels)),key=lambda s:sorted(s));index={x:i for i,cs in enumerate(classes) for x in cs}
def g(x):
 a,b,c,d=x
 return ((a+18*c)%36,(b+18*d)%36,(c+a)%2,(d+b)%2)
perm=[]
for cs in classes:
 im={index[g(x)] for x in cs};assert len(im)==1;perm.append(im.pop())
assert sorted(perm)==list(range(len(classes)))
def power(i,k):
 for _ in range(k):i=perm[i]
 return i
assert all(power(i,4)==i for i in range(len(classes)))
fix=[sum(power(i,k)==i for i in range(len(classes))) for k in [1,2,3]]
assert len(classes)==64 and fix==[0,0,0]
frob=[]
for cs in classes:
 im={index[(5*a%36,b,5*c%2,d)] for a,b,c,d in cs};assert len(im)==1;frob.append(im.pop())
left=set(range(len(classes)));frobdegrees=[]
while left:
 i=min(left);cy=[]
 while i not in cy:cy.append(i);left.remove(i);i=frob[i]
 assert i==cy[0];frobdegrees.append(len(cy))
assert Counter(frobdegrees)=={1:24,2:8,6:4}
left=set(range(len(classes)));orbits=[]
while left:
 i=min(left);o={power(i,k) for k in range(4)};left-=o;orbits.append(sorted(o))
ram=sum(fix);genus=(2*41-2-ram+8)//8
assert 2*41-2==4*(2*genus-2)+ram
assert genus==11
out={'source_genus':41,'geometric_cusps':len(classes),'closed_cusp_degrees_at_5':dict(Counter(frobdegrees)),'g':'(P2,Q36) -> (P2+18Q36,Q36+P2)','fixed_geometric_cusps_g_g2_g3':fix,'cusp_orbits':dict(Counter(map(len,orbits))),'noncuspidal_fixed_points':0,'noncuspidal_freeness_reason':'every power fixes 4Q of order9; any nontrivial elliptic automorphism has degree(u-1)<=4, so cannot fix it','ramification_sum':ram,'quotient_genus':genus,'warning':'This is a quotient-genus calculation, not a quadratic-point classification.'}
p=Path(__file__).resolve().parent/'data';(p/'M2N36_C4_quotient_cusps.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2));print('M2N36_C4_QUOTIENT_CUSPS_DONE')
