#!/usr/bin/env python3
"""Independent exact Tate-label audit for X1(75)/<26>, in characteristic 2."""
from math import gcd
from collections import Counter, defaultdict
from pathlib import Path
import json
N=75
D=[d for d in range(1,N+1) if N%d==0]
def U(m): return [a for a in range(m) if gcd(a,m)==1]
def neg(c):
 d,a,b=c;return (d,(-a)%d,(-b)%(N//d))
def can(c): return min(c,neg(c))
S={can((d,a,b)) for d in D for a in U(d) for b in U(N//d)}
def diamond(c,u):
 d,a,b=c;return can((d,u*a%d,u*b%(N//d)))
def frob(c):
 d,a,b=c;return can((d,2*a%d,b))
def qcan(c): return min(c,diamond(c,26))
Q={qcan(c) for c in S}
assert len(S)==112 and len(Q)==56
assert all(diamond(c,26)!=c and diamond(diamond(c,26),26)==c for c in S)
def orbits(A,act):
 left=set(A);out=[]
 while left:
  x=min(left);o=[];y=x
  while y not in o:
   o.append(y);y=act(y)
  assert y==x and set(o)<=left
  left.difference_update(o);out.append(tuple(o))
 return out
SO=orbits(S,frob);QO=orbits(Q,lambda c:qcan(frob(c)))
qindex={c:i for i,o in enumerate(QO) for c in o}
steps=defaultdict(set)
for o in SO:
 i=qindex[qcan(o[0])]
 assert all(qindex[qcan(c)]==i for c in o)
 assert len(o)%len(QO[i])==0
 steps[i].add(len(o)//len(QO[i]))
rows=[]
for i,o in enumerate(QO):
 widths={N//c[0] for c in o};assert len(widths)==1
 assert len(steps[i])==1
 rows.append({'degree':len(o),'width':widths.pop(),'source_push_step':next(iter(steps[i])), 'representative':list(min(o))})
expected={(1,25,2):10,(1,75,1):10,(4,5,1):2,(4,15,1):2,(10,1,2):1,(10,3,2):1}
actual=Counter((r['degree'],r['width'],r['source_push_step']) for r in rows)
assert actual==expected,(actual,expected)
# Degree-eight pushforwards of source entirely cuspidal divisors.
coef=[1]+[0]*8
for r in rows:
 w=r['degree']*r['source_push_step']
 for j in range(w,9):coef[j]+=coef[j-w]
assert coef[8]==131790
# Width75 source cusps are rational: every cyclotomic action fixes them.
assert all(can((c[0],u*c[1]%c[0],c[2]))==c for c in S if N//c[0]==75 for u in U(N))
# Diamond2 on closed target places (it commutes with geometric Frobenius).
perm={i:qindex[qcan(diamond(o[0],2))] for i,o in enumerate(QO)}
for i,o in enumerate(QO):
 assert all(qindex[qcan(diamond(c,2))]==perm[i] for c in o)
small={i for i,o in enumerate(QO) if len(o)<=8}
cycles=orbits(small,lambda i:perm[i])
cycle_summary=Counter((len(QO[o[0]]),N//QO[o[0]][0][0],len(o)) for o in cycles)
assert cycle_summary==Counter({(1,25,10):1,(1,75,10):1,(4,5,2):1,(4,15,2):1})
assert len([o for o in SO if len(o)==1])==20
assert {N//o[0][0] for o in SO if len(o)==1}=={75}
result={'N':N,'H_generator':26,'H_order_on_projective_labels':2,'source_geometric_cusps':len(S),'target_geometric_cusps':len(Q),'source_F2_signature':dict(Counter(map(len,SO))),'target_F2_signature':dict(Counter(map(len,QO))),'admissible_degree8_divisors':coef[8],'by_degree_width_step':[[*k,v] for k,v in sorted(actual.items())],'target_closed_cusps':rows,'diamond2_closed_cusp_permutation':perm,'diamond2_low_cusp_cycles':[list(o) for o in cycles],'diamond2_cycles_by_degree_width_length':[[*k,v] for k,v in sorted(cycle_summary.items())],'source_F2_rational_cusps_all_Q_rational':True}
path=Path(__file__).resolve().parent/'data'/'audit75_tate_parities.json';path.write_text(json.dumps(result,indent=2)+'\n')
print('PASS C75 free involution; source112,target56; signature1^20,4^4,10^2; 10width25 even,10width75 arbitrary,all4quartics arbitrary; admissibleD8=131790; diamond2cycles10,10,2,2; all20sourceF2rationalcuspsareQrational')
