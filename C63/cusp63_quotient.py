#!/usr/bin/env python3
"""Exact primitive-column cusp descent for Y=X1(63)/<8>."""
from collections import Counter
from math import gcd
from pathlib import Path
import json
H=Path(__file__).resolve().parent/'data';n=63
labels={(a,b) for a in range(n) for b in range(n) if gcd(gcd(a,b),n)==1}
def orbit(x):
 a,b=x
 return frozenset(((s*(a+k*b))%n,s*b%n) for s in (-1,1) for k in range(n))
classes=sorted(set(map(orbit,labels)),key=lambda x:sorted(x))
index={x:i for i,c in enumerate(classes) for x in c};reps=[min(c) for c in classes]
units=[u for u in range(n) if gcd(u,n)==1]
gal={u:[index[(u*a%n,b)] for a,b in reps] for u in units}
diamond=[index[(8*a%n,8*b%n)] for a,b in reps]
assert all(diamond[i]!=i and diamond[diamond[i]]==i for i in range(len(classes)))
qo=sorted({tuple(sorted((i,diamond[i]))) for i in range(len(classes))})
qi={i:j for j,o in enumerate(qo) for i in o}
qgal={u:[qi[gal[u][o[0]]] for o in qo] for u in units}
assert all(qi[gal[u][o[0]]]==qi[gal[u][o[1]]] for u in units for o in qo)
left=set(range(len(qo)));go=[]
while left:
 i=min(left);o={qgal[u][i] for u in units};left-=o;go.append(o)
gd={i:len(o) for o in go for i in o}
left=set(range(len(qo)));fo=[]
while left:
 i=min(left);o=[]
 while i not in o:
  o.append(i);left.remove(i);i=qgal[2][i]
 assert i==o[0];fo.append(o)
assert Counter(map(len,fo))=={1:9,2:3,3:9,6:1}
assert all(gd[o[0]]==1 for o in fo if len(o)==1)
out={'level':63,'quotient_diamond':8,'source_geometric_cusps':len(classes),
 'involution_free_on_all_cusps':True,'quotient_geometric_cusps':len(qo),
 'global_signature':dict(Counter(map(len,go))),
 'F2_signature':dict(Counter(map(len,fo))),
 'all9_F2_rational_cusps_lift_toQ':True}
(H/'cusp63_quotient.json').write_text(json.dumps(out,indent=2)+'\n')
print('OCTIC63_QUOTIENT_CUSPS_DONE',out)
