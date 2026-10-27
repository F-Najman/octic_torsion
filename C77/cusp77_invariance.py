#!/usr/bin/env python3
"""Exact primitive-column proof of the special77cuspidal lift property."""
from collections import Counter
from math import gcd
from pathlib import Path
import json
H=Path(__file__).resolve().parent/'data';n=77
labels={(a,b) for a in range(n) for b in range(n) if gcd(gcd(a,b),n)==1}
def orbit(x):
    a,b=x
    return frozenset(((s*(a+k*b))%n,s*b%n) for s in (-1,1) for k in range(n))
classes=sorted(set(map(orbit,labels)),key=lambda x:sorted(x))
index={x:i for i,c in enumerate(classes) for x in c};reps=[min(c) for c in classes]
units=[u for u in range(n) if gcd(u,n)==1]
gal={u:[index[(u*a%n,b)] for a,b in reps] for u in units}
diamond=[index[(pow(23,-1,n)*a%n,23*b%n)] for a,b in reps]
assert all(diamond[diamond[diamond[i]]]==i for i in range(len(classes)))
left=set(range(len(classes)));global_orbits=[]
while left:
    i=min(left);o={gal[u][i] for u in units};left-=o;global_orbits.append(sorted(o))
global_degree={i:len(o) for o in global_orbits for i in o}
left=set(range(len(classes)));local_orbits=[]
while left:
    i=min(left);o=[]
    while i not in o:
        o.append(i);left.remove(i);i=gal[2][i]
    assert i==o[0];local_orbits.append(o)
needed=[o for o in local_orbits if len(o)<=8]
assert Counter(map(len,needed))=={1:30,3:10}
rows=[]
for o in needed:
    gd=global_degree[o[0]]
    if len(o)==1:
        assert gd==1
    else:
        assert gd==6
        assert {diamond[i] for i in o}==set(o)
        # This equality holds on the actual characteristic-zero cusp labels.
        # Consequently the three-cusp lift is diamond-invariant as a divisor.
        assert all(diamond[i]==gal[4][i] for i in o)
    rows.append({'F2_degree':len(o),'global_Galois_degree':gd,
                 'label_indices':o,'representatives':[reps[i] for i in o],
                 'diamond23_preserves_divisor':{diamond[i] for i in o}==set(o)})
out={'level':77,'degree':8,'total_geometric_cusps':len(classes),
     'global_signature':dict(Counter(map(len,global_orbits))),
     'F2_signature':dict(Counter(map(len,local_orbits))),
     'cusps_through_degree8':rows,
     'nonQ_lifted_cuspidal_divisors_are_diamond23_invariant':True,
     'rational_remainder_cusps_lift_toQ':True}
(H/'cusp77_invariance.json').write_text(json.dumps(out,indent=2)+'\n')
print('OCTIC77_SPECIAL_CUSP_LIFTS_DONE',out['global_signature'],out['F2_signature'])
