#!/usr/bin/env python3
"""Independent primitive-vector cusp enumeration, without orbit formulas."""
from collections import Counter
from math import gcd
from pathlib import Path
import json

N = 69
def cusp(a,c):
    return min(((s*(a+k*c)) % N,(s*c) % N)
               for s in (1,-1) for k in range(N))
classes = {cusp(a,c) for a in range(N) for c in range(N) if gcd(gcd(a,c),N) == 1}
assert len(classes) == 88
units = [u for u in range(N) if gcd(u,N) == 1]
def act(u,v): return cusp(u*v[0],v[1])
qorbits = set()
for v in classes: qorbits.add(frozenset(act(u,v) for u in units))
forbit = set()
for v in classes:
    orbit = set(); w = v
    while w not in orbit: orbit.add(w); w = act(2,w)
    forbit.add(frozenset(orbit))
assert Counter(map(len,qorbits)) == {1:22,2:11,22:2}
assert Counter(map(len,forbit)) == {1:22,2:11,11:2,22:1}
assert all(all(act(u,v) == v for u in units) for v in classes if act(2,v) == v)
def legendre(a,p):
    b = pow(a % p,(p-1)//2,p)
    return -1 if b == p-1 else b
H = [u for u in units if legendre(u,3)*legendre(u,23) == 1]
assert len(H) == 22 and 5 in H and 68 in H
report = {'level':69,'geometric_cusps':88,
          'Q_orbit_degrees':dict(Counter(map(len,qorbits))),
          'F2_orbit_degrees':dict(Counter(map(len,forbit))),
          'every_F2_rational_source_cusp_Q_rational':True,
          'quadratic_kernel':H,'diamond5_trivial':True}
(Path(__file__).resolve().parent/'data'/'audit69_cusps.json').write_text(json.dumps(report,indent=2)+'\n')
print('C69_INDEPENDENT_PRIMITIVE_VECTOR_CUSP_AUDIT_DONE total88=true rational_source_cusps22=true quadratic_kernel22=true diamond5=true')
