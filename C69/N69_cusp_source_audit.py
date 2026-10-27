#!/usr/bin/env python3
"""Cusp-complement audit for the octic69 quotient computation."""
from math import gcd
from pathlib import Path
from sympy import divisors, totient, kronecker_symbol
import json
N=69
units=[a for a in range(1,N) if gcd(a,N)==1]
rows=[]
for d in divisors(N):
    m=N//d
    stab={a for a in units if a%m in ({1%m,(-1)%m} if d<=2 else {1%m})}
    size=len(units)//len(stab)
    number=int(totient(d)*totient(m))//2
    assert number%size==0
    f=1
    while pow(2,f,N) not in stab:f+=1
    rows.append({'cusp_type':int(d),'Q_orbit_degree':size,
                 'Q_orbits':number//size,'F2_degree':f})
assert rows==[
 {'cusp_type':1,'Q_orbit_degree':22,'Q_orbits':1,'F2_degree':22},
 {'cusp_type':3,'Q_orbit_degree':22,'Q_orbits':1,'F2_degree':11},
 {'cusp_type':23,'Q_orbit_degree':2,'Q_orbits':11,'F2_degree':2},
 {'cusp_type':69,'Q_orbit_degree':1,'Q_orbits':22,'F2_degree':1}]
assert all(r['Q_orbit_degree']==1 for r in rows if r['F2_degree']==1)
H={a for a in units if kronecker_symbol(69,a)==1}
assert len(H)==22 and 5 in H and 68 in H
out={'degree':8,'level':69,'source_cusp_orbits':rows,
     'all_F2_rational_source_cusps_Q_rational':True,
     'quadratic_quotient_H':sorted(H),'diamond5_trivial_on_quotient':True,
     'consequence':'For a source degree7 noncuspidal part, its degree1 cuspidal complement maps to a rational cusp killed as a divisor by T5-6 on the quotient.'}
(Path(__file__).resolve().parent/'data'/'N69_cusp_source_audit.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
print('OCTIC69_SOURCE_CUSP_COMPLEMENT_AUDIT_DONE')
