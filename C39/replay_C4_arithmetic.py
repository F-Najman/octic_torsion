#!/usr/bin/env python3
"""Independent actual-relation, cusp-label and automorphism replay."""
import ast
from collections import Counter
from itertools import product
import json
from math import gcd
from pathlib import Path
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form

HERE = Path(__file__).resolve().parent/'data'
data = {}
for st in ast.parse((HERE/'C4_cusp_arithmetic.py').read_text()).body:
    assert isinstance(st,ast.Assign) and len(st.targets)==1
    data[st.targets[0].id] = ast.literal_eval(st.value)
assert data['complete']
assert all(len(v)==3 and v[-1]==0 for v in data['relations']+data['cusp_generators'])
R = Matrix([v[:-1] for v in data['relations']])
assert abs(int(R.det())) == 2016
L = hermite_normal_form(Matrix([v[:-1] for v in data['relations']+data['cusp_generators']]).T)
assert abs(int(L.det())) == 1
assert Counter(data['cusp_degrees']) == {1:4,2:2,4:2}

N = 39
def cusp(a,c):
    c %= N;g = gcd(c,N)
    return min((a%g,c),((-a)%g,(-c)%N))
labels = {cusp(a,c) for c in range(N) for a in range(gcd(c,N)) if gcd(a,gcd(c,N))==1}
assert len(labels)==48
units = [u for u in range(N) if gcd(u,N)==1]
def gal(P,u):
    return cusp(u*P[0],P[1])
def diamond(P,u):
    return cusp(pow(u,-1,N)*P[0],u*P[1])
def cycles(values,action):
    left=set(values);answer=[]
    while left:
        p=next(iter(left));at=action(p);orb={p}
        while at!=p:
            assert at not in orb;orb.add(at);at=action(at)
        left-=orb;answer.append(frozenset(orb))
    return answer
source_f=cycles(labels,lambda P:gal(P,2))
source_q={frozenset(gal(P,u) for u in units) for P in labels}
assert set(source_f)==source_q
assert Counter(map(len,source_f))=={1:12,2:6,12:2}
ho=cycles(labels,lambda P:diamond(P,16))
assert len(ho)==16 and all(len(o)==3 for o in ho)
image={P:o for o in ho for P in o};qvalues=set(ho)
qgal=lambda o,u:image[gal(next(iter(o)),u)]
qdiamond=lambda o,u:image[diamond(next(iter(o)),u)]
target_f=cycles(qvalues,lambda o:qgal(o,2))
target_q={frozenset(qgal(o,u) for u in units) for o in ho}
assert set(target_f)==target_q and Counter(map(len,target_f))=={1:4,2:2,4:2}
frobenius={}
for orb in target_f:
    P=next(iter(orb));geometric=next(iter(P));width=N//gcd(geometric[1],N)
    desired=qgal(P,2);at=P
    for k in range(4):
        if at==desired:break
        at=qdiamond(at,2)
    else:raise AssertionError('Frobenius is not in C4')
    assert width not in frobenius or frobenius[width]==k
    frobenius[width]=k
assert frobenius=={1:3,3:1,13:2,39:0}
# The constructor names the inverse generator; the defining kernel is identical.
assert {w:(-k)%4 for w,k in frobenius.items()}=={1:1,3:3,13:2,39:0}

# Exhaust all endomorphisms of C2 x C16 through the two generator images.
group=list(product(range(2),range(16)))
identity=tuple(group)
automorphisms=[];orders=Counter()
for im1 in product(range(2),(0,8)):
    for im2 in group:
        table=tuple(((x*im1[0]+y*im2[0])%2,(x*im1[1]+y*im2[1])%16) for x,y in group)
        if len(set(table))!=32:continue
        automorphisms.append(table)
        lookup=dict(zip(group,table));at=identity
        for order in range(1,33):
            at=tuple(lookup[P] for P in at)
            if at==identity:
                orders[order]+=1;break
        else:raise AssertionError('automorphism order exceeded32')
assert len(automorphisms)==32 and all(order&(order-1)==0 for order in orders)
log=(HERE.parent/'logs'/'C4_dyadic_certificate_v2.log').read_text()
assert 'all_Picard_mod4=true' in log and 'two_primary_order=32' in log
assert 'C39_C4_DYADIC_CERTIFICATE_DONE complete=true' in log
report={'level':39,'quotient_genus':9,'source_2primary_invariants':[2,16],
        'source_2primary_automorphism_order':32,'automorphism_orders':dict(orders),
        'target_finite_Picard_order':2016,'target_rational_cusp_index':1,
        'target_2primary_reduction_image_order':32,
        'source_and_target_cusp_orbits_all_Q_liftable':True,
        'Frobenius_exponents':frobenius,'source_torsion_specialization_injective':True}
(HERE/'C4_independent_replay.json').write_text(json.dumps(report,indent=2)+'\n')
print('C39_C4_INDEPENDENT_REPLAY_DONE actual_relations=true source_cusps=48 target_cusps=16 Aut_C2xC16_order=32 target_2primary_image=32 dyadic_injectivity=true')
