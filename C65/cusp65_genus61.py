#!/usr/bin/env python3
"""Exact Tate cusp action for X1(65)/<51>, independently of curve models."""
from collections import Counter
from math import gcd
import json
from pathlib import Path

N=65
def canonical(c,a):
    c %= N
    d=gcd(c,N)
    a %= d
    assert gcd(a,d)==1
    return min((c,a),((-c)%N,(-a)%d))

labels=sorted({canonical(c,a) for c in range(N) for a in range(gcd(c,N)) if gcd(a,gcd(c,N))==1})
assert len(labels)==96
def diamond(P,u):
    return canonical(u*P[0],u*P[1])
def frobenius(P):
    return canonical(P[0],2*P[1])
def cycles(S, action):
    remaining=set(S)
    out=[]
    while remaining:
        p=min(remaining)
        orb=[p]
        q=action(p)
        while q!=p:
            assert q not in orb
            orb.append(q)
            q=action(q)
        assert set(orb)<=remaining
        remaining-=set(orb)
        out.append(orb)
    return out

source_closed=cycles(labels,frobenius)
rational_source=[o[0] for o in source_closed if len(o)==1]
units=[u for u in range(N) if gcd(u,N)==1]
assert len(rational_source)==24 and all(canonical(c,u*a)==(c,a) for c,a in rational_source for u in units)
H_orbits=cycles(labels,lambda p:diamond(p,51))
assert len(H_orbits)==48 and all(len(o)==2 for o in H_orbits)
quotient={p:min(o) for o in H_orbits for p in o}
target_labels=sorted(set(quotient.values()))
target_closed=cycles(target_labels,lambda p:quotient[frobenius(p)])
place_of={p:i for i,o in enumerate(target_closed) for p in o}
source_place_of={p:i for i,o in enumerate(source_closed) for p in o}
permutations=[]
for u in [27,41]:
    perm=[place_of[quotient[diamond(o[0],u)]] for o in target_closed]
    assert sorted(perm)==list(range(len(target_closed)))
    permutations.append(perm)

places=[]
for i,o in enumerate(target_closed):
    above=sorted({source_place_of[p] for p in labels if place_of[quotient[p]]==i})
    source_degrees=[len(source_closed[k]) for k in above]
    degree=len(o)
    assert sum(source_degrees)==2*degree
    assert all(d%degree==0 for d in source_degrees)
    widths={N//gcd(p[0],N) for p in o}
    assert len(widths)==1
    allowed=[d//degree for d in source_degrees if d<=8]
    places.append({'index':i,'degree':degree,'width':next(iter(widths)),
                   'source_closed_degrees':source_degrees,
                   'allowed_positive_push_multiplicities':sorted(set(allowed))})

admissible=[p for p in places if p['allowed_positive_push_multiplicities']]
assert Counter(len(o) for o in source_closed)=={1:24,4:6,6:4,12:2}
assert Counter((p['degree'],p['width']) for p in places)=={(1,65):12,(2,13):6,(6,1):2,(6,5):2}
assert all(p['allowed_positive_push_multiplicities']==([2] if p['degree']==2 else [1]) for p in admissible)
assert Counter(p['degree'] for p in admissible)=={1:12,2:6,6:2}
dp=[1]+[0]*8
for p in admissible:
    weight=p['degree']*p['allowed_positive_push_multiplicities'][0]
    for d in range(weight,9):
        dp[d]+=dp[d-weight]
assert dp[8]==83949

remaining=set(p['index'] for p in admissible)
deck_orbits=[]
while remaining:
    orb={min(remaining)}
    while True:
        more=orb|{perm[i] for i in orb for perm in permutations}
        if more==orb:
            break
        orb=more
    assert orb<=remaining
    remaining-=orb
    deck_orbits.append(sorted(orb))
assert sorted(map(len,deck_orbits))==[2,6,12]
blocks=sum(8//(places[o[0]]['degree']*places[o[0]]['allowed_positive_push_multiplicities'][0]) for o in deck_orbits)
assert blocks==11
result={'level':65,'degree':8,'quotient_diamond':51,'source_geometric_cusps':96,'F2_rational_source_cusps_are_Q_rational':True,
        'target_geometric_cusps':48,'source_signature':dict(sorted(Counter(map(len,source_closed)).items())),
        'target_places':places,'diamond27_permutation':permutations[0],
        'diamond41_permutation':permutations[1],'admissible_deck_orbits':deck_orbits,
        'admissible_degree8_divisors':dp[8],'source_RR_blocks':blocks}
(Path(__file__).resolve().parent/'data'/'cusp65_genus61.json').write_text(json.dumps(result,indent=2)+'\n')
print({k:v for k,v in result.items() if k not in ['target_places','diamond27_permutation','diamond41_permutation','admissible_deck_orbits']})
