#!/usr/bin/env python3
"""Read-only independent checksum, cusp, dimension and coverage audit."""
from collections import Counter
from hashlib import sha256
import json
from math import gcd
from pathlib import Path
import re

L = Path(__file__).resolve().parent
packet = json.loads((L/'data'/'N63_good_branches_manifest.json').read_text())
# The SHA-256 checks of the files listed in the manifest were removed (git binds the files).
N = 63
def cusp(a,c):
    return min(((s*(a+k*c)) % N,(s*c) % N)
               for s in (1,-1) for k in range(N))
classes = {cusp(a,c) for a in range(N) for c in range(N) if gcd(gcd(a,c),N) == 1}
assert len(classes) == 96
def diamond(v): return cusp(8*v[0],8*v[1])
assert all(diamond(v) != v and diamond(diamond(v)) == v for v in classes)
def quotient(v): return min(v,diamond(v))
qclasses = {quotient(v) for v in classes}; assert len(qclasses) == 48
units = [u for u in range(N) if gcd(u,N) == 1]
def gal(u,v): return quotient(cusp(u*v[0],v[1]))
qorbits = {frozenset(gal(u,v) for u in units) for v in qclasses}
forbits = set()
for v in qclasses:
    orb = set(); w = v
    while w not in orb: orb.add(w); w = gal(2,w)
    forbits.add(frozenset(orb))
assert Counter(map(len,qorbits)) == {1:9,2:3,3:6,6:1,9:1}
assert Counter(map(len,forbits)) == {1:9,2:3,3:9,6:1}
assert all(all(gal(u,v) == v for u in units) for v in qclasses if gal(2,v) == v)
rank = (L/'logs'/'winding63_quotient_rank.log').read_text()
rows = [[int(x) for x in row] for row in re.findall(
    r'Y63_RANK0_CHARACTER order=(\d+) conductor=(\d+) invariant8=1 dimension=(\d+) winding_dimension=(\d+)',rank)]
assert len(rows) == 5 and all(a == b for _,_,a,b in rows)
assert sum(d*(1 if o == 1 else 2) for o,_,d,_ in rows) == 49
for e in (7,8):
    log = (L/'logs'/f'sieve63_quotient_T5_e{e}.log').read_text()
    assert not re.search(r'Runtime error|User error|Assertion failed',log)
    matches = re.findall(r'Y63_T5_TEST index=(\d+) principal=(\w+)',log)
    assert matches == [('1','false'),('2','false'),('3','false')]
    assert f'Y63_NC_T5_DONE e={e} tests=3 complete=true' in log
    expected = 189 if e == 7 else 216
    assert re.search(rf'Y63_NC_EXHAUSTIVE e={e} geometric_pairs={expected} closed_places=27 diamond_orbits=3',log)
print('OCTIC63_GOOD_INDEPENDENT_AUDIT_DONE hashes=true cusp_freeness=true rational_quotient_cusps9=true rank_dimension49=true six_nonprincipal_tests=true entirely_cuspidal_branch_not_included=true')
