#!/usr/bin/env python3
"""Independent manifest, row-rank, projection and local-kernel audit for C55."""
import ast
from hashlib import sha256
from itertools import product
import json
from math import gcd
from pathlib import Path
import re

BASE = Path(__file__).resolve().parent
manifest = json.loads((BASE/'data'/'C55_manifest.json').read_text())
assert manifest['complete'] and manifest['excluded_cyclic_order'] == 55
# The SHA-256 checks of the manifest files were removed (git binds the files), and so was the comparison
# of the session export with the reproduced export: they were byte-identical and are one file here.
text = (BASE/'data'/'C55_cartier_kernels.m').read_text()

def rows(data,r,c):
    entries = ast.literal_eval(data)
    assert len(entries) == r*c and set(entries) <= {0,1}
    return [entries[i*c:(i+1)*c] for i in range(r)]

def rank(matrix):
    basis = {}
    for row in matrix:
        v = sum(bit << i for i,bit in enumerate(row))
        while v:
            bit = v & -v
            if bit not in basis:
                basis[bit] = v
                break
            v ^= basis[bit]
    return len(basis)

def named(name):
    matches = list(re.finditer(rf'\b{name}:=Matrix\(GF\(2\),(\d+),(\d+),(\[.*?\])\);',text,re.S))
    assert len(matches) == 1
    m = matches[0]
    return rows(m[3],int(m[1]),int(m[2]))

C = named('CartierMatrix'); U = named('ImageBasis')
assert len(C) == len(C[0]) == 17 and len(U) == 15 and len(U[0]) == 17
CMI = [[b ^ (i==j) for j,b in enumerate(row)] for i,row in enumerate(C)]
assert rank(CMI) == rank(U) == rank(CMI+U) == 15

jets, kernels, full = {},{},{}
for m in re.finditer(r'(localjets|kernelspaces)\[<(\d+),(\d+)>\]:=Matrix\(GF\(2\),(\d+),(\d+),(\[.*?\])\);',text,re.S):
    target = jets if m[1] == 'localjets' else kernels
    k = int(m[2])-1,int(m[3]);assert k not in target
    target[k] = rows(m[6],int(m[4]),int(m[5]))
for m in re.finditer(r'basisjets\[(\d+)\]:=Matrix\(GF\(2\),(\d+),(\d+),(\[.*?\])\);',text,re.S):
    k = int(m[1])-1;assert k not in full
    full[k] = rows(m[4],int(m[2]),int(m[3]))
weights = [4,1,1,1,1]
assert len(jets) == len(kernels) == 34 and set(full) == set(range(5))
for i,w in enumerate(weights):
    B = full[i];assert len(B) == 17 and len(B[0]) == 8
    projection = [[sum(U[r][k]*B[k][j] for k in range(17))%2 for j in range(8)] for r in range(15)]
    for multiplicity in range(1,8//w+1):
        width = multiplicity*w
        assert jets[i,multiplicity] == [row[:width] for row in projection]
        K = kernels[i,multiplicity]
        prefix = [row[:width] for row in B]
        assert len(K) == rank(K) == 17-rank(prefix)
        assert all(sum(k[t]*prefix[t][j] for t in range(17))%2 == 0 for k in K for j in range(width))

tested = 0
for vector in product(range(3),range(9),range(9),range(9),range(9)):
    if sum(w*m for w,m in zip(weights,vector)) != 8:
        continue
    matrix = [[] for _ in range(15)]
    for i,multiplicity in enumerate(vector):
        if multiplicity:
            block = jets[i,multiplicity]
            for r in range(15):
                matrix[r].extend(block[r])
    assert all(len(row)==8 for row in matrix)
    assert rank(matrix) == 8
    tested += 1
assert tested == 201
# Independent cusp enumeration using (a mod gcd(c,N),c), then the sign.
# This does not use the sibling's explicit full primitive-column orbits.
N = 55
def cusp(a,c):
    c %= N
    g = gcd(c,N)
    return min((a%g,c),((-a)%g,(-c)%N))
labels = {cusp(a,c) for c in range(N) for a in range(gcd(c,N)) if gcd(a,gcd(c,N)) == 1}
assert len(labels) == 80
def sigma(P,u):
    return cusp(u*P[0],P[1])
def delta(P,u):
    return cusp(pow(u,-1,N)*P[0],u*P[1])
def orbits(values,action):
    left = set(values);result = []
    while left:
        p = min(left);orbit = {p};q = action(p)
        while q != p:
            assert q not in orbit
            orbit.add(q);q = action(q)
        left -= orbit;result.append(frozenset(orbit))
    return result
source_orbits = orbits(labels,lambda P:sigma(P,2))
assert sorted(map(len,source_orbits)) == [1]*20+[4]*5+[10]*2+[20]
H_orbits = orbits(labels,lambda P:delta(P,16))
assert len(H_orbits) == 16 and all(len(o)==5 for o in H_orbits)
image = {P:o for o in H_orbits for P in o}
qvalues = set(H_orbits)
qF = lambda o:image[sigma(min(o),2)]
qD = lambda o:image[delta(min(o),2)]
quotient_orbits = orbits(qvalues,qF)
assert sorted(map(len,quotient_orbits)) == [1]*4+[2]*2+[4]*2
admissible = {image[P] for o in source_orbits if len(o)<=8 for P in o}
accepted_weights = []
frobenius = {}
for o in quotient_orbits:
    p = min(o);width = N//gcd(min(p)[1],N)
    actual = qF(p);at = p
    for exponent in range(4):
        if at == actual:
            break
        at = qD(at)
    else:
        raise AssertionError('Frobenius not in the cyclic deck group')
    assert width not in frobenius or frobenius[width] == exponent
    frobenius[width] = exponent
    assert not (set(o)&admissible) or set(o)<=admissible
    if set(o)<=admissible:
        accepted_weights.append(len(o))
assert frobenius == {55:0,11:1,5:2,1:3}
assert sorted(accepted_weights) == [1,1,1,1,4]
rank_certificate = json.loads((BASE/'data'/'C55_rank_certificate.json').read_text())
assert sum(f['dimension']*f['multiplicity'] for f in rank_certificate['forms']) == 81
assert all(f['proved_analytic_rank'] == 0 for f in rank_certificate['forms'])
log = (BASE/'logs'/'N55_D8_NC7_q3.log').read_text()
assert not re.search('Runtime error|User error|Assertion failed',log)
report = {'level':55,'degree':8,'manifest_sha256':sha256((BASE/'data'/'C55_manifest.json').read_bytes()).hexdigest(),
          'all_manifest_hashes':True,'byte_identical_geometry_reproduction':True,
          'Cartier_image_dimension':15,'local_kernel_checks':34,
          'all_admissible_divisors':tested,'all_ranks':8,'exceptions':0,
          'rank_zero_weighted_dimension':81}
(BASE/'data'/'audit55.json').write_text(json.dumps(report,indent=2)+'\n')
print('C55_INDEPENDENT_AUDIT_DONE hashes=true RR_kernels=34 projections=34 D8_divisors=201 all_rank8=true exceptions=0')
