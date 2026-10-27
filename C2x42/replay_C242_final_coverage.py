#!/usr/bin/env python3
"""Check the complete finite coverage and the remaining C2x42 decisions.

The preceding independent positive-point and actual-relation HNF replays
establish 41 remaining supports. This checker recomputes every simultaneous
source automorphism witness, checks all 30 Hecke decisions and every kernel
orbit/residue-degree multiplier, and verifies the final deck symmetry.
Magma supplies the exact IsPrincipal decisions; the serialized-divisor
Magma replays independently reconstruct those obstructions and check h0.
"""
import ast
import hashlib
import json
import time
from collections import Counter
from pathlib import Path
import sys;sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'common'))
from finite_field5_checked import FF, EC, pack

P = Path(__file__).resolve().parent
t0 = time.monotonic()
bound = {}

def raw(name):
    b = (P / ('logs' if name.endswith('.log') else 'data') / name).read_bytes()
    bound[name] = hashlib.sha256(b).hexdigest()
    return b.decode()

def read(name):
    return ast.literal_eval(raw(name).replace('[*', '[').replace('*]', ']')
                            .replace('true', 'True').replace('false', 'False'))

def mul(E, Q, n):
    R = None
    while n:
        if n & 1:
            R = E.add(R, Q)
        Q = E.add(Q, Q)
        n //= 2
    return R

fields, original = read('M2N42_41_remaining_supports.m')
Fs = {d: FF(f) for d, f in fields}
old = {i: (desc, sp) for i, desc, sp in original}
assert len(old) == len(original) == 41
# Reconstruct every saved support directly from the complete source records.
masters={d:read(f'M2N42_D8_direct_full{d}_arithmetic.txt') for d in range(4,9)}
assert fields==[[d,masters[d][1]] for d in range(4,9)]
all4=read('M2N42_D8_direct_small_all_points.txt')
assert len(all4)==1 and all4[0][0]==4 and len(all4[0][2])==126
for idx,(desc,sp) in old.items():
    if desc[0]==1:
        _,d,i=desc;dd,ky=masters[d][8][i-1];assert dd==d
        expected=[[d,ky]]
    else:
        assert desc[0]==2 and desc[1]==desc[3]==4
        own=[r for r in masters[4][8] if r[0]==4]
        expected=[[4,own[desc[2]-1][1]],[4,sum(all4[0][2][desc[4]-1],[])]]
    assert sp==expected,(idx,desc)
prior = json.loads(raw('replay_M2N42_D8_quotient_filters.json'))
assert prior['three_two_isogeny_filter_survivors'] == 41
assert {i: desc for i, desc in prior['remaining']} == {i: d for i, (d, _) in old.items()}
orbits = json.loads(raw('M2N42_41_full_Aut_orbits.json'))
rf, reps = read('M2N42_full_Aut_representative_supports.m')
assert rf == fields and len(reps) == orbits['orbit_count'] == 30
assert orbits['source_Aut_order'] == 36
repdict = {i: (d, sp) for i, d, sp in reps}
assert len(repdict) == 30

def key(sp):
    return Counter((d, pack(ky[:d]), pack(ky[d:])) for d, ky in sp)

def act(sp, par):
    a, b, u = par
    assert a != b and a in range(3) and b in range(3)
    assert u in [1, 2, 4, 5, 8, 10]
    out = []
    for d, ky in sp:
        F = Fs[d]
        bb, cc = pack(ky[:d]), pack(ky[d:])
        E = EC(F, bb, cc)
        Q = (bb, bb)
        ps = [(0, 0), mul(E, Q, 21)]
        ps.append(E.add(*ps))
        pn = ps[a]
        qn = E.add(ps[b], mul(E, Q, 2*u))
        out.append((d, *F.key(E.bc(pn, qn))))
    return Counter(out)

seen = set()
for g in orbits['orbits']:
    i = g['representative_input_index']
    assert repdict[i] == (g['representative_descriptor'], g['representative_supports'])
    for member in g['members']:
        j = member['input_index']
        assert j not in seen and old[j][0] == member['descriptor']
        seen.add(j)
        assert act(old[j][1], member['automorphism_parameter_to_representative']) == key(repdict[i][1])
assert seen == set(old)

prefix = read('C242_g31_bounded_Hecke_v3_arithmetic.txt')
assert len(prefix) == 3 and len(prefix[1]) == 4 and prefix[2] is False
assert 'C242_G31_BOUNDED_HECKE_DONE 4 false' in raw('probe_C242_genus31_Hecke_v3.log')
datasets = [prefix]
for sh in range(2):
    d = read(f'C242_g31_remaining_shard{sh}_arithmetic.txt')
    assert len(d) == 4 and len(d[1]) == len(d[3]) == 13 and d[2] is True
    assert [row[:3] for row in d[1]] == [row[:3] for row in d[3]]
    assert f'C242_G31_REMAINING_SHARD_DONE {sh} 13 true' in raw(f'sieve_C242_genus31_remaining_shard{sh}.log')
    datasets.append(d)

done = {}
weights = 0
for ds in datasets:
    metadata = {}
    for ell, d, ky, e, w, ed, wd, kernels in ds[0]:
        ck = ell, d, tuple(ky)
        assert ck not in metadata
        assert d == e*w == ed*wd and sum(k for k, _, _ in kernels) == ell+1
        for k, er, wr in kernels:
            assert k > 0 and er > 0 and d*k == er*wr
            weights += 1
        assert sum(er*wr for _, er, wr in kernels) == d*(ell+1)
        metadata[ck] = True
    for i, desc, ell, principal in ds[1]:
        assert i not in done and desc == repdict[i][0]
        sp = repdict[i][1]
        ncdeg = sum(d for d, _ in sp)
        assert ncdeg in [4, 7, 8]
        assert ell == (17 if ncdeg == 4 else 11)
        assert all((ell, d, tuple(ky)) in metadata for d, ky in sp)
        done[i] = principal
assert set(done) == set(repdict)
assert [i for i, yes in done.items() if yes] == [118,193]
assert sum(not yes for yes in done.values()) == 28
assert [sum(sum(d for d, _ in sp) == n for _, sp in repdict.values()) for n in [4,7,8]] == [2,1,27]

# One survivor is pure degree8, a pair of distinct quartics, and is excluded
# by the literal forgetful deck transformation. The other needs refinement.
desc, sp = repdict[193]
assert desc == [2, 4, 5, 4, 38] and [d for d, _ in sp] == [4, 4]
assert len(key(sp)) == 2
# This parameter is literally (P2,Q) -> (P2+21Q,-Q), and [-1] identifies
# it with the deck involution (P2,Q) -> (P2+21Q,Q) of X242 -> X142.
assert act(sp, (2,1,10)) == key(sp)
for d, ky in sp:
    F = Fs[d]
    b, c = pack(ky[:d]), pack(ky[d:])
    E = EC(F,b,c)
    Q = (b,b)
    T = mul(E,Q,21)
    assert T is not None and T != (0,0) and E.add(T,T) is None
    pn = E.add((0,0),T)
    assert pn != (0,0) and E.add(pn,T) == (0,0)

def image_degree(d, ky, scalar=1):
    F=Fs[d];bb,cc=pack(ky[:d]),pack(ky[d:]);E=EC(F,bb,cc)
    Q=mul(E,(bb,bb),scalar);Q14=mul(E,Q,3);Q3=mul(E,Q,14)
    xx,yy=Q14;uu=F.div(yy,xx)
    bn=F.div(F.pow(xx,3),F.mul(yy,yy))
    cn=F.div(F.mul(cc,F.mul(xx,xx)),F.mul(yy,yy))
    zz=F.div(Q3[0],F.mul(uu,uu));dn=F.mul(bn,F.sub(F.sub(1,bn),cn))
    check=F.sub(F.add(F.add(F.sc(3,F.pow(zz,4)),F.sc(4,F.mul(cn,F.pow(zz,3)))),F.sc(6,F.mul(dn,F.mul(zz,zz)))),F.mul(dn,dn))
    assert check==0
    return F.degree([bn,cn,zz])
for ds in datasets:
    for ell,d,ky,e,w,ed,wd,kernels in ds[0]:
        assert image_degree(d,ky)==e and image_degree(d,ky,ell)==ed

extra = []
refinement = P/'data'/'C242_g31_pair118_refinement_arithmetic.txt'
assert refinement.exists()
if refinement.exists():
    ref = read(refinement.name)
    assert len(ref) == 3 and len(ref[2]) == 1
    for j, rr in enumerate(ref[1]):
        assert rr[:3] == [118,repdict[118][0],[13,17,19][j]]
    assert ref[1]==[[118,repdict[118][0],13,True],[118,repdict[118][0],17,False]]
    metadata={}
    for ell,d,ky,e,w,ed,wd,kernels in ref[0]:
        ck=ell,d,tuple(ky);assert ck not in metadata
        assert ell in [13,17] and d==e*w==ed*wd
        assert image_degree(d,ky)==e and image_degree(d,ky,ell)==ed
        assert sum(k for k,_,_ in kernels)==ell+1
        assert all(k>0 and er>0 and wr>0 and d*k==er*wr for k,er,wr in kernels)
        assert sum(er*wr for _,er,wr in kernels)==d*(ell+1)
        weights+=len(kernels);metadata[ck]=True
    assert set(metadata)=={(ell,d,tuple(ky)) for ell in [13,17] for d,ky in repdict[118][1]}
    extra = ref[1]
    assert all(r[3] for r in extra[:-1])
resolved118 = bool(extra) and extra[-1][3] is False
assert resolved118
out = {
    'status': 'complete exclusion' if resolved118 else 'complete coverage; pair118 remains unresolved',
    'source_group': [2,42], 'first_filter_inputs': 2405,
    'small_quotient_survivors': 288, 'five_map_survivors': 41,
    'full_source_Aut_order': 36, 'full_source_Aut_orbits': 30,
    'nonprincipal_genus31_initial_Hecke_cases': 28,
    'deck_excluded_index': 193, 'deck_excluded_descriptor': desc,
    'deck_excluded_support_degrees': [4,4],
    'deck_parameter': [2,1,10], 'deck_fixes_remaining_divisor': True,
    'deck_global_obstruction': 'rank0 and odd5 injectivity lift h-invariance; gonality>8 makes it literal; forgetting P2 gives an order42 point over degree dividing4, excluded by DN Theorem1.1',
    'pair118_refinement': extra,
    'remaining': 0 if resolved118 else 1, 'checked_kernel_orbit_weights': weights,
    'exact_Hecke_results': [[i,done[i]] for i in sorted(done)],
    'verification_boundary': __doc__, 'sha256': bound,
    'seconds': round(time.monotonic()-t0,3)
}
(P/'data'/'C242_final_coverage.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
print('C242_FINAL_COVERAGE_DONE 2405 288 41 30 28 2',out['remaining'])
