"""Independent F5 arithmetic replay of all 7,862 source-to-target Aut16 keys."""
import hashlib
import json
import runpy
from pathlib import Path
from replay_source_partitions import Field, Curve

HERE = Path(__file__).resolve().parent / 'data'
source = HERE / 'source_partition_inputs_complete.py'
orbit_file = HERE / 'clean_run/high_target_aut_orbits.py'
a = runpy.run_path(str(source))
o = runpy.run_path(str(orbit_file))
assert o['complete'] is True
assert o['master_sha'] == hashlib.sha256(source.read_bytes()).hexdigest()
assert o['target_automorphisms'] == [[e,h,u] for e in (0,1) for h in (0,1) for u in (1,5,7,11)]

def target(F, raw):
    b,c = map(F.encode, raw)
    E = Curve(F,b,c)
    r,ty = E.multiple(24)
    assert ty == 0
    x,y = E.P
    x = F.sub(x,r)
    A = F.add(E.A2,F.mul(3,r))
    B = F.add(F.add(F.mul(3,F.mul(r,r)),F.mul(2,F.mul(E.A2,r))),E.A4)
    Z = object.__new__(Curve)
    Z.F = F
    Z.A2 = F.neg[F.mul(2,A)]
    Z.A4 = F.sub(F.mul(A,A),F.mul(4,B))
    Z.A6 = 0
    Z.P = (F.div(F.mul(y,y),F.mul(x,x)), F.div(F.mul(y,F.sub(B,F.mul(x,x))),F.mul(x,x)))
    assert Z.on_curve(Z.P) and Z.on_curve((0,0))
    assert Z.multiple(24) is None and Z.multiple(12) not in (None,(0,0)) and Z.multiple(8) is not None
    return Z

def normalized(F,Z,P,Q):
    r,py = P
    assert py == 0
    x,y = Q
    x = F.sub(x,r)
    b = F.div(F.mul(x,F.mul(x,x)),F.mul(y,y))
    c = F.div(F.mul(F.add(Z.A2,F.mul(3,r)),F.mul(x,x)),F.mul(y,y))
    d = F.div(F.mul(F.add(F.add(Z.A4,F.mul(2,F.mul(Z.A2,r))),F.mul(3,F.mul(r,r))),F.mul(F.mul(x,x),F.mul(x,x))),F.mul(F.mul(y,y),F.mul(y,y)))
    assert d == F.mul(b,F.sub(F.sub(1,b),c))
    return b,c

def frobenius_key(F,b,c):
    # Magma compares low-first coefficient lists lexicographically. This
    # also permits target points in the proper degree4 subfield of F5^8.
    keys = []
    for _ in range(F.d):
        keys.append(F.digits[b]+F.digits[c])
        b,c=F.frob(b),F.frob(c)
    return min(keys)

seen, representatives, reports = {}, [], {}
records = iter(o['cases'])
for d in (6,7,8):
    F=Field(a['field_polynomials'][d])
    before=len(seen)
    for idx,raw in enumerate(a['representative_pairs'][d],1):
        Z=target(F,raw)
        q2=Z.multiple(12)
        keys=[]
        for eps,eta,u in o['target_automorphisms']:
            P=Z.add((0,0),q2) if eps else (0,0)
            Q=Z.add(Z.multiple(u),(0,0)) if eta else Z.multiple(u)
            keys.append(frobenius_key(F,*normalized(F,Z,P,Q)))
        key=(d,)+min(keys)
        if key not in seen:
            seen[key]=len(seen)+1
            representatives.append([d,idx])
        row=next(records)
        assert row == [d,idx,seen[key],list(key)], (d,idx,row,key)
        if idx % 1000 == 0:
            print('C48_TARGET_AUT_INDEPENDENT_PROGRESS',d,idx,flush=True)
    reports[d]={'source_representatives':len(a['representative_pairs'][d]),'target_clusters':len(seen)-before}
    print('C48_TARGET_AUT_INDEPENDENT_PASS',d,reports[d],flush=True)
assert next(records,None) is None
assert representatives == o['representatives']
assert len(seen)==3507 and sum(r['source_representatives'] for r in reports.values())==7862
result={'complete':True,'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
        'orbit_export_sha256':hashlib.sha256(orbit_file.read_bytes()).hexdigest(),
        'all_7862_keys_independently_recomputed':True,'target_clusters':3507,'degrees':reports}
(HERE/'clean_run/high_target_aut_replay.json').write_text(json.dumps(result,indent=2)+'\n')
print('C48_TARGET_AUT_INDEPENDENT_ALL_PASS',flush=True)
