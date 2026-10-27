"""Independent integer-lattice and exact-support replay of all 2,076 tests."""
from collections import Counter
import hashlib
import json
import runpy
from pathlib import Path
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form
from replay_source_partitions import Field, Curve

P=Path(__file__).resolve().parent/'data'
path=P/'clean_run/lower_support_arithmetic.py'
a=runpy.run_path(str(path))
setup=runpy.run_path(str(P/'clean_run/target_2_24_arithmetic.py'))
source=runpy.run_path(str(P/'source_partition_inputs_complete.py'))
assert a['complete'] is True
for name in ('relations','rational_cusp_generators','finite_cusp_classes','finite_cusp_degrees'):
    assert a[name] == setup[name]
R=Matrix(a['relations'])
assert R.shape==(8,9) and R[:,8] == Matrix.zeros(8,1)
H=hermite_normal_form(R[:,:8].T.row_join(Matrix([v[:8] for v in a['rational_cusp_generators']]).T))
assert abs(int(H.det()))==4608
def red(v):
    assert len(v)==9 and v[8]==0
    z=list(v[:8])
    for i in reversed(range(8)):
        q=z[i]//int(H[i,i])
        for j in range(i+1):
            z[j]-=q*int(H[j,i])
    return tuple(z)
def add(v,w):
    return red([x+y for x,y in zip(v,w)]+[0])
comp=[{(0,)*8}]+[set() for _ in range(8)]
for degree,v in zip(a['finite_cusp_degrees'],a['finite_cusp_classes']):
    v=red(v)
    for j in range(degree,9):
        comp[j].update(add(t,v) for t in list(comp[j-degree]))
assert list(map(len,comp))==[1,5,23,75,178,378,661,1065,1487]

for d,count in ((3,16),(4,64),(5,656)):
    F=Field(source['field_polynomials'][d])
    full={}
    for raw,v in zip(a['all_source_pairs'][d],a['all_source_classes'][d]):
        b,c=map(F.encode,raw)
        key=F.key(b,c)
        assert key not in full
        full[key]=v
    assert len(full)==count
    seen=set()
    for raw,v in zip(source['representative_pairs'][d],a['representative_classes'][d]):
        b,c=map(F.encode,raw)
        assert full[F.key(b,c)]==v
        E=Curve(F,b,c)
        for unit in (1,5,7,11,13,17,19,23):
            seen.add(F.key(*E.tate_pair(E.multiple(unit))))
    assert set(full)==seen

expected=[]
for d in (3,4,5):
    for i,v in enumerate(a['representative_classes'][d],1):
        expected.append(([d],[i],v))
for d,e in ((3,3),(3,4),(3,5),(4,4)):
    for i,v in enumerate(a['representative_classes'][d],1):
        for j,w in enumerate(a['all_source_classes'][e],1):
            expected.append(([d,e],[i,j],[x+y for x,y in zip(v,w)]))
assert len(expected)==len(a['tests'])==2076
counts=Counter();admitted=[]
for row,(saved,(degs,idx,v)) in enumerate(zip(a['tests'],expected),1):
    assert saved[:2]==[degs,idx]
    assert red(saved[2])==red(v)
    passed=red([-z for z in v]) in comp[8-sum(degs)]
    assert saved[3]==int(passed)
    counts[tuple(degs)]+=1
    if passed:
        admitted.append([row,degs,idx])
assert len(admitted)==18
assert Counter(tuple(z[1]) for z in admitted)=={(3,5):2,(4,4):16}
out={'complete':True,'arithmetic_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
     'tested':2076,'support_counts':{','.join(map(str,k)):v for k,v in counts.items()},
     'admitted':admitted,'every_finite_cusp_complement_included':True,
     'all_second_supports_including_repetitions_checked':True,'exclusion_claimed':False}
(P/'clean_run/low_filter_replay.json').write_text(json.dumps(out,indent=2)+'\n')
print('C48_LOW_FILTER_LATTICE_SUPPORT_REPLAY_PASS tested=2076 admitted=18')
