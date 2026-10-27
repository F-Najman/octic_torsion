"""Independent lattice/weight/coverage replay of the completed high filter."""
import ast
import hashlib
import json
import runpy
from collections import Counter
from pathlib import Path
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form
from replay_source_partitions import Field

P=Path(__file__).resolve().parent/'data';D=P/'clean_run_v2'
setup=runpy.run_path(str(P/'clean_run/target_2_24_arithmetic.py'))  # identical to the clean_run_v2 copy, merged
orbit=runpy.run_path(str(P/'clean_run/high_target_aut_orbits.py'))
master=runpy.run_path(str(P/'source_partition_inputs_complete.py'))
assert orbit['complete'] is True
assert orbit['master_sha']==hashlib.sha256((P/'source_partition_inputs_complete.py').read_bytes()).hexdigest()
R=Matrix(setup['relations']);n=8
assert R.shape==(n,n+1) and R[:,n]==Matrix.zeros(n,1)
H=hermite_normal_form(R[:,:n].T.row_join(Matrix([v[:n] for v in setup['rational_cusp_generators']]).T))
assert abs(int(H.det()))==4608
def red(v):
    assert len(v)==n+1 and v[n]==0
    z=list(v[:n])
    for i in reversed(range(n)):
        q=z[i]//int(H[i,i])
        for j in range(i+1):z[j]-=q*int(H[j,i])
    return tuple(z)
def add(v,w):return red([x+y for x,y in zip(v,w)]+[0])
comp=[{(0,)*n}]+[set() for _ in range(8)]
for degree,v in zip(setup['finite_cusp_degrees'],setup['finite_cusp_classes']):
    for j in range(degree,9):
        comp[j].update(add(z,red(v)) for z in list(comp[j-degree]))
assert list(map(len,comp))==[1,5,23,75,178,378,661,1065,1487]
rows={};shas={}
for shard in range(8):
    path=D/f'high_shard{shard}_arithmetic.py'
    a={x.targets[0].id:ast.literal_eval(x.value) for x in ast.parse(path.read_text()).body}
    assert a['complete'] is True and a['shard']==shard and a['master_sha']==orbit['master_sha']
    assert a['relations']==setup['relations'] and a['rational_cusp_generators']==setup['rational_cusp_generators']
    ids=[cid for cid in range(1,3508) if (cid-1)%8==shard]
    assert a['expected_cluster_ids']==ids==[r[0] for r in a['tests']]
    for row in a['tests']:
        cid,d,idx,e,v,passed=row
        assert cid not in rows and [d,idx]==orbit['representatives'][cid-1]
        assert d%e==0 and d//e in (1,2) and e>0
        assert passed==int(red([-z for z in v]) in comp[8-d])
        rows[cid]=row
    shas[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
assert sorted(rows)==list(range(1,3508))
keys={r[2]:r[3] for r in orbit['cases']}
for d in (6,7,8):
    F=Field(master['field_polynomials'][d])
    for cid,row in rows.items():
        if row[1]!=d:continue
        key=keys[cid];assert key[0]==d and len(key)==1+2*d
        b,c=F.encode(key[1:1+d]),F.encode(key[1+d:])
        bb,cc=b,c
        for e in range(1,d+1):
            bb,cc=F.frob(bb),F.frob(cc)
            if (bb,cc)==(b,c):break
        assert e==row[3]
survivors=[r[:3] for r in orbit['cases'] if rows[r[2]][-1]]
assert Counter(r[0] for r in survivors)=={6:4,7:2,8:2}
assert survivors==[[6,47,46],[6,85,79],[6,171,121],[6,243,147],
                   [7,348,453],[7,1102,453],[8,2286,2486],[8,2705,2699]]
assert len(orbit['cases'])==7862
out={'complete':True,'class_tests':3507,'source_cases_covered':7862,
     'all_3507_cluster_ids_unique_and_exhaustive':True,'proper_pushforward_weights_recomputed':True,
     'all_cusp_complements_included':True,'target_Aut_truth_values_expanded_to_all_source_cases':True,
     'passing_source_cases':survivors,'arithmetic_sha256':shas,'exclusion_claimed':False}
(D/'high_filter_replay.json').write_text(json.dumps(out,indent=2)+'\n')
print('C48_HIGH_FILTER_INDEPENDENT_PASS clusters=3507 source_cases=7862 admitted=8')
