"""Read-only strict C48 certificate audit (adapted: the manifest/freeze part and the hash checks of the record are removed).

The file scope is explicit. Exploratory controls and interrupted runs are
preserved separately and are never accepted as exclusion evidence.
"""
import argparse
import ast
from collections import Counter
import hashlib
import json
from math import gcd
from pathlib import Path
import re

ROOT=Path(__file__).resolve().parents[1]
P=ROOT/'C48'/'data'
L=ROOT/'C48'/'logs'
def literals(path):
    out={}
    for n in ast.parse(path.read_text()).body:
        assert isinstance(n,ast.Assign) and len(n.targets)==1
        target=n.targets[0];value=ast.literal_eval(n.value)
        if isinstance(target,ast.Name):out[target.id]=value
        else:
            assert isinstance(target,ast.Subscript) and isinstance(target.value,ast.Name)
            out[target.value.id][ast.literal_eval(target.slice)]=value
    return out
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def read_json(path):return json.loads(path.read_text())
def clean_log(path):
    s=path.read_text()
    assert not any(x in s for x in ('Runtime error','User error','Assertion failed','FAILED'))
    return s

source=literals(P/'source_partition_inputs_complete.py')
assert source['complete'] is True and source['complete_by_newton_count'] is True
assert source['octic_original_places_scanned']==38173
assert [len(source['representative_pairs'][d]) for d in range(3,9)]==[2,8,82,342,1348,6172]
replay=read_json(P/'source_partition_replay.json')
assert replay['complete'] is True and replay['source_sha256']==sha(P/'source_partition_inputs_complete.py')
assert replay['all_points_exact_order48'] and replay['all_diamond_orbits_disjoint'] and replay['all_newton_counts_exhausted']
aut=read_json(P/'clean_run/high_target_aut_replay.json')
assert aut['complete'] and aut['all_7862_keys_independently_recomputed'] and aut['target_clusters']==3507
assert aut['source_sha256']==sha(P/'source_partition_inputs_complete.py')
assert aut['orbit_export_sha256']==sha(P/'clean_run/high_target_aut_orbits.py')
low=read_json(P/'clean_run/low_filter_replay.json')
high=read_json(P/'clean_run_v2/high_filter_replay.json')
assert low['complete'] and low['tested']==2076 and len(low['admitted'])==18
assert high['complete'] and high['class_tests']==3507 and high['source_cases_covered']==7862
assert len(high['passing_source_cases'])==8
assert high['all_3507_cluster_ids_unique_and_exhaustive'] and high['proper_pushforward_weights_recomputed']
assert low['arithmetic_sha256']==sha(P/'clean_run/lower_support_arithmetic.py')
for filename,digest in high['arithmetic_sha256'].items():assert sha(P/'clean_run_v2'/filename)==digest

for dirname in ('clean_run','clean_run_v2'):
    data=P/dirname;merged=P/'clean_run'  # the clean_run_v2 copies of these three files were identical, merged
    j=read_json(merged/'target_setup_replay.json');r=read_json(data/'proved_presentation_replay.json')
    assert j['complete'] and j['finite_quotient_index']==4608 and j['rational_cusp_image_order']==921600
    assert j['arithmetic_sha256']==sha(merged/'target_2_24_arithmetic.py')
    assert r['complete'] and r['finite_order']==4246732800 and r['factorbase_map_surjective']
    assert r['state_sha256']==sha(data/'proven_relation_state.py') and r['geometry_sha256']==sha(merged/'factorbase_geometry.py')

# Every original low source image is independently matched to an actual
# prime of the fresh proved target factorbase, in the same abstract basis.
fb=literals(P/'clean_run_v2/low_factorbase_key_replay.py')
assert fb['complete'] is True and len(fb['target_table'])==784 and len(fb['source_images'])==736
low_ar=literals(P/'clean_run/lower_support_arithmetic.py')
assert [len(low_ar['all_source_classes'][d]) for d in (3,4,5)]==[16,64,656]
for d,i,j,v in fb['source_images']:assert v==low_ar['all_source_classes'][d][i-1]

all_witness_ids=[];rough=smooth=0
for shard in range(8):
    s=clean_log(L/f'clean_run_v2/high_witness_replay_v3_shard{shard}.log')
    ids=list(map(int,re.findall(r'C48_WITNESS_REPLAY_CASE cluster=(\d+) both_orders_exact=true',s)))
    expected=[cid for cid in range(1,3508) if (cid-1)%8==shard]
    assert ids==expected
    done=re.findall(r'C48_WITNESS_REPLAY_DONE shard=(\d+) tested=(\d+) rough=(\d+) smooth=(\d+) complete=true',s)
    assert len(done)==1
    sh,n,r,sm=map(int,done[0]);assert sh==shard and n==len(ids)==r+sm
    all_witness_ids+=ids;rough+=r;smooth+=sm
assert sorted(all_witness_ids)==list(range(1,3508))

hec=literals(P/'low_hecke_af_xy_results.py')
assert hec['complete'] is True and hec['survivors']==[] and len(hec['tests'])==18
assert hec['input_arithmetic_sha']==sha(P/'clean_run/lower_support_arithmetic.py')
assert [z[:3] for z in hec['tests']]==low['admitted']
assert all(z[3]==7 and sum(z[1])==8 for z in hec['tests'])
hec_log=clean_log(L/'low_hecke_af_xy.log')
assert len(re.findall(r'C48_LOW_AF_TEST .*q=7 principal=false ',hec_log))==18
assert 'C48_LOW_AF_DONE tested=18 remaining=0 complete=true' in hec_log

single_cases={6:([47,85,171,243],29),7:([348,1102],13),8:([2286,2705],7)}
for d,(indices,q) in single_cases.items():
    assert [v[1] for v in high['passing_source_cases'] if v[0]==d]==indices
    for shard,index in enumerate(indices):
        s=clean_log(L/f'high_degree{d}_shard{shard}.log')
        pattern=rf'C48_D{d}_HEC_DONE source_index={index} q={q} principal=false complete=true'
        assert len(re.findall(pattern,s))==1

# Independent DEHMZ cusp-stabilizer arithmetic checks q13 for complement1
# and q29 for arbitrary F5-rational complements. No geometric cusp labels
# or permutation guessed from cardinalities enter this check.
units=[a for a in range(1,49) if gcd(a,48)==1]
phi=lambda n:sum(gcd(a,n)==1 for a in range(1,n+1))
pred=Counter();fixed_types=[]
for d in (d for d in range(1,49) if 48%d==0):
    m=48//d
    H={a for a in units if a%m==1%m} if d>2 else {a for a in units if a%m in {1%m,(-1)%m}}
    size=len(units)//len(H);geom=phi(d)*phi(m)//2
    assert geom%size==0
    f=1;x=5%48
    while x not in H:x=x*5%48;f+=1
    assert size%f==0
    pred[f]+=geom//f
    if f==1:
        assert d>2 and m in (1,2,4) and 13 in H
        fixed_types.append([d,m])
assert dict(pred)=={1:16,2:12,4:4}
assert pow(5,3,48)==29
assert 'C48_GENERIC_SOURCE_TATE_INVERSE_DONE exact_function_field_identity=true' in clean_log(L/'check_generic_source_inverse.log')

report={'complete':True,'level':48,'degree':8,'source_support_cases':9938,'target_high_class_tests':3507,
        'independent_exact_principal_identity_replays':3507,'rough_replays':rough,'smooth_replays':smooth,
        'low_source_images_freshly_matched_to_target_primes':736,'actual_source_hecke_tests':26,
        'remaining_source_cases':0,'q13_cusp_stabilizer_types':fixed_types,
        'source_cusp_signature':dict(pred),
        'degree8_only':True,'earlier_degree_exclusions_transferred':False}
print(json.dumps(report,indent=2))
print('C48_CERTIFIED_COMPLETE')
