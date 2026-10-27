#!/usr/bin/env python3
"""Independent checksum and six-support-type coverage audit for C46."""
import ast
from hashlib import sha256
import json
from pathlib import Path
import re

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'C46'
data=json.loads((BASE/'data'/'N46_combined_audit.json').read_text())
assert data['complete'] and data['survivors']==0
counts=[0,0,0,55,44,55,264,781]
assert data['closed_place_counts']==counts and data['independent_counts_agree']
assert 9406793==23*408991
logs={}
for item in data['details']:
    path=BASE/'logs'/item['log'];assert sha256(path.read_bytes()).hexdigest()==item['sha256']
    text=path.read_text();assert not re.search(r'Runtime error|User error|Assertion failed|ERROR ',text)
    logs[item['log']]=text
single=logs['N46_p3_single_completed.log'];mixed=logs['N46_p3_mixed_quotient.log']
expect={4:5,5:4,6:5,7:24,8:71}
for d,n in expect.items():
    entries=re.findall(rf'SMART_CASE_RESULT type=\[ {d} \] index=(\d+) killed=(\w+)',single)
    assert entries==[(str(i),'true') for i in range(1,n+1)]
    assert f'RESULT type=[ {d} ] total={n} survivors=0' in single
assert 'SINGLE_BRANCH_CHECKPOINT_DONE level=46 p=3 D=8 completed_degrees=4,5,6,7,8' in single
indices=list(map(int,re.findall(r'QUOTIENT_REJECT type=\[ 4, 4 \] index=(\d+) base=23',mixed)))
assert len(indices)==len(set(indices))==263 and set(indices)<=set(range(1,276))
assert 'RESULT type=[ 4, 4 ] total=275 survivors=0' in mixed
assert 'SUMMARY N=46 M=1 p=3 D=8 shard=0/1 total_survivors=0' in mixed
assert 'classgroup=[ 9406793, 0 ] rational_cusp_image=[ 408991 ] quotient=[ 23, 0 ]' in mixed
assert 'degrees=[ '+', '.join(['1']*22+['11']*2)+' ]' in mixed
assert sum(expect.values())==109 and 5*55==275 and 109+275==384
def partitions(n,m=None):
    if n==0:yield ();return
    for a in range(min(n,m or n),0,-1):
        for rest in partitions(n-a,a):yield (a,)+rest
assert sum(sum(1 for _ in partitions(n)) for n in range(1,9))==66
possible={a for n in range(1,9) for a in partitions(n) if all(counts[d-1] for d in a)}
assert possible=={(4,),(5,),(6,),(7,),(8,),(4,4)}
paths=['common/quotient_class_prefilter.m',
       'common/class_rational_reduction.m','common/hecke_sieve_quotient.m',
       'C46/logs/N46_p3_single_completed.log','C46/logs/N46_p3_mixed_quotient.log','C46/data/N46_combined_audit.json']
bindings=[{'path':p,
           'sha256':sha256((ROOT/p).read_bytes()).hexdigest()} for p in paths]
shared=ROOT/'mdmagma/v2/modular_curves/X1.m'
bindings.append({'path':str(shared.relative_to(ROOT)),'sha256':sha256(shared.read_bytes()).hexdigest()})
report={'level':46,'degree':8,'independent_math_audit':'passes',
        'source_singles':109,'mixed_pairs':275,'quotient_rejections':263,
        'remaining_pairs_killed_by_Hecke':12,'total_tests':384,
        'nonempty_support_types':sorted(map(list,possible)),
        'all_recorded_hashes_pass':True,'files':bindings}
(BASE/'data'/'audit46.json').write_text(json.dumps(report,indent=2)+'\n')
print('C46_INDEPENDENT_AUDIT_DONE theory=passes source_singles=109 mixed_pairs=275 quotient_rejections=263 remaining_Hecke_kills=12 total=384 survivors=0 hashes=true')
