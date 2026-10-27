#!/usr/bin/env python3
"""Check completion, coverage and exact integer arithmetic for octic C54."""
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
p=Path(__file__).resolve().parent
def log(name):
    s=(p/'logs'/name).read_text()
    assert not re.search(r'Runtime error|System error|Assertion failed|FAILED|Interrupted',s),name
    return s
s=log('winding54_octic.log')
assert 'weighted_dimension=52 bad_orders=[]' in s
assert s.count('OCTIC54_WINDING order=')==3
s=log('C54_via27_preflight.log')
assert 'OCTIC54_F2_CANDIDATES 531' in s and 'cuspidal_relations_equal=true' in s
assert 'OCTIC54_LOWER_SUPPORTS_DONE checked=54849 members=144 passing=0' in s
assert 'OCTIC54_DEGREE7_REVALIDATED representatives=828 weighted_source_places=11412 checked=62301 members=153 passing=0 complete=true' in s
assert 'OCTIC54_CLASS_ARITHMETIC_EXPORTED tests=62301 members=153' in s
assert 'OCTIC54_QUARTIC_SPLITTING 126 45' in s
s=log('C54_good8_enumeration_v2.log')
assert 'representatives=3579 weighted_source_places=47961 other_degree4_lifts=126 j_orbits=48915 curves=97836' in s
assert 'proper_j_records=0' in s and 'complete=true' in s
s=log('C54_good8_positive_verify_v2.log')
assert 'primitive27=true rational2counts=true nine_distinct_diamonds=true all_orbits_disjoint=true weighted=47961 roots1=2704 roots3=875 complete=true' in s
for shard,n in enumerate([895,895,895,894]):
    s=log(f'C54_good8_class_shard{shard}.log')
    assert f'OCTIC54_GOOD8_CLASS_DONE shard={shard} tested={n} members=0 passing=0 complete=true' in s
counts=(p.parent/'common'/'logs'/'independent_counts.log').read_text()
lines=[v for v in counts.splitlines() if 'm=1 n=54 p=5 ' in v]
assert len(lines)==1 and 'closed_place_counts=[ 0, 0, 54, 261, 612, 3033, 11412, 48087 ]' in lines[0]
assert 9*(2704+3*875)+126==48087
assert sum([7344,12084,8448,10176,7452,3579,48,5760,6156,1584,3249])==65880
run=subprocess.run([sys.executable,str(p/'replay_C54_classes.py')],check=True,text=True,capture_output=True)
print(run.stdout.strip())
manifest=p/'C54_manifest.json'
if manifest.exists():
    for entry in json.loads(manifest.read_text())['files']:
        q=p/entry['path']
        assert hashlib.sha256(q.read_bytes()).hexdigest()==entry['sha256'],str(q)
print('OCTIC_C54_COMPLETE: semistable2/odd5 reductions, fullsource coverage and all65880 integer-class tests passed.')
