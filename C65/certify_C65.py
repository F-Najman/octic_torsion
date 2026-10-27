#!/usr/bin/env python3
"""Check the exact completed octic C65 evidence and optional manifest."""
import ast
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

p=Path(__file__).resolve().parent
def loc(name):
    # the files in the layout of this repository
    if name.startswith('C65_inputs/'):
        return p.parent/'common'/'logs'/name.split('/')[-1]
    if name.endswith('.log'):
        return p/'logs'/name
    if name=='C65_genus61_rr_kernels.m':
        return p/'data'/'C65_root_allbad_rr_kernels.m'
    return p/'data'/name
def read(name):
    return loc(name).read_text()
def no_errors(name):
    s=read(name)
    assert not re.search(r'Runtime error|System error|Assertion failed|FAILED',s),name
    return s

counts=read('C65_inputs/independent_counts.log')
line=[s for s in counts.splitlines() if 'm=1 n=65 p=2 ' in s]
assert len(line)==1 and 'closed_place_counts=[ 0, 0, 0, 0, 0, 2, 48, 48 ]' in line[0]
assert 'DEGREE8_FAST_COUNTS_DONE complete=true' in counts
tate=json.loads(read('cusp65_genus61.json'))
assert tate['F2_rational_source_cusps_are_Q_rational']
assert tate['admissible_degree8_divisors']==83949 and tate['source_RR_blocks']==11

for e,pairs,degrees in [(7,56,'[ 28 ]'),(8,64,'[ 16, 16 ]')]:
    s=no_errors(f'C65_good{e}_D8.log')
    assert f'OCTIC65_GOOD{e}_EXHAUSTIVE pairs={pairs} closed_places=8 diamond_orbits=2' in s
    assert s.count(f'OCTIC65_GOOD{e}_A3_DIFFERENCE')==2
    assert s.count('principal=false')==2 and 'principal=true' not in s
    assert s.count(f'image_degrees={degrees}')==2
    assert f'OCTIC65_GOOD{e}_A3_DONE' in s
s=no_errors('C65_good6_actual13_exhaustive.log')
assert 'generators=48 input_z_infinity=48 isogeny_images=192 finite_z=192 infinity_z=0 full_source_degrees={ 12 }' in s
s=no_errors('C65_good6_D8.log')
assert 'cases=30 s_tests=30 u_tests=2 hecke2_tests=2 survivors=0 data=[] complete=true' in s
records=[ast.literal_eval(s) for s in read('C65_good6_D8_records.txt').splitlines()]
expected={(i,a,b) for i in range(1,4) for a in range(6,10) for b in range(a,10)}
assert len(records)==30 and {tuple(r[:3]) for r in records}==expected
assert all(len(r)==6 and r[3:6] in ([0,0,0],[1,1,0]) for r in records)
assert sum(r[3] for r in records)==sum(r[4] for r in records)==2
assert not any(all(r[3:6]) for r in records)

s=no_errors('C65_model_retained_completed.log')
assert 'irreducible12=true genus61=true index1008=true' in s
assert 'joint_image=48' in s
s=no_errors('C65_allbad_retained_completed.log')
assert s.count('C65_RR_SOURCE_DONE')==11 and s.count('C65_RR_DIRECT_TARGET_CHECK')==3
assert 'C65_GENUS61_FI_DONE tested=83949 total=83949 exceptions=0 complete=true' in s
original=read('C65_genus61_preflight_v3.log')
assert original.count('Runtime error')==1 and 'No valid universe containing all elements' in original
assert 'C65_GENUS61_RETAINED_MATRICES_EXPORTED' in original and 'Total time: 1501.190 seconds' in original
assert no_errors('winding65_hecke2_certificate.log').count('OCTIC65_WINDING')==12
assert 'OCTIC65_FINITE_IMAGE_T2_PLUS1_AND_SIXTH_POWER_AND_DIAMONDS_CERTIFIED 121 [ 1 ]' in read('winding65_hecke2_certificate.log')

for script,args in [('verify_winding65_hecke2.py',['winding65_hecke2_witness.m','winding65_hecke2_certificate.log']),
                    ('verify_C65_genus61.py',['C65_genus61_rr_kernels.m','cusp65_genus61.json'])]:
    command=[sys.executable,str(p/script)]+[str(loc(a)) for a in args]
    completed=subprocess.run(command,check=True,text=True,capture_output=True,cwd=p.parent)
    print(completed.stdout.strip())
manifest=p/'C65_manifest.json'
if manifest.exists():
    for entry in json.loads(manifest.read_text())['files']:
        f=p/entry['path']
        assert hashlib.sha256(f.read_bytes()).hexdigest()==entry['sha256'],str(f)
print('OCTIC_C65_COMPLETE: all four reduction branches excluded; exact rank121 support; retained and independent replay checks passed.')
