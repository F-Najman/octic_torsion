#!/usr/bin/env python3
"""Remove an exactly repeated initialization from the raw executed export.

The live source is preserved. Its explicit JetFile argument caused the
initial metadata and Cartier matrices to be written twice, identically.
No kernel, restriction, or mathematical input is changed by this step.
"""
from pathlib import Path
import hashlib,json,re
H=Path(__file__).resolve().parent/'data'
raw=H/'N63_cartier_kernels.m'
s=raw.read_text()
log=(H.parent/'logs'/'N63_cartier_rr.log').read_text()
assert not re.search(r'Runtime error|Assertion failed|SystemError|User error',log)
assert 'RANK0_CARTIER_DONE N=63 D=8 tested=44046 total=44046 exceptions=0 complete=true' in log
marker='// Direct differential-space restriction certificate.'
assert s.count(marker)==1 and s.count('print "FI_JET_DATA_READY";')==1
k=s.index(marker);prefix=s[:k];tail=s[k:]
for name in ['N','As','D','gx','o','udim','JetFactor','ds','perms']:
    a=re.findall(rf'\b{name}:=(.*?);',prefix,re.S)
    b=re.findall(rf'\b{name}:=(.*?);',tail,re.S)
    assert len(a)==len(b)==1 and a==b,name
for name in ['CartierMatrix','ImageBasis']:
    a=re.findall(rf'\b{name}:=Matrix\(GF\(2\),(.*?)\);',prefix,re.S)
    b=re.findall(rf'\b{name}:=Matrix\(GF\(2\),(.*?)\);',tail,re.S)
    assert len(a)==len(b)==1 and a==b,name
assert not re.search(r'kernelspaces\[|basisjets\[|localjets\[',prefix)
assert len(re.findall(r'kernelspaces\[<',tail))==103
assert len(re.findall(r'localjets\[<',tail))==103
assert len(re.findall(r'basisjets\[',tail))==22
out=H/'N63_cartier_kernels_canonical.m'
out.write_text('// Same mathematical export; duplicate initialization removed.\n'+tail)
record={'raw':raw.name,'raw_sha256':hashlib.sha256(s.encode()).hexdigest(),
 'canonical':out.name,'canonical_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),
 'discarded_prefix_sha256':hashlib.sha256(prefix.encode()).hexdigest(),
 'identical_repeated_metadata_checked':True,'identical_Cartier_and_image_matrices_checked':True,
 'discarded_prefix_contains_no_restriction_data':True,
 'all103kernels_and_restrictions_preserved':True}
(H/'N63_export_normalization.json').write_text(json.dumps(record,indent=2)+'\n')
print('OCTIC63_EXPORT_NORMALIZED all103kernels_preserved=true')
