#!/usr/bin/env python3
"""Strict octic certificate coverage, exact input and gcd/hash audit."""
from fractions import Fraction
import hashlib
import json
from math import gcd
from pathlib import Path
import re
import sympy as s

HERE=Path(__file__).resolve().parent
data=json.loads((HERE/'data/vanhoeij_octic_input.json').read_text())
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(HERE/'data/LowDegreePlaces')==data['source_sha256']
examples=data['examples'];expected=[33,34,34,34,35,39,40,42,42,42]
assert [e['N'] for e in examples]==expected
t,x,y=s.symbols('t x y')
def poly(co):
    return sum(s.Rational(c)*t**i for i,c in enumerate(co))
source_lines=(HERE/'data/LowDegreePlaces').read_text().splitlines()
for i,e in enumerate(examples,1):
    assert e['index']==i
    assert source_lines[e['source_line']-1]==e['source_text']
    f=poly(e['f']);xp=poly(e['x']);yp=poly(e['y'])
    assert s.degree(f,t)==8 and s.Poly(f,t).is_irreducible
    match=re.search(r'\[(.*) = 0, (.*) = 0\]',e['source_text'])
    for eq in match.groups():
        z=s.sympify(eq.replace('^','**')).subs({x:xp,y:yp})
        assert s.rem(z,f,t)==0

log=(HERE/'logs/verify_octic_examples.log').read_text()
assert not re.search(r'Runtime error|User error|Assertion failed',log)
assert 'OCTIC_VANHOEIJ_EXAMPLES_DONE count=10 complete=true' in log
full=re.findall(r'OCTIC_FULL_TORSION_CERTIFIED index=(\d+) invariants=\[ (\d+) \] exact=true',log)
assert full==[(str(i),str(n)) for i,n in enumerate(expected,1)]
records=re.findall(r'OCTIC_GOOD_REDUCTION index=(\d+) p=(\d+) residue_degree=(\d+) order=(\d+) gcd=(\d+)',log)
bounds={i:0 for i in range(1,11)}
for idx,p,d,order,reported in records:
    idx,p,d,order,reported=map(int,(idx,p,d,order,reported))
    assert s.isprime(p) and p>2 and 1<=d<=8 and order%expected[idx-1]==0
    bounds[idx]=gcd(bounds[idx],order);assert bounds[idx]==reported
    f=poly(examples[idx-1]['f']);disc=s.discriminant(f,t)
    assert int(s.numer(disc))%p and int(s.denom(disc))%p
    assert all(Fraction(c).denominator%p for key in ('f','x','y') for c in examples[idx-1][key])
assert [bounds[i] for i in range(1,11)]==expected
isog=(HERE/'logs/verify_octic_two_isogenies.log').read_text()
assert 'OCTIC_TWO_ISOGENIES_DONE targets=7 complete=true' in isog
assert isog.count('OCTIC_TWO_ISOGENY_CERTIFIED')==7
assert not re.search(r'Runtime error|User error|Assertion failed',isog)
names=['audit_octic_realizations.py','extract_octic_source.py',
       'data/LowDegreePlaces','data/vanhoeij_octic_input.json',
       'data/vanhoeij_octic_input.m','logs/verify_octic_examples.log',
       'verify_octic_examples.m','logs/verify_octic_two_isogenies.log',
       'verify_octic_two_isogenies.m']
if (HERE/'logs/verify_octic_examples_root.log').exists():
    names.append('logs/verify_octic_examples_root.log')
files=[HERE/n for n in names]
assert all(p.is_file() for p in files)
report={'degree':8,'complete':True,'new_cyclic_groups':[33,34,35,39,40,42],
        'examples':10,'independent_symbolic_input_checks':True,
        'exact_field_degrees':[8]*10,'full_torsion_orders':expected,
        'two_isogeny_targets':7,'scope':'realizations only, no completeness inferred from the source search',
        'primary_url':'https://www.math.fsu.edu/~hoeij/files/X1N/LowDegreePlaces',
        'reproduction_command':'magma -b realizations/verify_octic_examples.m',
        'files':{str(p.relative_to(HERE)):sha(p) for p in sorted(files)}}
(HERE/'data/realizations_manifest.json').write_text(json.dumps(report,indent=2)+'\n')
print('OCTIC_REALIZATIONS_AUDIT_DONE examples=10 fields_irreducible=true exact_full_torsions=true isogeny_targets=7 complete=true')
