#!/usr/bin/env python3
"""Bind the octic85 quotient-descent and99 local-diamond exclusions."""
import ast
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import re

H=Path(__file__).resolve().parent
cache=json.loads((H.parent/'common'/'data'/'lmfdb_newforms_cache.json').read_text())
cmtext=(H.parent/'common'/'data'/'dcm_list_all_min_orders_X1_10k.m').read_text()
cm={int(n):int(d) for n,d,k in re.findall(r'\[\*\s*(\d+),\s*(\d+),\s*(\d+),',cmtext)}


def divs(n):return [m for m in range(1,n+1) if n%m==0]
def phi(n):return sum(math.gcd(a,n)==1 for a in range(1,n+1))
def order(n,a):
    p,k=a%n,1
    while p not in (1,n-1):p,k=p*a%n,k+1
    return k


total,positive=0,[]
for m in divs(85):
    for f in cache[str(m)]:
        assert f['analytic_rank_proved'] is True and f['analytic_rank'] is not None
        total+=len(divs(85//m))*f['dim']
        if f['analytic_rank']>0:positive.append({'level':m,**f})
assert total==225 and {f['char_orbit_label'] for f in positive if f['char_order']>1}=={'j'}
chlog=(H/'logs'/'check_85_characters.log').read_text()
assert chlog.rstrip().endswith('OCTIC85_CHARACTER_CERTIFICATE_DONE')
assert 'value_at2=1 verified=true' in chlog
assert order(85,2)==8 and cm[85]==32 and Fraction(325*3456,32768)>32
ff=(H/'logs'/'finite_fields_2_octic.log').read_text()
assert ff.rstrip().endswith('DONE')
assert re.search(r'TORSION n=85 p=2 : k with Z/85 over F_p\^k = \[\s*\]',ff)
# The check that the copy of arXiv:2001.05616v3 states Theorem 4.6 (rational points on X0(85)) was removed:
# the paper is not part of this repository.
data85={'degree':8,'excluded_cyclic_order':85,'status':'complete','a':2,'diamond_order':8,
        'rank_dimension_sum':225,'positive_rank_forms':positive,'minimum_CM_degree':32,
        'gonality_lower_bound':str(Fraction(325*3456,32768)),
        'principal_function_degree_bound':32,'orbit_size':8,'quotient_image_degree':1,
        'final_obstruction':'X0(85)(Q) has no noncuspidal point',
        'reference':{'url':'https://arxiv.org/pdf/2001.05616v3','theorem':'4.6','PDF_page':8}}

s=(H/'logs'/'winding99_octic.log').read_text()
assert s.rstrip().endswith('WINDING_CHARACTER_PROBE_DONE')
assert not re.search(r'Runtime error|Assertion failed|SystemError|User error',s)
rows=re.findall(r'WINDING N=99 order=(\d+) conductor=(\d+) relative_dimension=(\d+) winding_dimension=(\d+) coefficient_degree=(\d+) time=',s)
assert sum(int(d)*int(c) for o,cond,d,w,c in rows)==281
assert all(o=='1' or d==w for o,cond,d,w,c in rows)
assert 'bad_orders=[ 1 ]' in s and order(99,2)==15 and cm[99]>8
ff=(H/'logs'/'finite_fields_99_octic.log').read_text();assert ff.rstrip().endswith('DONE')
assert re.search(r'TORSION n=99 p=2 : k with Z/99 over F_p\^k = \[\s*\]',ff)
hasse=[]
for e in range(1,9):
    q=2**e;t=math.isqrt(4*q)
    multiples=[m for m in range(99,q+t+2,99) if abs(m-q-1)<=t]
    assert multiples==[]
    hasse.append({'e':e,'q':q,'maximum_trace_absolute_value':t,'multiples99_in_interval':multiples})
assert Fraction(325*4320,32768)>32
data99={'degree':8,'excluded_cyclic_order':99,'status':'complete','a':2,'diamond_order':15,
        'rank_dimension_sum':281,'minimum_CM_degree':cm[99],
        'gonality_lower_bound':str(Fraction(325*4320,32768)),
        'principal_function_degree_bound':32,'good_reduction_degrees_through8':[],
        'independent_Hasse_check':hasse}

common=['certify_85_99.py','inputs/dcm_list_all_min_orders_X1_10k.m']
files85=common+['N85_exclusion.tex','inputs/lmfdb_newforms_cache.json',
    'check_85_characters.m','check_85_characters.log','finite_fields_2_octic.m',
    'finite_fields_2_octic.log','extra_inputs/isogeny_graphs_2001.05616v3.pdf',
    'extra_inputs/isogeny_graphs_2001.05616v3.txt','../shared/mdmagma/Magma/magma.spec']
files99=common+['N99_exclusion.tex','winding_probe.m','winding99_octic.log',
    'finite_fields_2_octic.m','finite_fields_99_octic.log']
for data,files in [(data85,files85),(data99,files99)]:
    # The SHA-256 binding of the files of the original record (data['files']) was removed.
    n=data['excluded_cyclic_order'];(H/'data'/f'N{n}_manifest.json').write_text(json.dumps(data,indent=2)+'\n')
    print('OCTIC_ADDITIONAL_EXCLUSION_DONE',n)
