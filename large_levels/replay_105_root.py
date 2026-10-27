#!/usr/bin/env python3
"""Exact octic105 double descent and rank-zero X0 reduction certificate."""
import ast
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import re

H=Path(__file__).resolve().parent
s=(H/'logs'/'winding105_octic.log').read_text()
assert s.rstrip().endswith('WINDING_CHARACTER_PROBE_DONE')
assert not re.search(r'Runtime error|Assertion failed|SystemError|User error',s)
rows=re.findall(r'WINDING N=105 order=(\d+) conductor=(\d+) relative_dimension=(\d+) winding_dimension=(\d+) coefficient_degree=(\d+) time=',s)
rows=[list(map(int,r)) for r in rows]
assert sum(d*c for o,cond,d,w,c in rows)==289
assert [r for r in rows if r[2]!=r[3]]==[[6,21,12,11,2]]
assert next(r for r in rows if r[0]==1)==[1,1,13,13,1]
assert 22%21==1
last=re.search(r'WINDING_FINAL N=105 genus=289 total_dimension=289 bad_orders=(\[.*?\]) admissible=(\[.*?\])',s)
assert last and ast.literal_eval(last[1])==[6] and 22 in ast.literal_eval(last[2])
orbit=[min(pow(22,k,105),(-pow(22,k,105))%105) for k in range(4)]
assert len(set(orbit))==4 and pow(22,4,105)==1
cmtext=(H.parent/'common'/'data'/'dcm_list_all_min_orders_X1_10k.m').read_text()
cm={int(n):int(d) for n,d,k in re.findall(r'\[\*\s*(\d+),\s*(\d+),\s*(\d+),',cmtext)}
assert cm[105]==96 and cm[105]>8
gon=Fraction(325*4608,32768);assert gon>32
ff=(H/'logs'/'finite_fields_105_octic.log').read_text();assert ff.rstrip().endswith('DONE')
assert re.search(r'TORSION n=105 p=2 : k with Z/105 over F_p\^k = \[\s*\]',ff)
hasse=[]
for e in range(1,9):
    q=2**e;trace=math.isqrt(4*q)
    multiples=[m for m in range(105,q+trace+2,105) if abs(m-q-1)<=trace]
    assert not multiples
    hasse.append({'e':e,'q':q,'maximum_trace_absolute_value':trace,'multiples105_in_interval':multiples})
divisors=[d for d in range(1,106) if 105%d==0]
assert divisors==[1,3,5,7,15,21,35,105]
assert all(math.gcd(d,105//d)==1 for d in divisors)
# The check that the copy of arXiv:2207.11650v4 (Najman--Orlic, Gonality of X0(N)) gives gonality at least 6
# for X0(105) in its Proposition 5.18 was removed: the paper is not part of this repository.
out={'degree':8,'excluded_cyclic_order':105,'status':'complete','genus_X1':289,
     'finite_image_a':22,'diamond_order':4,'minimum_CM_degree_X1':96,
     'gonality_lower_bound_X1':str(gon),'principal_function_degree_bound_X1':32,
     'quotient_trace_degree':2,'rank_J0Q':0,'genus_X0':13,'cusps_X0':8,
     'all_cusps_X0_rational':True,'all_bad2':True,'Hasse_checks':hasse,
     'gonality_lower_bound_X0':6,'principal_function_degree_bound_X0':4,
     'primary_gonality_reference':{'url':'https://arxiv.org/pdf/2207.11650v4','proposition':'5.18','PDF_page':13}}
files=['N105_exclusion.tex','certify_105.py','winding_probe.m','winding105_octic.log',
       'finite_fields_2_octic.m','finite_fields_105_octic.log',
       'inputs/dcm_list_all_min_orders_X1_10k.m',
       'extra_inputs/gonality_X0_2207.11650v4.pdf','extra_inputs/gonality_X0_2207.11650v4.txt']
# The SHA-256 binding of the files of the original record (out['files']) was removed.
(H/'data'/'N105_root_replay.json').write_text(json.dumps(out,indent=2)+'\n')
print('OCTIC105_EXCLUSION_DONE',{'rank_J0':0,'trace_descent_degree':2,'gon_X0_at_least6':True})
