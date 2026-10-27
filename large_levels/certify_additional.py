#!/usr/bin/env python3
"""Audit and bind the separate completed octic exclusions81 and175."""
import ast
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import re

H = Path(__file__).resolve().parent


def winding(n):
    text=(H/'logs'/f'winding{n}_octic.log').read_text()
    assert text.rstrip().endswith('WINDING_CHARACTER_PROBE_DONE')
    assert not re.search(r'Runtime error|Assertion failed|SystemError|User error',text)
    rows=re.findall(r'WINDING N=(\d+) order=(\d+) conductor=(\d+) relative_dimension=(\d+) winding_dimension=(\d+) coefficient_degree=(\d+) time=([0-9.]+)',text)
    data=[list(map(int,row[:6])) for row in rows]
    final=re.findall(r'WINDING_FINAL N=(\d+) genus=(\d+) total_dimension=(\d+) bad_orders=(\[.*?\]) admissible=(\[.*?\])',text)
    assert len(final)==1
    nn,g,total=map(int,final[0][:3]);assert nn==n and g==total
    assert sum(d*c for nn,o,cond,d,w,c in data)==g
    return g,ast.literal_eval(final[0][3]),ast.literal_eval(final[0][4]),data


g,bad,adm,rows=winding(81)
assert g==190 and bad==[] and all(d==w for nn,o,cond,d,w,c in rows)
replay=ast.literal_eval((H/'logs'/'two_torsion_81_replay.log').read_text())
assert replay['level']==81 and replay['integral_homology_rank']==380
assert replay['kernel_dimension']==0 and replay['independent_binary_replay']
log=(H/'logs'/'rational_two_torsion_81_v2.log').read_text()
assert 'RATIONAL_TWO_DONE N=81 kernel_dimension=0 diamond_trivial=true' in log
assert 'q=5 kernel_dimension=8 ' in log and 'q=7 kernel_dimension=0 ' in log
assert not re.search(r'Runtime error|Assertion failed|SystemError|User error|System Error',log)
cusps=next(x for x in json.loads((H/'data'/'cusp_lift_checks.json').read_text()) if x['N']==81)
assert all(c['global_degree']==c['residue_degree'] for c in cusps['F2_closed_cusps'] if c['residue_degree']<=8)
ff=(H/'logs'/'finite_fields_2_octic.log').read_text();assert ff.rstrip().endswith('DONE')
assert re.search(r'TORSION n=81 p=2 : k with Z/81 over F_p\^k = \[\s*\]',ff)
assert Fraction(325*2916,32768)>8
data81={'degree':8,'excluded_cyclic_order':81,'status':'complete','genus':190,
        'rankQ':0,'rational2torsion':0,'integral_homology_rank':380,
        'upper2torsion_kernel_dimensions':[8,0],'good_reduction_degrees_through8':[],
        'all_cuspidal_degree8_divisors_liftQ':True,'gonality_lower_bound':str(Fraction(325*2916,32768)),
        'principal_function_degree_bound':8,'CM_exclusion_needed':False}

g,bad,adm,rows=winding(175)
assert g==1033 and set(bad)=={1,5} and 24 in adm
z=24%175;order=1
while z not in (1,174):z=z*24%175;order+=1
assert order==3 and 8%order
cmtext=(H.parent/'common'/'data'/'dcm_list_all_min_orders_X1_10k.m').read_text()
cm={int(n):int(d) for n,d,k in re.findall(r'\[\*\s*(\d+),\s*(\d+),\s*(\d+),',cmtext)}
assert cm[175]==120 and cm[175]>8
assert Fraction(325*14400,32768)>64
data175={'degree':8,'excluded_cyclic_order':175,'status':'complete','genus':1033,
         'a':24,'diamond_order':3,'incomplete_winding_character_orders':bad,
         'q':3,'kq':8,'index':14400,'gonality_lower_bound':str(Fraction(325*14400,32768)),
         'principal_function_degree_bound':64,'minimum_CM_degree':120}

files81=['N81_exclusion.tex','certify_additional.py','winding_probe.m','winding81_octic.log',
         'rational_two_torsion_8gb.m','rational_two_torsion_81_v2.log',
         'two_torsion_81_v2_matrices.m','replay_two_torsion_81.py','two_torsion_81_replay.log',
         'cusp_lift_checks.py','cusp_lift_checks.json','cusp_lift_checks.log',
         'finite_fields_2_octic.m','finite_fields_2_octic.log']
files175=['N175_exclusion.tex','certify_additional.py','winding_probe.m','winding175_octic.log',
          'inputs/dcm_list_all_min_orders_X1_10k.m','inputs/least_cm_degree_README.md',
          'inputs/DN_2412.16016_main.tex']
for data,files in [(data81,files81),(data175,files175)]:
    # The SHA-256 binding of the files of the original record (data['files']) was removed.
    n=data['excluded_cyclic_order']
    (H/'data'/f'N{n}_manifest.json').write_text(json.dumps(data,indent=2)+'\n')
    print('OCTIC_ADDITIONAL_EXCLUSION_DONE',n)
