#!/usr/bin/env python3
"""Exact octic91 rank, automorphism and mixed-difference certificate."""
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import re

H=Path(__file__).resolve().parent
cache=json.loads((H.parent/'common'/'data'/'lmfdb_newforms_cache.json').read_text())
total,positive=0,[]
for m in [1,7,13,91]:
    mult=sum((91//m)%d==0 for d in range(1,91//m+1))
    for f in cache[str(m)]:
        assert f['analytic_rank_proved'] is True and f['analytic_rank'] is not None
        total+=mult*f['dim']
        if f['analytic_rank']>0:positive.append({'level':m,'multiplicity':mult,**f})
assert total==265
assert {f['char_orbit_label'] for f in positive if f['char_order']>1}=={'g','u'}
char=(H/'logs'/'check_91_characters.log').read_text()
assert char.rstrip().endswith('OCTIC91_CHARACTER_CERTIFICATE_DONE')
assert 'orbit=g ' in char and 'orbit=u ' in char
assert char.count('value_at3=1 verified=true')==2

def canonical(a):return min(a%91,(-a)%91)
group={(i,j):canonical(pow(16,i,91)*pow(3,j,91)) for i in range(3) for j in range(6)}
assert len(set(group.values()))==18 and canonical(pow(16,3,91))==1
assert canonical(pow(3,6,91))==1
# Integer permutation matrices, row-vector action on the regular orbit.
S=[[0]*18 for _ in range(18)];T=[[0]*18 for _ in range(18)]
for i in range(3):
    for j in range(6):
        a=6*i+j;S[a][6*((i+1)%3)+j]=1;T[a][6*i+(j+1)%6]=1
I=[[int(i==j) for j in range(18)] for i in range(18)]
A=[[S[i][j]-I[i][j] for j in range(18)] for i in range(18)]
B=[[T[i][j]-I[i][j] for j in range(18)] for i in range(18)]
mixed=[[sum(A[i][k]*B[k][j] for k in range(18)) for j in range(18)] for i in range(18)]
witness=[(i*j)%3 for i in range(3) for j in range(6)]
assert all(sum(row[j]*witness[j] for j in range(18))%3==1 for row in mixed)
# If D*mixed=0, this equality forces sum(D)=0 modulo3.
group_ring={'group':'C3 x C6','elements':[[i,j,a] for (i,j),a in group.items()],
            'mixed_difference_matrix':mixed,'augmentation_witness_mod3':witness,
            'mixed_times_witness_mod3':[1]*18,'verified':True}
(H/'data'/'N91_group_ring_certificate.json').write_text(json.dumps(group_ring,indent=2)+'\n')

auto=(H/'logs'/'check_good8_91_certified.log').read_text()
assert auto.rstrip().endswith('OCTIC_GOOD91_AUTOMORPHISM_DONE complete=true')
assert 'qualifying_twists=2' in auto
assert auto.count('stabilizer_multipliers=[ 1, 16, 17, 74, 75, 90 ]')==2
assert not re.search(r'Runtime error|Assertion failed|System Error|SystemError|User error',auto)
ff=(H/'logs'/'finite_fields_2_octic.log').read_text();assert ff.rstrip().endswith('DONE')
assert re.search(r'TORSION n=91 p=2 : k with Z/91 over F_p\^k = \[ 8 \]',ff)
cmtext=(H.parent/'common'/'data'/'dcm_list_all_min_orders_X1_10k.m').read_text()
cm={int(n):int(d) for n,d,k in re.findall(r'\[\*\s*(\d+),\s*(\d+),\s*(\d+),',cmtext)}
assert cm[91]==24 and cm[91]>8
gon=Fraction(325*4032,32768);assert gon>32
out={'degree':8,'excluded_cyclic_order':91,'status':'complete','genus':265,
     'rank_dimension_sum':265,'finite_image_diamond':3,'finite_image_diamond_order':6,
     'residue_automorphism_diamond':16,'residue_automorphism_diamond_order':3,
     'combined_diamond_group_order':18,'minimum_CM_degree':24,
     'gonality_lower_bound':str(gon),'principal_function_degree_bound':32,
     'good_reduction_degrees_through8':[8],
     'branches':{'entirely_cuspidal':'usual finite-image diamond argument',
                 'pure_good8':'mixed diamond relation and independent mod3 augmentation witness'},
     'positive_rank_forms':positive}
files=['N91_exclusion.tex','certify_91.py','N91_group_ring_certificate.json',
       'check_91_characters.m','check_91_characters.log',
       'check_good8_91_certified.m','check_good8_91_certified.log',
       'finite_fields_2_octic.m','finite_fields_2_octic.log',
       'inputs/lmfdb_newforms_cache.json','inputs/dcm_list_all_min_orders_X1_10k.m',
       '../shared/mdmagma/Magma/magma.spec']
# The SHA-256 binding of the files of the original record (out['files']) was removed.
(H/'data'/'N91_manifest.json').write_text(json.dumps(out,indent=2)+'\n')
print('OCTIC91_EXCLUSION_DONE',{'group_order':18,'mod3_witness':True,'both_branches':True})
