#!/usr/bin/env python3
"""Assemble all exact checks of the octic C2x28 exclusion."""
import ast,hashlib,json,runpy
from pathlib import Path
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form
P=Path(__file__).resolve().parent/'data'
for script in ['replay_N28_D8_Hecke.py','replay_N28_D8_quotients.py','replay_N28_D8_exact_point.py']:
 runpy.run_path(str(P.parent/script))
def read(name):
 return ast.literal_eval((P/name).read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
sigs,ds,base,lifts,inv,rels,cusp,tests=read('N28_D8_quotient_cusp_lifts_certify_arithmetic.txt')
old=read('N28_D8_quotient_cusp_RR_arithmetic.txt')
assert [sigs,ds,base,lifts]==old[:4]
assert base==1 and sigs[0]==[-24,1] and ds[0]==1
assert lifts[0]==[0,0,1,0,0,-1,-1,-1,0,0,0,0,0,-1,-1]
assert 8-sum(v*d for v,d in zip(lifts[0],ds))==22
assert all(row[-1]==0 for row in rels+cusp)
R=hermite_normal_form(Matrix([r[:-1] for r in rels]).T);Ri=R.inv()
H=hermite_normal_form(Matrix([r[:-1] for r in rels+cusp]).T)
assert abs(R.det())==359424 and abs(H.det())==4
assert len(tests)==4
for kind,e,mult,cl in tests:
 assert e*mult==8 and cl[-1]==0
 predicted=[sum(lifts[kind][j]*cusp[j][k] for j in range(15)) for k in range(4)]
 assert all(x.q==1 for x in Ri*Matrix([x-y for x,y in zip(predicted,cl[:-1])]))
log=(P.parent/'logs'/'N28_D8_quotient_charzero_certify_v2.log').read_text()
for marker in ['N28_D8_QUOTIENT_CHARZERO_RR_DIMENSION 1',
 'N28_D8_QUOTIENT_CHARZERO_EXACT_PLACE_MATCH true',
 'N28_D8_QUOTIENT_CHARZERO_RR_DONE',
 'N28_D8_QUOTIENT_PSI2_FACTOR_DEGREES [ 1, 2 ]',
 'N28_D8_QUOTIENT_DIRECT_TORSION_DONE']:
 assert marker in log
assert 'Runtime error' not in log and 'Assertion failed' not in log
assert (P/'N28_D8_quotient_direct_torsion_arithmetic.txt').read_bytes()==(P/'N28_D8_quotient_charzero_certify_arithmetic.txt').read_bytes()
hh=json.loads((P/'replay_N28_D8_Hecke.json').read_text())
qq=json.loads((P/'replay_N28_D8_quotients.json').read_text())
assert hh['unresolved_representatives']==[35,46] and qq['remaining_full_Aut24_orbits']==[35]
report={'degree':8,'excluded_groups':[[2,28],[1,56]],'status':'complete; awaiting independent root review',
 'full_aut_representatives':49,'two_prime_Hecke_survivors':[35,46],
 'quotient_class_obstruction_kills':[46],'exact_global_cusp_lift_replayed':True,
 'final_orbit':35,'global_RR_dimension':1,'unique_effective_divisor_degree':8,
 'psi2_factor_degrees':[1,2],'independent_good_prime_2torsion_obstruction':23,
 'final_source_orbits_remaining':0,
 'verification_boundary':'Python independently replays field/count/partition/integer-class/local2 arithmetic. Native exact RR dimension and divisor identification are checked by the completed fresh Magma end-to-end source; no independent Python RR implementation is claimed.'}
(P/'replay_N28_D8_final.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2));print('N28_D8_COMPLETE_EXCLUSION_REPLAY_DONE')
