#!/usr/bin/env python3
"""Independent Newton and finite orbit-coverage arithmetic replay."""
import ast,json,hashlib
from pathlib import Path
from collections import Counter
P=Path(__file__).resolve().parent
shas={}
def read(name):
 t=(P/'data'/name).read_text();shas[name]=hashlib.sha256(t.encode()).hexdigest()
 return ast.literal_eval(t.replace('[*','[').replace('*]',']'))
g,f,recorded=read('preflight_3_21_Frob_arithmetic.txt');assert g==25 and len(f)==51 and f[-1]==1
trace=[None]
for n in range(1,9):trace.append(-n*f[2*g-n]-sum(f[2*g-j]*trace[n-j] for j in range(1,n)))
closed=[None]
for n in range(1,5):
 z=4**n+1-trace[2*n]-sum(d*closed[d] for d in range(1,n) if n%d==0)
 assert z%n==0;closed.append(z//n)
assert closed[1:]==recorded==[24,0,8,144]
cusps=[24,0,8,0];nc=[a-b for a,b in zip(recorded,cusps)];assert nc==[0,0,0,144]
keys,labels,reps,tests,remaining=read('sieve_3_21_D4_arithmetic.txt')
assert len(keys)==len(set(map(tuple,keys)))==len(labels)==144
assert Counter(labels)=={1:72,2:72}
assert len(reps)==2 and all(1<=i<=144 for i in reps) and [labels[i-1] for i in reps]==[1,2]
assert tests==[[1,13,0],[2,13,0]] and not remaining
log=(P/'logs'/'sieve_3_21_D4_p2.log').read_text()
for marker in ['M3N21_D4_FULLAUT_REFERENCE_EQUAL 72','M3N21_D4_HECKE_REFERENCE_EQUAL 13 56','M3N21_D4_SIEVE_DONE tested=2 remaining=0 complete=true']:assert marker in log
report={'degree':8,'group':[3,21],'relative_degree':4,'genus':25,'finite_field':4,
 'closed_places_degrees1_to4':recorded,'noncuspidal_places_degrees1_to4':nc,
 'distinct_actual_place_keys':144,'full_automorphism_group_order':72,'orbits':[72,72],
 'Hecke_prime':13,'tests':tests,'survivors':remaining,'fresh_moduli_and_Hecke_reference_checks':True,
 'independent_Newton_and_coverage_replay':True,'input_sha256':shas}
(P/'data'/'replay_good_3_21_D4.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2));print('M3N21_GOOD_D4_NEWTON_COVERAGE_REPLAY_DONE')
