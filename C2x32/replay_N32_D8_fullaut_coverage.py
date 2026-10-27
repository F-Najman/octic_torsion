#!/usr/bin/env python3
"""Independent Newton, exact full32-key partition and input-field replay."""
import ast,hashlib,json,importlib.util
from collections import Counter
from pathlib import Path
from math import gcd
from sympy import Poly,Symbol
P=Path(__file__).resolve().parent/'data'
shas={}
def read(name):
 t=(P/name).read_text();shas[name]=hashlib.sha256(t.encode()).hexdigest()
 return ast.literal_eval(t.replace('[*','[').replace('*]',']'))
g,f,closed=read('preflight_2_32_Frob_arithmetic.txt');assert g==37 and len(f)==75 and f[-1]==1
tr=[None];cc=[None]
for n in range(1,9):
 tr.append(-n*f[74-n]-sum(f[74-j]*tr[n-j] for j in range(1,n)))
 a=3**n+1-tr[n]-sum(d*cc[d] for d in range(1,n) if n%d==0)
 assert a%n==0;cc.append(a//n)
assert cc[1:]==closed==[16,8,0,2,0,192,320,962]
records,coverage=read('N32_D8_fullaut_Hecke_input.txt');assert len(records)==50
assert Counter(r[0] for r in records)=={6:8,7:10,8:32}
D=[1,3,5,7,9,11,13,15];pars={(b,c,d) for b in (0,1) for c in (0,1) for d in D}
for b,c,d in pars:
 for B,C,E in pars:
  dn=(E*d+16*C*b)%32
  if dn not in D:dn=(-dn)%32
  assert ((b+B)%2,(c+C)%2,dn) in pars
x=Symbol('x');covermaps={}
for dd,keys,labels,orbs,reps in coverage:
 expected={6:(192,8),7:(320,10),8:(960,32)}[dd]
 assert (len(keys),len(orbs))==expected and len(labels)==len(keys) and len(reps)==len(orbs)
 allkeys=set(map(tuple,keys));assert len(allkeys)==len(keys)
 assert all(len(k)==2*dd and all(a in (0,1,2) for a in k) for k in keys)
 seen={}
 for j,orb in enumerate(orbs,1):
  assert len(orb)==32
  ss=set(map(tuple,orb));assert ss<=allkeys
  for k in ss:assert k not in seen;seen[k]=j
 assert set(seen)==allkeys
 assert [seen[tuple(k)] for k in keys]==labels
 assert [labels[i-1] for i in reps]==list(range(1,len(orbs)+1))
 covermaps[dd]=(keys,labels,reps)
for dd,i,ff,vs in records:
 assert len(ff)==dd+1 and ff[-1]==1 and all(a in (0,1,2) for a in ff)
 assert Poly.from_list(list(reversed(ff)),x,modulus=3).is_irreducible
 assert len(vs)==2 and all(len(v)==dd for v in vs)
 def mul(a,b):
  z=[0]*(2*dd-1)
  for k,u in enumerate(a):
   for l,v in enumerate(b):z[k+l]=(z[k+l]+u*v)%3
  for k in range(2*dd-2,dd-1,-1):
   for l in range(dd):z[k-dd+l]=(z[k-dd+l]-z[k]*ff[l])%3
  return z[:dd]
 vals=[list(v) for v in vs];ks=[]
 for _ in range(dd):
  ks.append(tuple(vals[0]+vals[1]));vals=[mul(mul(v,v),v) for v in vals]
 assert len(set(ks))==dd and vals==vs
 keys,labels,reps=covermaps[dd];assert min(ks)==tuple(keys[reps[i-1]-1])
spec=importlib.util.spec_from_file_location('tate',P.parents[1]/'common'/'cusp_conditions_relative.py');t=importlib.util.module_from_spec(spec);spec.loader.exec_module(t)
a=t.compute(2,32,3,8);gal=[u for u in range(1,32) if gcd(u,32)==1]
assert a['allowed_residues_by_complement_degree'][0]==gal
assert 5 in a['allowed_residues_by_complement_degree'][0]
assert 11 in a['allowed_residues_by_complement_degree'][1]
assert 19 in a['allowed_residues_by_complement_degree'][1]
assert Counter(a['closed_cusp_degrees'])=={1:16,2:8,4:2,8:2}
report={'degree':8,'group':[2,32],'finite_field':3,'genus':37,'full_aut_group_order':32,'NC_counts':[0,0,0,0,0,192,320,960],
 'representatives':[8,10,32],'total_representatives':50,'all_actual_keys':1472,'all_finite_field_record_roundtrips_replayed':True,
 'independent_Newton_full_partition_and_parameter_group_replay':True,'q11_q19_admissible_degree2_complements':True,
 'all_Galois_stabilizes_degree1_cusp_orbits':True,'input_sha256':shas}
(P/'replay_N32_D8_fullaut_coverage.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2));print('N32_D8_FULL32_AUT_COVERAGE_REPLAY_DONE')
