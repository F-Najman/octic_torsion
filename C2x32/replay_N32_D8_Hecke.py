#!/usr/bin/env python3
"""Independently assemble every full-Aut32 record and exact Hecke decision."""
import ast,hashlib,json,re,runpy
from collections import Counter
from pathlib import Path
P=Path(__file__).resolve().parent/'data'
runpy.run_path(str(P.parent/'replay_N32_D8_fullaut_coverage.py'))
bound={}
def content(name):
 raw=(P.parent/'logs'/name if name.endswith('.log') else P/name).read_bytes();bound[name]=hashlib.sha256(raw).hexdigest();return raw.decode()
def read(name):return ast.literal_eval(content(name).replace('[*','[').replace('*]',']'))
records,coverage=read('N32_D8_fullaut_Hecke_input.txt')
decisions={}
for shard in [0,1]:
 text=content(f'N32_D8_fullaut_Hecke_shard{shard}.log')
 assert 'Runtime error' not in text and 'Assertion failed' not in text
 seen={}
 for s,n,d,i,q,flag in re.findall(r'N32_D8_FULLAUT_HIGH_HECKE_CASE (\d+) (\d+) (\d+) (\d+) (\d+) (true|false)',text):
  s,n,d,i,q=map(int,[s,n,d,i,q]);assert s==shard and (n-1)%2==shard
  assert records[n-1][:2]==[d,i]
  seen.setdefault(n,[]).append((q,flag=='true'))
 expected={n for n in range(1,51) if (n-1)%2==shard}
 assert set(seen)==expected
 for n,flags in seen.items():
  qs=[11,19] if records[n-1][0]==6 else [5,11]
  assert flags in [[(qs[0],False)],[(qs[0],True),(qs[1],False)],[(qs[0],True),(qs[1],True)]]
 survivors=[n for n,flags in seen.items() if flags[-1][1]]
 m=re.search(r'N32_D8_FULLAUT_HIGH_HECKE_SURVIVORS \d+ (\[[^\]]*\])',text)
 assert m and ast.literal_eval(m[1])==survivors
 assert re.search(rf'N32_D8_FULLAUT_HIGH_HECKE_DONE shard={shard} tested=25 remaining={len(survivors)} complete=true',text)
 decisions.update(seen)
pending=[n for n,flags in sorted(decisions.items()) if flags[-1][1]]
assert pending==[2] and records[1][:2]==[6,2]
sigs,pairs,vs,tests,remaining=read('N32_D8_degree6_correction_v2_arithmetic.txt')
assert len(sigs)==8 and all(len(s)==3 for s in sigs)
assert sorted(Counter(map(tuple,sigs)).values())==[1,1,1,1,2,2]
wantpairs=[[i+1 for i,s in enumerate(sigs) if tuple(s)==sig] for sig in sorted(set(map(tuple,sigs))) if sum(tuple(s)==sig for s in sigs)==2]
assert pairs==wantpairs and sorted(map(len,pairs))==[2,2]
wantvs=[]
for pair in pairs:
 for sign in [-1,1]:
  v=[0]*8;v[pair[0]-1]=-sign;v[pair[1]-1]=sign;wantvs.append(v)
assert vs==wantvs and len(tests)==1 and remaining==[]
n,q,corr,flags=tests[0];assert [n,q]==[2,5]
want=[[0]*8]+[[v[i]+5*w[i] for i in range(8)] for v in vs for w in vs]
assert corr==want and len(set(map(tuple,corr)))==17 and flags==[0]*17
text=content('N32_D8_degree6_cusp_corrections_v2.log')
assert 'N32_D8_FRESH_T5_GENERIC_REFERENCE_EQUAL 36' in text
assert re.search(r'N32_D8_DEGREE6_CORRECTIONS_DONE \[\s*\]',text)
assert 'Runtime error' not in text and 'Assertion failed' not in text
report={'degree':8,'excluded_groups':[[2,32],[1,64]],'genus':37,
 'NC_counts':[0,0,0,0,0,192,320,960],'full_Aut_order':32,
 'support_records':50,'records_by_degree':[8,10,32],
 'all_shard_decisions_replayed':True,'two_prime_survivors':[2],
 'degree6_cusp_unit_profile_sizes':[2,1,2,1,1,1],
 'final_q5_corrections':17,'generic_q5_reference_degree':36,
 'remaining_candidates':0,'status':'complete; awaiting independent root review',
 'verification_boundary':'Independent Python checks Frobenius/place counts, full orbit partitions, input field arithmetic, Tate complement admissibility and complete assembly of exact Hecke/principality decisions. Divisor arithmetic and IsPrincipal are the exact native Magma computations; no independent Python divisor-arithmetic implementation is claimed.',
 'checksums':bound}
(P/'replay_N32_D8_Hecke.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2));print('N32_D8_COMPLETE_EXCLUSION_REPLAY_DONE')
