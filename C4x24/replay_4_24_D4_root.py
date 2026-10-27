#!/usr/bin/env python3
from pathlib import Path
import ast,json
from collections import Counter
P=Path(__file__).resolve().parent/'data'
def read(s):return ast.literal_eval((P/s).read_text().replace('[*','[').replace('*]',']'))
ex,rem=read('sieve_4_24_D4_fullaut_root_arithmetic.txt');assert rem==[] and len(ex)==1
d,keys,labels,rep,tests=ex[0];assert d==4 and len(keys)==len(labels)==96 and len({tuple(k) for k in keys})==96
assert set(labels)=={1,2,3,4} and rep==[labels.index(i)+1 for i in range(1,5)]
assert tests==[[i,13,0] for i in range(1,5)]
f,cl=read('control_4_24_D4_root_Frobenius.txt');assert len(f)==67
tr=[0];cc=[]
for n in range(1,5):
 tr.append(-n*f[66-n]-sum(f[66-j]*tr[n-j] for j in range(1,n)))
 v=5**n+1-tr[n]-sum(e*cc[e-1] for e in range(1,n) if n%e==0);assert v%n==0;cc.append(v//n)
assert cl==cc==[16,24,0,96]
G={(pow(d,-1,4)*(1+6*b*c)%4,b,c,d) for b in range(4) for c in range(4) for d in [1,5,7,11]};assert len(G)==64
def mul(x,y):
 a,b,c,d=x;aa,bb,ccc,dd=y
 r=((a*aa+6*b*ccc)%4,(a*bb+b*dd)%4,(c*aa+d*ccc)%4,(6*c*bb+d*dd)%24)
 if r[3] not in [1,5,7,11]:r=tuple((-t)%m for t,m in zip(r,[4,4,4,24]))
 return r
assert all(mul(x,y) in G for x in G for y in G)
for name,marker in [('sieve_4_24_D4_fullaut_root.log','M4N24_ROOT_D4_FULLAUT_HECKE_DONE'),('control_4_24_D4_root.log','M4N24_ROOT_CONTROL_NEWTON_COUNTS')]:
 log=(P.parent/'logs'/name).read_text();assert marker in log and 'Runtime error' not in log and 'Assertion failed' not in log
report={'degree':8,'group':[4,24],'prime':5,'genus':33,'NC_counts':[0,0,0,96],'full_aut_order':64,'full_aut_orbit_sizes':dict(Counter(labels)),'representatives':4,'all_q13_nonprincipal':True,'independent_Newton_counts':True,'group_closure_verified':True,'native_generic_Hecke_reference':False,'additional_generic_reference_status':'900second bounded run incomplete; not used','exact_fast_kernel_orbit_assertions':True,'complete':True}
(P/'replay_4_24_D4_root.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2));print('M4N24_D4_INDEPENDENT_REPLAY_COMPLETE')
