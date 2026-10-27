#!/usr/bin/env python3
"""Exact simultaneous source-Aut36 orbits, with repeated supports retained.

No Hecke propagation through the quotient is assumed. A global source
trace divisor can be replaced by its source automorphism image before
applying the necessary quotient condition to the selected representative.
"""
import ast,json,time
from collections import Counter
from pathlib import Path
import sys;sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'common'))
from finite_field5_checked import FF,EC,pack
p=Path(__file__).resolve().parent;t0=time.monotonic()
def rd(n):return ast.literal_eval((p/'data'/n).read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
def mul(E,Q,n):
 R=None
 while n:
  if n&1:R=E.add(R,Q)
  Q=E.add(Q,Q);n//=2
 return R
fields,rows=rd('M2N42_41_remaining_supports.m');Fs={d:FF(f) for d,f in fields}
pars=[(a,b,u) for a in range(3) for b in range(3) if a!=b for u in [1,2,4,5,8,10]];cache={}
def images(D,key):
 kk=D,tuple(key)
 if kk not in cache:
  F=Fs[D];b,c=pack(key[:D]),pack(key[D:]);E=EC(F,b,c);Q=(b,b);ps=[(0,0),mul(E,Q,21)];ps.append(E.add(*ps))
  cache[kk]=[F.key(E.bc(ps[a],E.add(ps[b],mul(E,Q,2*u)))) for a,b,u in pars]
 return cache[kk]
def base(sp):return tuple(sorted((D,pack(k[:D]),pack(k[D:])) for D,k in sp))
def divisor_images(sp):
 ims=[images(D,key) for D,key in sp]
 return [tuple(sorted((sp[j][0],*ims[j][i]) for j in range(len(sp)))) for i in range(36)]
groups={};seen=set();records=[]
for k,desc,sp in rows:
 assert k not in seen;seen.add(k);ims=divisor_images(sp);canonical=min(ims)
 if canonical not in groups:groups[canonical]={'representative_input_index':k,'representative_descriptor':desc,'representative_supports':sp,'representative_key':base(sp),'members':[]}
 g=groups[canonical];j=ims.index(g['representative_key']);assert ims[j]==g['representative_key']
 g['members'].append({'input_index':k,'descriptor':desc,'automorphism_parameter_to_representative':pars[j]})
assert len(seen)==41
for g in groups.values():g['representative_key']=[list(v) for v in g['representative_key']]
out={'source_Aut_order':36,'input_divisors':41,'orbit_count':len(groups),'group_sizes':[len(g['members']) for g in groups.values()],'orbits':list(groups.values()),'verification_boundary':__doc__,'seconds':round(time.monotonic()-t0,3)}
(p/'data'/'M2N42_41_full_Aut_orbits.json').write_text(json.dumps(out,indent=2)+'\n')
def mag(x):return '[* '+', '.join(map(mag,x))+' *]' if isinstance(x,list) else str(x)
reps=[[g['representative_input_index'],g['representative_descriptor'],g['representative_supports']] for g in groups.values()]
(p/'data'/'M2N42_full_Aut_representative_supports.m').write_text(mag([fields,reps])+'\n')
print(json.dumps({k:v for k,v in out.items() if k!='orbits'},indent=2));print([(g['representative_input_index'],g['representative_descriptor'],len(g['members'])) for g in groups.values()]);print('M2N42_41_FULL_AUT_ORBITS_DONE')
