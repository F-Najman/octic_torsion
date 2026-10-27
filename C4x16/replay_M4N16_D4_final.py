#!/usr/bin/env python3
"""Independent integer, coverage and parity replay of the octic M4N16 proof."""
import ast, json, runpy
from collections import Counter
from pathlib import Path
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form
P=Path(__file__).resolve().parent/'data'
def read(n):
 return ast.literal_eval((P/n).read_text().replace('[*','[').replace('*]',']'))
runpy.run_path(str(P.parent/'replay_M4N16_D4_diagonal.py'))
a=read('diagonal_4_16_D4_p5_p3_arithmetic.txt')
sigs,ds,inv5,r5,inv9,r9,cusp,source,target,matches,nc,reps,coverage,rat9,pc,counts,fps=a
# Matrix group with the correct simultaneous-sign quotient.
G={(pow(d,-1,4),b,c,d) for d in [1,3,5,7] for b in range(4) for c in range(4)}
assert len(G)==64
for A in G:
 for B in G:
  aa,b,c,d=A; e,f,g,h=B
  out=((e*aa)%4,(e*b+f*d)%4,(g*aa+h*c)%4,(4*g*b+h*d)%16)
  if out[3]>8:out=tuple((-x)%m for x,m in zip(out,[4,4,4,16]))
  assert out in G
full=read('M4N16_D4_Aut64_coverage_arithmetic.txt')
for row,old in zip(full,coverage):
 d,keys,labels,images,idx=row
 assert [d,keys,labels,idx]==old
 assert len(images)==len(idx)
 for j,orb in enumerate(images,1):
  assert len(orb)==64
  assert set(map(tuple,orb))=={tuple(k) for k,l in zip(keys,labels) if l==j}
 assert len(keys)=={3:320,4:1808}[d]
# Frobenius trace check is independent of the place enumeration.
def closed(fp,q,n,step=1):
 fp=list(reversed(fp)) # Magma Eltseq is in ascending powers.
 assert fp[0]==1
 S=[0]*(step*n+1)
 for k in range(1,len(S)):
  S[k]=-sum(fp[j]*S[k-j] for j in range(1,k) if j<len(fp))
  if k<len(fp):S[k]-=k*fp[k]
 out=[]
 for d in range(1,n+1):
  v=q**d+1-S[step*d]-sum(e*out[e-1] for e in range(1,d) if d%e==0)
  assert v%d==0;out.append(v//d)
 return out
assert closed(fps[0],5,4)==[16,4,64,164]
assert closed(fps[1],9,4,2)==[24,8,320,1808]
# The refined Hecke output uses EXACTLY the same support keys/partitions.
old=read('sieve_4_16_D4_F9_cusp_corrections_arithmetic.txt')
new=read('sieve_4_16_D4_F9_cusp_corrections_refined_arithmetic.txt')
assert new[:3]==old[:3]==coverage[0][1:]
assert old[4]==[2,4] and new[4]==[2] and new[5]==[4,3,2,1]
allowed=[0]+[1+8*v+w for v in range(8) for w in range(8) if w//2+1==new[5][v//2]]
assert len(allowed)==17
for (i,flags),(j,oldflags) in zip(new[3],old[3]):
 assert i==j and len(flags)==17 and len(oldflags)==65
 assert flags==[oldflags[k] for k in allowed]
assert [i for i,f in new[3] if any(f)]==[2]
last=[r for r in matches if target[r[0]-1][1] in new[4]]
assert last==[[29,[1055]]]
assert target[28][:3]==[3,2,5] and source[1054][:3]==[4,31,0]
# Replay the short global lift against the ACTUAL, potentially non-diagonal relations.
lift=read('M4N16_D4_doubled_cusp_lift.txt')
assert lift[:3]==[4,16,8] and lift[3:6]==[sigs,ds,1]
v=lift[6];assert lift[7:]==[1055,source[1054]] and len(v)==24
R=hermite_normal_form(Matrix([r[:-1] for r in r5]).T)
pred=[sum(v[j]*cusp[0][j][k] for j in range(24)) for k in range(len(inv5)-1)]
assert all(x.q==1 for x in R.inv()*Matrix([u-2*w for u,w in zip(pred,source[1054][3][:-1])]))
assert 8-sum(x*d for x,d in zip(v,ds))==9
rr=read('M4N16_D4_final_RR13_arithmetic.txt')
assert rr==[[13,1,[4,4],[1,1]]]
log=(P.parent/'logs'/'M4N16_D4_final_RR13.log').read_text()
for marker in ['M4N16_D4_FINAL_RR_FUNCTION_WITNESS','M4N16_D4_FINAL_EFFECTIVE_DIVISOR','M4N16_D4_DOUBLED_CUSP_FINITE_RR_DONE']:
 assert marker in log
assert 'Runtime error' not in log and 'Assertion failed' not in log
report={'degree':8,'excluded_group':[4,16],'relative_degree':4,'Aut_order':64,
 'NC_counts_at_5':[0,0,64,160],'NC_counts_at_9':[0,0,320,1808],
 'source_divisors':1184,'target_divisors':153,'diagonal_matches':17,
 'refined_Hecke_matches':last,'global_doubled_cusp_lift_replayed':True,
 'RR_prime':13,'RR_dimension':1,'effective_degrees':[4,4],'effective_multiplicities':[1,1],
 'even_multiplicity_obstruction':True,'remaining_candidates':0,
 'verification_boundary':'Independent Python checks actual integer relations, cusp graph, exhaustive support partitions, all64 orbit images, Frobenius counts, correction subsetting, exact doubled lift and parity. Exact native RR dimension/function/divisor are certified by the fresh Magma source; no independent Python RR implementation is claimed.'}
(P/'replay_M4N16_D4_final.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2));print('M4N16_D4_COMPLETE_EXCLUSION_REPLAY_DONE')
