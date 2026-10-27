#!/usr/bin/env python3
"""Independent exhaustive support-index and Newton coverage for the 60 pairs."""
import ast,hashlib,json,runpy
from collections import Counter
from math import comb,gcd
from functools import reduce
from pathlib import Path
P=Path(__file__).resolve().parent/'data'
def read(name):return ast.literal_eval((P/name).read_text().replace('[*','[').replace('*]',']'))
frobs=read('M3N15_necessity_frobenius_arithmetic.txt');polys={p:f for p,f,o in frobs}
orders={p:o for p,f,o in frobs}
def closed(p,base_degree):
 f=polys[p];n=len(f)-1;s=[0]
 for k in range(1,4*base_degree+1):s.append(-k*f[n-k]-sum(f[n-j]*s[k-j] for j in range(1,k)))
 cc=[]
 for d in range(1,5):
  z=p**(base_degree*d)+1-s[base_degree*d]-sum(e*cc[e-1] for e in range(1,d) if d%e==0)
  assert z%d==0;cc.append(z//d)
 return cc
assert closed(2,2)==[16,8,0,102] and closed(7,1)==[16,24,144,400]
h=1990656;assert reduce(gcd,[orders[p] for p in [7,13,19]])==3*h
assert orders=={2:5971968,7:149299200,13:12230590464,19:247669456896}
a=read('diagonal_3_15_D4_arithmetic.txt');s4,s7,pairs=a[12:15]
assert len(s4)==5102 and len(s7)==6264 and len(pairs)==60
cuspal=comb(19,4)+8*comb(17,2)+comb(9,2);assert cuspal==5000
assert Counter((typ,i) for typ,i,cl,hcl,pred in s4)==Counter({**{(0,i):1 for i in range(1,5001)},**{(4,i):1 for i in range(1,103)}})
expected=Counter()
for i in range(1,25):
 for j in range(1,137):expected[(2,i,j)]+=1
for i in range(1,145):
 for j in range(1,17):expected[(3,i,j)]+=1
for i in range(1,397):expected[(4,i,1)]+=1
for i in range(1,25):
 for j in range(i,25):expected[(22,i,j)]+=1
assert Counter((typ,i,j) for typ,i,j,cl in s7)==expected
runpy.run_path(str(P.parent/'replay_diagonal_3_15_D4.py'))
runpy.run_path(str(P.parent/'replay_exact_3_15_D4.py'))
report={'group':[3,15],'relative_degree':4,'F4_all_cuspal_divisors':5000,'F4_NC_quartics':102,
 'F7_single2_plus_rational_cusp_pair':24*136,'F7_single3_plus_rational_cusp':144*16,
 'F7_single4':396,'F7_quadratic_pairs_including_repetition':300,'F7_total':6264,
 'Newton_closed_F4':[16,8,0,102],'Newton_closed_F7':[16,24,144,400],
 'gcd_odd_J_orders':3*h,'global_cusp_group_order':h,'quotient_order_divides':3,
 'all_60_matches_recomputed_from_actual_relations':True,'unmultiplied_global_classes':12,
 'all12_one_full_Aut48_orbit':True,
 'coarse_input_sha256':hashlib.sha256((P/'diagonal_3_15_D4_arithmetic.txt').read_bytes()).hexdigest(),
 'necessity_of_p13_p19_point_refinements':False,
 'boundary':'Native exact point-to-class maps are those computed in the bound Magma source. This checker independently verifies every support pattern/index, Newton count and integer-class match.'}
(P/'replay_M3N15_support_coverage.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2));print('M3N15_FULL_SUPPORT_AND_60_PAIR_COVERAGE_REPLAY_DONE')
