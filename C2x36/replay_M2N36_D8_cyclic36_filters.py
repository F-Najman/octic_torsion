#!/usr/bin/env python3
"""Independent actual-lattice replay of the complete461-state quotient filter."""
import ast,hashlib,json,sys,time
from collections import Counter
from pathlib import Path
from sympy import Matrix,ilcm
from sympy.matrices.normalforms import hermite_normal_form
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'common'))
from finite_field5_checked import FF,EC,pack
P=Path(__file__).resolve().parent/'data';t0=time.monotonic();bound={}
# Cases 1..180 come from the checkpoint of the interrupted run, kept in this repository as
# filter_2_36_D8_cyclic36_fast_interrupted_arithmetic.txt (its name in the hashes below is the
# original one).  With the argument 'complete', the complete rerun of all 461 cases
# (filter_2_36_D8_cyclic36_fast_arithmetic.txt) is replayed instead of the prefix and tail.
COMPLETE=sys.argv[1:]==['complete']
RENAMED={} if COMPLETE else {'filter_2_36_D8_cyclic36_fast_arithmetic.txt':'filter_2_36_D8_cyclic36_fast_interrupted_arithmetic.txt'}
def read(n):
 raw=(P/RENAMED.get(n,n)).read_bytes();bound[n]=hashlib.sha256(raw).hexdigest()
 return ast.literal_eval(raw.decode().replace('[*','[').replace('*]',']'))
def tatexy(F,c,d,Q):
 x,y=Q;aa3=F.sc(2,y);aa2=F.add(F.sc(3,x),c)
 aa4=F.add(F.add(F.sc(3,F.mul(x,x)),F.sc(2,F.mul(c,x))),d)
 slope=F.div(aa4,aa3);aaa1=F.sc(2,slope);aaa2=F.sub(aa2,F.mul(slope,slope))
 b=F.neg(F.div(F.pow(aaa2,3),F.mul(aa3,aa3)))
 cc=F.neg(F.div(F.sub(F.mul(aaa1,aaa2),aa3),aa3))
 r=F.div(b,cc);s=F.div(F.mul(cc,cc),F.sub(b,cc))
 tt=F.add(F.sub(F.mul(r,s),F.sc(2,r)),1)
 return F.div(F.sub(s,r),tt),F.div(tt,F.add(F.sub(F.sub(F.mul(s,s),s),r),1))
def isog2(F,E,Q,r):
 A=F.add(E.c,F.sc(3,r));B=F.add(F.add(F.sc(3,F.mul(r,r)),F.sc(2,F.mul(E.c,r))),E.d)
 x=F.sub(Q[0],r);xx=F.mul(x,x)
 return F.neg(F.sc(2,A)),F.sub(F.mul(A,A),F.sc(4,B)),(
  F.div(F.mul(Q[1],Q[1]),xx),F.div(F.mul(Q[1],F.sub(B,xx)),xx))
first=read('filter_2_36_D8_quotient_arithmetic.txt')
gd0,fields,singles,small,pairs,inputs=first;assert len(inputs)==461
own={D:[r for r in singles if r[0]==D] for D in range(3,9)}
sm={D:[r for r in small if r[0]==D] for D in [3,4,5]}
def supports(desc):
 if desc[0]==1:
  rows=[r for r in own[desc[1]] if r[1]==desc[2]];assert len(rows)==1;return rows
 return [own[desc[1]][desc[2]-1],sm[desc[3]][desc[4]-1]]
full=read('M2N36_cyclic36_full_group_complements.txt');gd,comps,qinv,qrel,qimages=full
inv,rel,osdegrees,hgens,hinv,native_comps=gd
assert inv==[4,4,12,12,36,7812,7812,0] and hinv==[12,252,7812]
assert rel==[[inv[i] if i==j else 0 for j in range(8)] for i in range(7)]
assert all(v[-1]==0 for v in rel+hgens)
A=Matrix([r[:-1] for r in rel+hgens]);H=hermite_normal_form(A.T)
assert H.shape==(7,7) and abs(int(H.det()))==214272
Hi=H.inv();den=ilcm(*[a.q for a in Hi]);MM=(den*Hi).applyfunc(int)
def key(v):
 assert len(v)==8 and v[-1]==0
 return tuple(int(a)%den for a in MM*Matrix(v[:-1]))
assert all(key(r)==(0,)*7 for r in rel+hgens)
compkeys=[{key(a) for a in rows} for rows in comps]
assert [len(a) for a in compkeys]==[1,7,27]
assert [len(a) for a in native_comps]==[1,7,27]
Fdict={D:FF(f) for D,f in fields}
def check_weight(mapno,D,sourcekey,e,w):
 F=Fdict[D];b,c=pack(sourcekey[:D]),pack(sourcekey[D:]);E=EC(F,b,c);Q=(b,b)
 if mapno==1:xy=tatexy(F,c,E.d,Q)
 else:
  if mapno==2:r=0
  else:
   Q2=None
   for _ in range(18):Q2=E.add(Q2,Q)
   r=F.sub(F.neg(c),Q2[0])
  cp,dp,qp=isog2(F,E,Q,r);xy=tatexy(F,cp,dp,qp)
 assert F.degree(xy)==e and D==e*w
legacy=read('filter_2_36_D8_cyclic36_arithmetic.txt')
prefix=read('filter_2_36_D8_cyclic36_fast_arithmetic.txt')
tail=read('filter_2_36_D8_cyclic36_fast_tail181_arithmetic.txt') if not COMPLETE else [prefix[0],[],[],prefix[3]]
assert legacy[0]==prefix[0]==tail[0]==gd
assert legacy[2]==prefix[2][:80] and len(legacy[2])==80 and len(prefix[2])==(461 if COMPLETE else 180) and len(tail[2])==(0 if COMPLETE else 281)
images={};checks=prefix[2]+tail[2]
for block in [prefix,tail]:
 for ky,e,w,cl in block[1]:
  m,D=ky[:2];sourcekey=ky[2:];check_weight(m,D,sourcekey,e,w)
  kk=tuple(ky)
  if kk in images:assert images[kk]==cl
  else:images[kk]=cl
assert [r[0] for r in checks]==list(range(1,462))
remaining=[];patterns=Counter()
for k,desc,cd,tests,flag in checks:
 assert desc==inputs[k-1];sp=supports(desc);assert cd==8-sum(r[0] for r in sp)
 assert cd in [0,1,2];results=[]
 for mapno,cl,hits in tests:
  assert mapno==len(results)+1
  sumcl=[sum(images[tuple([mapno,r[0]]+r[2])][j] for r in sp) for j in range(8)]
  assert all((sumcl[j]-cl[j])%inv[j]==0 for j in range(7)) and sumcl[-1]==cl[-1]==0
  wanted=tuple(-a%den for a in key(cl)) in compkeys[cd]
  assert bool(hits)==wanted;results.append(wanted)
  if not wanted:assert len(results)==len(tests)
 assert (len(results)==3 and all(results))==bool(flag)
 if flag:remaining.append(desc)
 patterns[tuple(r[0] for r in sp)]+=1
assert remaining==[[1,8,1606]] and prefix[3]==(remaining if COMPLETE else []) and tail[3]==remaining
for name,marker in [(('filter_2_36_D8_cyclic36_fast.log','M2N36_CYCLIC36_FILTER_DONE 461 1 complete=true') if COMPLETE else ('filter_2_36_D8_cyclic36_fast_tail181.log','M2N36_CYCLIC36_FILTER_DONE 281 1 complete=true')),
 ('reference_2_36_D8_cyclic36_maps.log','M2N36_CYCLIC36_NATIVE_REFERENCE_DONE')]:
 s=(P.parent/'logs'/name).read_text();assert marker in s and 'Runtime error' not in s and 'Assertion failed' not in s
out={'status':'complete necessary-condition filtering; one survivor remains',
 'source_group':[2,36],'source_genus':41,'target_genus':17,'tested':461,
 'all_source_cases_linked_to_first_stage':True,'support_patterns':{str(k):v for k,v in patterns.items()},
 'lattice_index':214272,'rational_Jacobian_order':12*252*7812,
 'target_cusp_complement_classes':[1,7,27],
 'all_weighted_moduli_degrees_independently_recomputed':True,
 'all_actual_group_memberships_replayed':True,'remaining':remaining,
 'verification_boundary':'Magma supplies the proved finite class presentation, exact divisor-class images and RR ideal identities. Python independently verifies source coverage assembly, finite-field isogeny moduli degrees, all class sums, actual integer-lattice membership, every allowed cusp complement and every decision.',
 'seconds':round(time.monotonic()-t0,3),'sha256':bound}
(P/('replay_M2N36_D8_cyclic36_filters_complete.json' if COMPLETE else 'replay_M2N36_D8_cyclic36_filters.json')).write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2));print('M2N36_D8_CYCLIC36_FILTER_REPLAY_DONE')
