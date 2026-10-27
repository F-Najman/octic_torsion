#!/usr/bin/env python3
"""Independent field/elliptic-group/Aut48/count/coverage and Hecke assembly replay."""
import ast,hashlib,importlib.util,itertools,json,re
from collections import Counter
from fractions import Fraction
from pathlib import Path
from sympy import Poly,Symbol
P=Path(__file__).resolve().parent/'data'
shas={}
def read(n):
 raw=(P/n).read_bytes();shas[n]=hashlib.sha256(raw).hexdigest()
 return ast.literal_eval(raw.decode().replace('[*','[').replace('*]',']'))
records,coverage=read('preflight_3_24_D4_Aut48_input.txt')
assert len(records)==14 and Counter(r[0] for r in records)=={3:1,4:13}
g,fps=read('preflight_3_24_D4_Frob_arithmetic.txt');assert g==25
for p,fp,cl,order in fps:
 assert len(fp)==51 and fp[-1]==1 and sum(fp)==order
 tr=[];closed=[]
 for n in range(1,5):
  tr.append(-n*fp[50-n]-sum(fp[50-j]*tr[n-j-1] for j in range(1,n)))
  v=p**n+1-tr[-1]-sum(d*closed[d-1] for d in range(1,n) if n%d==0)
  assert v%n==0;closed.append(v//n)
 assert closed==cl
assert fps[0][0]==7 and fps[0][2]==[24,12,48,456]
spec=importlib.util.spec_from_file_location('tate',P.parents[1]/'common'/'cusp_conditions_relative.py');tate=importlib.util.module_from_spec(spec);spec.loader.exec_module(tate)
cusps=tate.compute(3,24,7,4)
assert Counter(cusps['closed_cusp_degrees'])=={1:24,2:12}
assert cusps['allowed_residues_by_complement_degree'][0]==[1,7,13,19]
assert cusps['rational_residue_cusps_not_F_rational']==0
assert Fraction(325*576,32768)>5
mats=[t for t in itertools.product(range(3),repeat=4) if (t[0]*t[3]-t[1]*t[2])%3==1]
assert len(mats)==24
pars={(a,b,c,d,u) for a,b,c,d in mats for u in [1,3]}
for a,b,c,d,u in pars:
 for A,B,C,D,U in pars:
  out=((A*a+B*c)%3,(A*b+B*d)%3,(C*a+D*c)%3,(C*b+D*d)%3,U*u%8)
  if out[4] not in [1,3]:out=tuple((-x)%3 for x in out[:4])+((-out[4])%8,)
  assert out in pars
class FF:
 def __init__(self,f):self.f=f;self.d=len(f)-1;self.z=(0,)*self.d;self.o=self.c(1)
 def c(self,a):return (a%7,)+(0,)*(self.d-1)
 def add(self,a,b):return tuple((x+y)%7 for x,y in zip(a,b))
 def neg(self,a):return tuple(-x%7 for x in a)
 def sub(self,a,b):return self.add(a,self.neg(b))
 def mul(self,a,b):
  v=[0]*(2*self.d-1)
  for i,x in enumerate(a):
   for j,y in enumerate(b):v[i+j]=(v[i+j]+x*y)%7
  for i in range(len(v)-1,self.d-1,-1):
   for j in range(self.d):v[i-self.d+j]=(v[i-self.d+j]-v[i]*self.f[j])%7
  return tuple(v[:self.d])
 def pow(self,a,n):
  v=self.o
  while n:
   if n&1:v=self.mul(v,a)
   a=self.mul(a,a);n//=2
  return v
 def div(self,a,b):assert b!=self.z;return self.mul(a,self.pow(b,7**self.d-2))
 def sc(self,n,a):return self.mul(self.c(n),a)
 def key(self,qt):
  keys=[]
  for _ in range(self.d):keys.append(qt[0]+qt[1]);qt=tuple(self.pow(x,7) for x in qt)
  return min(keys)
class EC:
 def __init__(self,F,q,t):
  self.F=F;self.a=F.add(F.sub(F.mul(q,t),q),F.add(t,F.c(2)))
  self.b=F.add(F.mul(q,F.mul(t,t)),F.mul(F.sub(F.o,q),t))
  assert self.b!=F.z and F.sub(F.pow(self.a,3),F.sc(27,self.b))!=F.z
 def add(self,P,Q):
  F=self.F
  if P is None:return Q
  if Q is None:return P
  x,y=P;X,Y=Q
  if x==X:
   if F.add(F.add(y,Y),F.add(F.mul(self.a,x),self.b))==F.z:return None
   slope=F.div(F.sub(F.sc(3,F.mul(x,x)),F.mul(self.a,y)),F.add(F.sc(2,y),F.add(F.mul(self.a,x),self.b)))
  else:slope=F.div(F.sub(Y,y),F.sub(X,x))
  nu=F.sub(y,F.mul(slope,x));xx=F.sub(F.add(F.mul(slope,slope),F.mul(self.a,slope)),F.add(x,X))
  yy=F.sub(F.neg(F.mul(F.add(slope,self.a),xx)),F.add(nu,self.b));return xx,yy
 def mul(self,n,P):
  Q=None
  while n:
   if n&1:Q=self.add(Q,P)
   P=self.add(P,P);n//=2
  return Q
 def normal(self,P,Q):
  F=self.F;r,tt=P;xx,yy=Q
  sl=F.div(F.sub(F.sc(3,F.mul(r,r)),F.mul(self.a,tt)),F.add(F.sc(2,tt),F.add(F.mul(self.a,r),self.b)))
  aa=F.add(self.a,F.sc(2,sl));bb=F.add(self.b,F.add(F.mul(self.a,r),F.sc(2,tt)))
  x=F.sub(xx,r);y=F.sub(F.sub(yy,tt),F.mul(sl,x));scale=F.div(y,F.mul(x,x))
  t=F.neg(F.mul(F.mul(scale,scale),x));q=F.div(F.sub(F.sub(F.mul(scale,aa),t),F.c(2)),F.sub(t,F.o))
  assert F.mul(F.pow(scale,3),bb)==F.add(F.mul(q,F.mul(t,t)),F.mul(F.sub(F.o,q),t))
  return q,t
cm={}
for d,keys,labels,orbs,idx in coverage:
 assert d in [3,4] and len(keys)=={3:48,4:456}[d] and len(idx)=={3:1,4:13}[d]
 assert len(set(map(tuple,keys)))==len(keys) and len(labels)==len(keys)
 assert all(len(k)==2*d and set(k)<={0,1,2,3,4,5,6} for k in keys)
 seen={}
 for j,orb in enumerate(orbs,1):
  assert len(orb)==48
  for k in set(map(tuple,orb)):assert k not in seen;seen[k]=j
 assert set(seen)==set(map(tuple,keys)) and [seen[tuple(k)] for k in keys]==labels
 assert all(labels[i-1]==j for j,i in enumerate(idx,1))
 cm[d]=(keys,labels,orbs,idx)
paired_roots=[]
for d,i,f,qt in records:
 assert len(f)==d+1 and f[-1]==1 and Poly.from_list(f[::-1],Symbol('x'),modulus=7).is_irreducible
 F=FF(f);qt=tuple(map(tuple,qt));q,t=qt;E=EC(F,q,t);O=(F.z,F.z);Q=(F.neg(t),F.mul(t,t))
 assert E.mul(3,O) is None and E.mul(24,Q) is None and E.mul(12,Q) is not None and E.mul(8,Q) is not None
 R3=E.mul(8,Q);R8=E.mul(3,Q);assert O not in [R3,E.mul(2,R3)]
 x3,y3=R3;sl=F.div(F.sub(F.sc(3,F.mul(x3,x3)),F.mul(E.a,y3)),F.add(F.sc(2,y3),F.add(F.mul(E.a,x3),E.b)))
 wp=F.div(F.neg(y3),F.sub(F.mul(sl,x3),y3));assert wp in [F.c(2),F.c(4)];paired_roots.append(wp[0])
 qq=qt;conj=[]
 for _ in range(d):conj.append(qq);qq=tuple(F.pow(v,7) for v in qq)
 assert len(set(conj))==d and qq==qt
 assert E.add(E.mul(2,R3),E.mul(3,R8))==Q
 keys,labels,orbs,idx=cm[d];assert F.key(qt)==tuple(keys[idx[i-1]-1])
 images=[]
 for a,b,c,dd in mats:
  PP=E.add(E.mul(a,O),E.mul(b,R3));RR=E.add(E.mul(c,O),E.mul(dd,R3))
  for u in [1,3]:images.append(F.key(E.normal(PP,E.add(E.mul(2,RR),E.mul(3*u,R8)))))
 assert set(images)==set(map(tuple,orbs[i-1]))
assert len(set(paired_roots))==1
tests,remaining=read('sieve_3_24_D4_Aut48_arithmetic.txt');assert remaining==[]
groups={}
for n,d,i,q,flag in tests:
 assert records[n-1][:2]==[d,i] and q in [13,19] and flag in [0,1]
 groups.setdefault(n,[]).append((q,flag))
assert set(groups)==set(range(1,15))
assert all(v in [[(13,0)],[(13,1),(19,0)]] for v in groups.values())
log=(P.parent/'logs'/'sieve_3_24_D4_Aut48.log').read_text()
assert 'complete=true' in log and 'Runtime error' not in log and 'Assertion failed' not in log
logged=[tuple(map(int,[n,d,i,q]))+(flag=='true',) for n,d,i,q,flag in re.findall(r'M3N24_D4_FULLAUT_HECKE_TEST (\d+) (\d+) (\d+) (\d+) (true|false)',log)]
assert logged==[tuple(r[:4])+(bool(r[4]),) for r in tests]
ref=(P.parent/'logs'/'reference_3_24_D4_T13.log').read_text();assert 'M3N24_D4_T13_GENERIC_REFERENCE_EQUAL 42' in ref and 'complete=true' in ref
report={'degree':8,'group':[3,24],'relative_degree':4,'genus':25,'NC_counts':[0,0,48,456],
 'Aut_order':48,'support_representatives':[1,13],'actual_keys':504,
 'all48_elliptic_images_independently_recomputed':True,'rational_residue_cusps_F_rational':24,
 'Weil_pairing_P_8Q_in_F7':paired_roots[0],
 'all_records_excluded':14,'remaining':0,'generic_T13_reference_degree':42,
 'verification_boundary':'Python independently replays exact finite fields, elliptic additions, point orders/independence, all48 canonical automorphism images for every representative, full actual-key partitions, Frobenius counts, Tate global cusp action and every logged Hecke decision. Native exact divisor arithmetic/IsPrincipal are not reimplemented in Python.',
 'input_sha256':shas}
(P/'replay_M3N24_D4.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2));print('M3N24_D4_COMPLETE_EXCLUSION_REPLAY_DONE')
