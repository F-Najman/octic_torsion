#!/usr/bin/env python3
"""Positive finite-point certificate; no dependence on the producer's j search.

Each elliptic level structure is checked from scratch, all24 images are
expanded and counted, and their disjoint union has the independent
Frobenius cardinality in every support degree3..8.
"""
import ast,hashlib,json,time
from array import array
from collections import Counter
from pathlib import Path
from sympy import Poly,Symbol,divisors
P=Path(__file__).resolve().parent/'data';t0=time.monotonic()
def read(n):return ast.literal_eval((P/n).read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
def pack(v):return sum(int(a)*5**i for i,a in enumerate(v))
def digs(a,d):return [(a//5**i)%5 for i in range(d)]
# Two four-digit chunks give constant-time addition in all our fields.
AD=array('H',[0])*(625*625);NG=array('H',[0])*625
for a in range(625):
 av=digs(a,4);NG[a]=pack([-x%5 for x in av])
 for b in range(625):AD[625*a+b]=pack([(x+(b//5**i)%5)%5 for i,x in enumerate(av)])
class FF:
 def __init__(self,f):
  self.f=f;self.d=len(f)-1;self.q=5**self.d;self.n=self.q-1
  assert f[-1]==1 and Poly.from_list(f[::-1],Symbol('x'),modulus=5).is_irreducible
  self.ex=array('I',[0])*self.n;self.lg=array('i',[-1])*self.q;v=1;topbase=5**(self.d-1)
  for i in range(self.n):
   assert self.lg[v]==-1;self.ex[i]=v;self.lg[v]=i;top=v//topbase;v=(v%topbase)*5
   if top:v=pack([(x-top*f[j])%5 for j,x in enumerate(digs(v,self.d))])
  assert v==1 and min(self.lg[1:])>=0
  self.rv=array('I',(pack(digs(a,self.d)[::-1]) for a in range(self.q)))
 def add(self,a,b):return AD[(a%625)*625+b%625]+625*AD[(a//625)*625+b//625]
 def neg(self,a):return NG[a%625]+625*NG[a//625]
 def sub(self,a,b):return self.add(a,self.neg(b))
 def mul(self,a,b):return self.ex[(self.lg[a]+self.lg[b])%self.n] if a and b else 0
 def div(self,a,b):assert b;return self.ex[(self.lg[a]-self.lg[b])%self.n] if a else 0
 def pow(self,a,n):return self.ex[(self.lg[a]*n)%self.n] if a else int(n==0)
 def sc(self,n,a):return self.mul(n%5,a)
 def score(self,bc):return self.rv[bc[0]],self.rv[bc[1]]
 def key(self,bc):
  best=bc
  for _ in range(self.d-1):
   bc=(self.pow(bc[0],5),self.pow(bc[1],5))
   if self.score(bc)<self.score(best):best=bc
  return best
 def degree(self,bc):return next(d for d in divisors(self.d) if all(self.pow(a,5**d)==a for a in bc))
class EC:
 def __init__(self,F,b,c):
  self.F=F;self.c=c;self.d=F.mul(b,F.sub(F.sub(1,b),c))
  assert self.d and F.sub(F.mul(c,c),F.sc(4,self.d))
 def add(self,P,Q):
  if P is None:return Q
  if Q is None:return P
  F=self.F;x,y=P;X,Y=Q
  if x==X:
   if F.add(y,Y)==0:return None
   slope=F.div(F.add(F.add(F.sc(3,F.mul(x,x)),F.sc(2,F.mul(self.c,x))),self.d),F.sc(2,y))
  else:slope=F.div(F.sub(Y,y),F.sub(X,x))
  z=F.sub(F.sub(F.mul(slope,slope),self.c),F.add(x,X))
  return z,F.sub(F.mul(slope,F.sub(x,z)),y)
 def bc(self,P,Q):
  F=self.F;assert P[1]==0;x=F.sub(Q[0],P[0]);yy=F.mul(Q[1],Q[1]);assert x and yy
  xx=F.mul(x,x)
  return F.div(F.mul(xx,x),yy),F.div(F.mul(F.add(self.c,F.sc(3,P[0])),xx),yy)
src=read('preflight_2_36_D8_Frob_arithmetic.txt')[0];assert src[0]==[72,37] and src[1]==41
fp,closed=src[2:4];tr=[];calc=[]
for n in range(1,9):
 tr.append(-n*fp[82-n]-sum(fp[82-j]*tr[n-j-1] for j in range(1,n)))
 z=5**n+1-tr[-1]-sum(d*calc[d-1] for d in range(1,n) if n%d==0);assert z%n==0;calc.append(z//n)
assert calc==closed==[24,8,48,120,720,2908,10416,50334]
expected={1:0,2:0,3:48,4:120,5:720,6:2904,7:10416,8:50334}
report=[];small=[];shas={}
for D in range(3,9):
 name=f'M2N36_D8_direct_full{D}_arithmetic.txt';raw=(P/name).read_bytes();shas[name]=hashlib.sha256(raw).hexdigest();a=read(name)
 if D==8:f,complete,jcount,curves,g4,g9,counts,records=a
 else:dd,f,complete,jcount,curves,g4,g9,counts,records=a;assert dd==D
 assert complete and len(f)==D+1
 F=FF(f);allpoints={};orbitkeys=set();ownrecords=0
 for d,ky in records:
  assert len(ky)==2*D and set(ky)<={0,1,2,3,4}
  b,c=pack(ky[:D]),pack(ky[D:]);assert F.degree((b,c))==d
  E=EC(F,b,c);P2=(0,0);Q=(b,b);qs=[None]
  for _ in range(36):qs.append(E.add(qs[-1],Q))
  assert qs[36] is None and qs[18] is not None and qs[12] is not None and qs[18]!=P2
  ps=[P2,E.add(P2,qs[18])];assert all(p[1]==0 for p in ps)
  ims={F.key(E.bc(p,E.add(P2 if h else None,qs[u]))) for p in ps for h in [0,1] for u in [1,5,7,11,13,17]}
  assert min(ims,key=F.score)==(b,c) and (b,c) not in orbitkeys;orbitkeys.add((b,c))
  for im in ims:assert im not in allpoints;allpoints[im]=d
  if d==D:ownrecords+=1
 got=[sum(d==e for d in allpoints.values()) for e in divisors(D)]
 assert got==counts==[expected[d] for d in divisors(D)]
 assert all(len(ky)==2*D for d,ky in records)
 if D<=5:
  small.append([D,f,[[digs(b,D),digs(c,D)] for b,c in sorted(allpoints,key=F.score)]])
 report.append({'field_degree':D,'all_records':len(records),'own_degree_representatives':ownrecords,'all_closed_points':got,'checked_all24_images':True})
 print('M2N36_DIRECT_FIELD_REPLAY',D,len(records),ownrecords,got,round(time.monotonic()-t0,3))
def magma(x):
 if isinstance(x,list):return '[* '+', '.join(magma(y) for y in x)+' *]'
 return str(x)
(P/'M2N36_D8_direct_small_all_points.txt').write_text(magma(small)+'\n')
out={'degree':8,'group':[2,36],'finite_field_prime':5,'source_genus':41,'NC_counts':[expected[d] for d in range(1,9)],'fields':report,'total_support_orbits':sum(r['own_degree_representatives'] for r in report),'small_supports_expanded':[48,120,720],'status':'complete positive source support certificate; no global exclusion claimed','verification_boundary':'Every curve is non-singular, Q has exactorder36, P2 is independent, everyAut24/Frobenius image is recomputed with polynomial-defined finitefield logtables and no Magma arithmetic. Disjoint unions match independent modular-symbol Frobenius counts. The j/twist producer search is not trusted for completeness.','input_sha256':shas}
(P/'replay_M2N36_D8_direct.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2));print('M2N36_D8_DIRECT_POSITIVE_REPLAY_DONE')
