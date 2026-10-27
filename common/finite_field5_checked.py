"""Independent F5 extension-field and elliptic group-law arithmetic.
Copied locally from the already validated positive source checker; no tools,
no native Magma arithmetic, and no write side effects.
"""
import ast,hashlib,json,time
from array import array
from collections import Counter
from pathlib import Path
from sympy import Poly,Symbol,divisors
P=Path(__file__).resolve().parent;t0=time.monotonic()
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
