"""Independent Python finite fields and generalized elliptic group law."""
class FF:
 def __init__(self,p,f):self.p=p;self.f=f;self.d=len(f)-1;self.z=(0,)*self.d;self.o=self.c(1)
 def c(self,a):return (a%self.p,)+(0,)*(self.d-1)
 def add(self,a,b):return tuple((x+y)%self.p for x,y in zip(a,b))
 def neg(self,a):return tuple(-x%self.p for x in a)
 def sub(self,a,b):return self.add(a,self.neg(b))
 def mul(self,a,b):
  v=[0]*(2*self.d-1)
  for i,x in enumerate(a):
   for j,y in enumerate(b):v[i+j]=(v[i+j]+x*y)%self.p
  for i in range(len(v)-1,self.d-1,-1):
   for j in range(self.d):v[i-self.d+j]=(v[i-self.d+j]-v[i]*self.f[j])%self.p
  return tuple(v[:self.d])
 def pow(self,a,n):
  v=self.o
  while n:
   if n&1:v=self.mul(v,a)
   a=self.mul(a,a);n//=2
  return v
 def div(self,a,b):assert b!=self.z;return self.mul(a,self.pow(b,self.p**self.d-2))
 def sc(self,n,a):return self.mul(self.c(n),a)
 def key(self,qt):
  keys=[]
  for _ in range(self.d):keys.append(qt[0]+qt[1]);qt=tuple(self.pow(x,self.p) for x in qt)
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
