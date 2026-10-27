class BinaryField:
 def __init__(self,poly):
  self.d=len(poly)-1;self.q=1<<self.d;self.pol=sum(a<<i for i,a in enumerate(poly));self.square=[self.mul(a,a) for a in range(self.q)]
  class Inv(dict):
   def __missing__(s,a):
    assert a;z=self.pow(a,self.q-2);s[a]=z;return z
  self.inv=Inv()
 def mul(self,a,b):
  out=0
  while b:
   if b&1:out^=a
   b>>=1;a<<=1
   if a&self.q:a^=self.pol
  return out
 def div(self,a,b):return self.mul(a,self.inv[b])
 def pow(self,a,n):
  out=1
  while n:
   if n&1:out=self.mul(out,a)
   a=self.square[a];n>>=1
  return out
 def trace(self,a):
  out=0
  for _ in range(self.d):out^=a;a=self.square[a]
  assert out in [0,1];return out

class Elliptic:
    def __init__(self,F,a):self.F=F;self.a1,self.a2,self.a3,self.a4,self.a6=a
    def on(self,P):
        if P is None:return True
        F=self.F;m=F.mul;x,y=P
        return F.square[y]^m(self.a1,m(x,y))^m(self.a3,y)==m(x,F.square[x])^m(self.a2,F.square[x])^m(self.a4,x)^self.a6
    def add(self,P,Q):
        if P is None:return Q
        if Q is None:return P
        F=self.F;m=F.mul;x,y=P;u,v=Q
        if x==u:
            if y^v^m(self.a1,x)^self.a3==0:return None
            assert P==Q
            lam=F.div(F.square[x]^self.a4^m(self.a1,y),m(self.a1,x)^self.a3)
        else:lam=F.div(y^v,x^u)
        nu=y^m(lam,x);xx=F.square[lam]^m(self.a1,lam)^self.a2^x^u
        yy=m(lam^self.a1,xx)^nu^self.a3
        out=(xx,yy);assert self.on(out);return out
    def multiple(self,n,P):
        out=None
        while n:
            if n&1:out=self.add(out,P)
            P=self.add(P,P);n>>=1
        return out
    def tate(self,P):
        F=self.F;m=F.mul;x,y=P
        aa3=self.a3^m(self.a1,x);aa2=x^self.a2
        aa4=F.square[x]^self.a4^m(self.a1,y)
        c0=F.div(aa4,aa3);a2=aa2^F.square[c0]^m(self.a1,c0)
        b=F.div(m(a2,F.square[a2]),F.square[aa3]);c=F.div(m(self.a1,a2)^aa3,aa3)
        r=F.div(b,c);s=F.div(F.square[c],b^c);t=m(r,s)^1
        xy=(F.div(s^r,t),F.div(t,F.square[s]^s^r^1))
        return (b,c),xy
