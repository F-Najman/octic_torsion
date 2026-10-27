"""Independent symbolic, finite-field and level-group completion replay.

No Magma, external curve data, previous-degree output or retained state.
The two exact rational rank bounds are separately native descent proofs.
"""
from math import isqrt
import sympy as S

s=S.symbols('s');I=S.I
def zero(a):assert S.cancel(a)==0
a=2*s/(1+s*s);b=(1-s*s)/(1+s*s);lam=a*a
P=(a,I*a*(a-1));Q=(1+b,b*(1+b));a2=-(1+lam)
def xadd(P,Q):
    x,y=P;u,v=Q
    return S.cancel(((v-y)/(u-x))**2-a2-x-u)
xplus=xadd(P,Q);xminus=xadd(P,(Q[0],-Q[1]))
w=S.cancel(4*(xminus-lam));expected=8*I*s*(s*s-1)/(s*s+1)**2
zero(w-expected);zero(-4*(xplus-lam)-w)
j4=lambda w:16*(w*w+16)**3/w**4
zero(256*(lam*lam-lam+1)**3/(lam*lam*(lam-1)**2)-j4(w))
assert max(S.degree(S.fraction(w)[0],s),S.degree(S.fraction(w)[1],s))==4
print('C420_CARTAN_SYMBOLIC_PASS degree4 both_generators jidentity')

bb=S.symbols('b');t5=bb-11-1/bb
c4=bb**4-12*bb**3+14*bb*bb+12*bb+1;disc=bb**5*(bb*bb-11*bb-1)
j5=lambda t:(t*t+10*t+5)**3/t
zero(c4**3/disc-j5(t5));zero(t5.subs(bb,-1/bb)-t5)
# Exact normalization at2P5=(b,b^2), using the general Tate formulas.
a1,a2,a3=1-bb,-bb,-bb;x,y=bb,bb*bb
aa3=2*y+a3+a1*x;aa2=3*x+a2;aa4=3*x*x+2*x*a2-a1*y
aaa1=2*aa4/aa3+a1;aaa2=aa2-(aa4/aa3)**2-a1*aa4/aa3
b2=S.cancel(-aaa2**3/aa3**2);c2=S.cancel(-(aaa1*aaa2-aa3)/aa3)
zero(b2+1/bb);zero(c2-b2)
print('C420_CYCLIC5_SYMBOLIC_PASS actual2P normalization and jidentity')

x,y=S.symbols('x y');ee=(x-3)*(x+1)*(x+2)
t=(x-3)**2*(x+2)/(x+1)**2;w=y*(x+1)**2/(4*(x+2))
num=S.fraction(S.cancel(16*(w*w+16)**3*t-w**4*(t*t+10*t+5)**3))[0]
assert S.rem(num,y*y-ee,y)==0
w2=ee*(x+1)**4/(16*(x+2)**2)
u=S.cancel(w2*(t*t+10*t+5)/(4*(w2+16)))
ss=S.cancel((u*u*(2*t+10)-16)/(4*(u*u+2*u-4)));xx=S.cancel(2*ss+1)
zero(xx-x);zero(4*w*(xx+2)/(xx+1)**2-y)
jnum,jden=S.fraction(S.cancel(j5(t)))
assert S.degree(jnum,x)-S.degree(jden,x)==5
assert S.factor(jden)==(x-3)**2*(x+1)**10*(x+2)
print('C420_BIRATIONAL_SYMBOLIC_PASS both_inverses cusp_poles=[4,20,2,10]')

def ga(a,b):return ((a[0]+b[0])%3,(a[1]+b[1])%3)
def gn(a):return ((-a[0])%3,(-a[1])%3)
def gm(a,b):return ((a[0]*b[0]-a[1]*b[1])%3,(a[0]*b[1]+a[1]*b[0])%3)
def gi(a):
    n=(a[0]*a[0]+a[1]*a[1])%3;assert n
    z=pow(n,-1,3);return ((a[0]*z)%3,(-a[1]*z)%3)
def gc(n):return (n%3,0)
def add(P,Q):
    if P is None:return Q
    if Q is None:return P
    x,y=P;u,v=Q
    if x==u:
        if ga(y,v)==gc(0):return None
        slope=gm(gc(-7),gi(gm(gc(2),y)))
    else:slope=gm(ga(v,gn(y)),gi(ga(u,gn(x))))
    xx=ga(ga(gm(slope,slope),gn(x)),gn(u));yy=ga(gn(y),gm(slope,ga(x,gn(xx))))
    return xx,yy
elts=[(a,b) for a in range(3) for b in range(3)]
points=[None]+[(x,y) for x in elts for y in elts if gm(y,y)==ga(gm(gm(x,x),x),gm(gc(-7),x))]
assert len(points)==16 and sum(add(P,P) is None for P in points)==4
assert all(add(add(P,P),add(P,P)) is None for P in points)
print('C420_FINITE_GROUP_INDEPENDENT_PASS E(F9)=[4,4]')

pts=[(-2,0),(-1,0),(3,0)]+[(-1+2*ep*I,et*(4-2*ep*I)) for ep in [-1,1] for et in [-1,1]]
assert len(set(pts))==7
for px,py in pts:
    zero(py*py-ee.subs(x,px))
    if py!=0:
        tx=S.simplify(t.subs(x,px));wx=S.simplify(w.subs({x:px,y:py}))
        assert tx in [-11-2*I,-11+2*I] and wx in [-2*I,2*I]
        zero(j5(tx)-1728);zero(j4(wx)-1728)
assert isqrt(5)**2!=5
print('C420_GLOBAL_POINT_INDEPENDENT_PASS four_cusps four_j1728 torsion_halving_criterion')

def norm(M,u):
    M=[a%4 for a in M];u%=5
    if u in [3,4]:M=[-a%4 for a in M];u=-u%5
    return tuple(M+[u])
def mul(p,q):
    a,b,c,d,u=p;A,B,C,D,U=q
    return norm([a*A+b*C,a*B+b*D,c*A+d*C,c*B+d*D],u*U)
H1={(1,0,0,1,1),(3,2,2,3,1),(0,3,1,0,2),(2,1,3,2,2)}
H2={(1,0,0,1,1),(1,2,0,1,1),(1,3,2,3,2),(1,1,2,3,2)}
T=(1,0,1,1,1);Ti=(1,0,3,1,1)
assert {mul(mul(T,h),Ti) for h in H1}==H2
assert {mul(a,b) for a in H1 for b in H1}==H1
print('C420_SUBGROUP_CONJUGACY_INDEPENDENT_PASS T=[[1,0],[1,1]]')
print('C420_COMPACT_REPLAY_DONE complete=true')
