#!/usr/bin/env python3
from pathlib import Path
from fractions import Fraction
import ast,json,hashlib
from sympy import Symbol,Poly,sqrt
P=Path(__file__).resolve().parent/'data'
raw=(P/'realization_3_15_direct_v4_root_arithmetic.txt').read_text()
def exact(n):
 if isinstance(n,ast.Expression):return exact(n.body)
 if isinstance(n,ast.List):return [exact(x) for x in n.elts]
 if isinstance(n,ast.Constant) and type(n.value)==int:return Fraction(n.value)
 if isinstance(n,ast.UnaryOp) and isinstance(n.op,ast.USub):return -exact(n.operand)
 if isinstance(n,ast.BinOp) and isinstance(n.op,ast.Div):return exact(n.left)/exact(n.right)
 raise ValueError(ast.dump(n))
f,ai,pts,reduced,disc,order=exact(ast.parse(raw.replace('[*','[').replace('*]',']'),mode='eval'))
x=Symbol('x');ff=Poly(x**4-5*x**3+8*x*x-6*x+3,x,extension=sqrt(-3));assert ff.is_irreducible and ff.discriminant()==-507
assert f==[[3,0],[-6,0],[8,0],[-5,0],[1,0]]
p=43
def rat(c):return int(c.numerator)*pow(int(c.denominator),-1,p)%p
def ev(v,z,w):return sum(sum(rat(c)*pow(z,j,p) for j,c in enumerate(a))*pow(w,i,p) for i,a in enumerate(v))%p
a1,a2,a3,a4,a6=map(int,reduced);assert [ev(a,6,4) for a in ai]==[a1,a2,a3,a4,a6]
assert (6*6+6+1)%p==0 and (4**4-5*4**3+8*4**2-6*4+3)%p==0 and (-507)%p and 3%p
points={(xx,yy) for xx in range(p) for yy in range(p) if (yy*yy+a1*xx*yy+a3*yy-xx**3-a2*xx*xx-a4*xx-a6)%p==0};assert len(points)+1==int(order)==45
b2=(a1*a1+4*a2)%p;b4=(a1*a3+2*a4)%p;b6=(a3*a3+4*a6)%p;b8=(a1*a1*a6+4*a2*a6-a1*a3*a4+a2*a3*a3-a4*a4)%p
assert (-b2*b2*b8-8*b4**3-27*b6*b6+9*b2*b4*b6)%p==int(disc)!=0
def add(A,B):
 if A is None:return B
 if B is None:return A
 x,y=A;u,v=B
 if x==u and (y+v+a1*x+a3)%p==0:return None
 if A==B:
  d=(2*y+a1*x+a3)%p;la=(3*x*x+2*a2*x+a4-a1*y)*pow(d,-1,p)%p;nu=(-x**3+a4*x+2*a6-a3*y)*pow(d,-1,p)%p
 else:
  d=(u-x)%p;la=(v-y)*pow(d,-1,p)%p;nu=(y*u-v*x)*pow(d,-1,p)%p
 xx=(la*la+a1*la-a2-x-u)%p;yy=(-(la+a1)*xx-nu-a3)%p;assert (xx,yy) in points;return xx,yy
def mult(n,A):
 R=None
 for _ in range(n):R=add(R,A)
 return R
pp=(ev(pts[0],6,4),ev(pts[1],6,4));qq=(ev(pts[2],6,4),ev(pts[3],6,4));assert pp in points and qq in points
assert mult(3,pp) is None and mult(15,qq) is None and mult(5,qq) is not None and mult(3,qq) is not None
assert len({add(mult(i,pp),mult(j,qq)) for i in range(3) for j in range(15)})==45
report={'degree':8,'group':[3,15],'relative_quartic_irreducible_over_Q_zeta3':True,'relative_polynomial':[3,-6,8,-5,1],'relative_polynomial_discriminant':-507,'good_unramified_split_prime':43,'zeta3_maps_to':6,'quartic_generator_maps_to':4,'reduced_a_invariants':[a1,a2,a3,a4,a6],'reduced_discriminant':int(disc),'independently_counted_reduction_order':45,'independent_reduction_subgroup_3_15_order':45,'full_global_torsion_exact':[3,15],'native_exact_level_check_required':True,'input_sha256':hashlib.sha256(raw.encode()).hexdigest()}
(P/'replay_3_15_direct_root.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2));print('M3N15_DIRECT_OCCURRENCE_INDEPENDENT_REPLAY_DONE')
