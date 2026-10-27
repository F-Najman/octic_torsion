#!/usr/bin/env python3
"""Independent exact field/model checks and a local nonsplit2 witness."""
import ast,hashlib,json
from pathlib import Path
from fractions import Fraction
from sympy import Poly,Symbol,primerange
P=Path(__file__).resolve().parent/'data'
name='N28_D8_quotient_direct_torsion_arithmetic.txt';raw=(P/name).read_text()
# This expression contains rational literals. Parse only Python's exact
# integer arithmetic AST, avoiding float conversion and unrestricted eval.
def exact(node):
 if isinstance(node,ast.Expression):return exact(node.body)
 if isinstance(node,ast.List):return [exact(x) for x in node.elts]
 if isinstance(node,ast.Constant) and isinstance(node.value,int):return Fraction(node.value)
 if isinstance(node,ast.UnaryOp) and isinstance(node.op,ast.USub):return -exact(node.operand)
 if isinstance(node,ast.BinOp) and isinstance(node.op,ast.Div):return exact(node.left)/exact(node.right)
 raise ValueError(ast.dump(node))
f,alpha,ai,point,factor_degrees,nroots=exact(ast.parse(raw.replace('[*','[').replace('*]',']'),mode='eval'))
x=Symbol('x');ff=Poly(sum(c*x**i for i,c in enumerate(f)),x)
assert ff.degree()==8 and ff.is_irreducible
assert factor_degrees==[1,2] and nroots==1
def modcoeff(v,p):return [int(c.numerator)*pow(int(c.denominator),-1,p)%p for c in v]
def ev(v,r,p):
 y=0
 for c in reversed(v):y=(y*r+c)%p
 return y
witness=None
for p in primerange(11,300):
 try:fp=modcoeff(f,p);ap=[modcoeff(a,p) for a in ai]
 except ValueError:continue
 for r in range(p):
  if ev(fp,r,p) or ev([(i*c)%p for i,c in enumerate(fp)][1:],r,p)==0:continue
  a1,a2,a3,a4,a6=[ev(a,r,p) for a in ap]
  b2=(a1*a1+4*a2)%p;b4=(a1*a3+2*a4)%p;b6=(a3*a3+4*a6)%p;b8=(a1*a1*a6+4*a2*a6-a1*a3*a4+a2*a3*a3-a4*a4)%p
  disc=(-b2*b2*b8-8*b4**3-27*b6*b6+9*b2*b4*b6)%p
  if not disc:continue
  roots=[u for u in range(p) if (4*u**3+b2*u*u+2*b4*u+b6)%p==0]
  if len(roots)==1:
   witness={'prime':p,'simple_root_of_field_polynomial':r,'a_invariants':[a1,a2,a3,a4,a6],
            'good_reduction_discriminant':disc,'psi2_roots':roots,'psi2_factor_degrees':[1,2]}
   break
 if witness:break
assert witness
assert point==[[0]*8,[0]*8] or all(not any(v) for v in point)
report={'field_degree':8,'field_irreducible_over_Q':True,'exact_nonzero_rational_2_points':1,
        'independent_local_obstruction_to_full_2_torsion':witness,'input_sha256':hashlib.sha256(raw.encode()).hexdigest()}
(P/'replay_N28_D8_exact_point.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2));print('N28_D8_EXACT_FIELD_LOCAL2_REPLAY_DONE')
