#!/usr/bin/env python3
"""Independent binary and GF8 replay of the C57 U3 Fourier certificate."""
import ast
from pathlib import Path
import re
import sys
s=Path(sys.argv[1]).read_text()
def val(name):
 x=re.findall(rf'\b{name}:=(.*?);',s,re.S);assert len(x)==1,name
 return ast.literal_eval(x[0])
def matrix(pattern):
 x=re.findall(pattern+r':=Matrix\(GF\(2\),(\d+),(\d+),(\[.*?\])\);',s,re.S);assert len(x)==1,pattern
 nr,nc=map(int,x[0][:2]);a=ast.literal_eval(x[0][2]);assert len(a)==nr*nc and set(a)<={0,1}
 return [sum(a[i*nc+j]<<j for j in range(nc)) for i in range(nr)],nc
def apply(x,rows):
 out=0
 while x:
  b=x&-x;out^=rows[b.bit_length()-1];x^=b
 return out
def mul(a,b):return [apply(x,b) for x in a]
def pivots(a):
 p={}
 for x in a:
  while x:
   k=x.bit_length()-1
   if k in p:x^=p[k]
   else:p[k]=x;break
 return p
def contains(p,x):
 while x:
  k=x.bit_length()-1
  if k not in p:return False
  x^=p[k]
 return True
assert val('N')==57 and val('gx')==25
n=val('CheckedCoefficients');orig=val('OriginalCoefficients');assert n==20 and orig==3*n and 3*n>48
mod=val('F8ModulusCoefficients');assert len(mod)==4 and mod[0]==mod[3]==1 and sum(mod)%2==1
modbits=sum(x<<i for i,x in enumerate(mod))
def fm(a,b):
 out=0
 while b:
  if b&1:out^=a
  b>>=1;a<<=1
  if a&8:a^=modbits
 return out
def sqrt(a):return fm(fm(a,a),fm(a,a))
full,nfull=matrix(r'\bFullSourceQExp');B,nb=matrix(r'\bSourceQExp');T,nt=matrix(r'\bTargetQExp')
C,nc=matrix(r'\bCartierMatrix');H,nh=matrix(r'\bU3Matrix');U,nu=matrix(r'\bImageBasis')
assert len(full)==len(B)==len(T)==len(C)==len(H)==25
assert nfull==3*orig and nb==nt==3*n and nc==nh==nu==25
assert B==[a&((1<<(3*n))-1) for a in full] and len(pivots(B))==25
expectedT=[sum(((a>>(3*(3*k-1)))&7)<<(3*(k-1)) for k in range(1,n+1)) for a in full]
expectedC=[sum(sqrt((a>>(3*(2*k-1)))&7)<<(3*(k-1)) for k in range(1,n+1)) for a in full]
assert T==expectedT and mul(H,B)==T
assert mul(C,B)==expectedC
assert mul(C,H)==mul(H,C)
CI=[r^(1<<i) for i,r in enumerate(C)];HI=[r^(1<<i) for i,r in enumerate(H)]
old=mul(C,CI);extra=mul(CI,HI);expected=pivots(old+extra);up=pivots(U)
assert len(pivots(old))==22 and len(pivots(extra))==23
assert len(expected)==len(up)==len(U)==23 and all(contains(expected,r) for r in U)
print({'level':57,'Fourier_coefficients_checked':n,'original_coefficients':orig,'closed_cusp_degree':3,
       'canonical_degree':48,'F2_operator_uniquely_identified':True,'GF8_U3_selection':True,
       'Cartier_Fourier_normalization':True,'U3_commutes_Cartier':True,
       'old_image_dimension':22,'joint_image_dimension':23,'independent_U3_qexp_replay':True})
