#!/usr/bin/env python3
"""Independent rational linear-algebra replay of the C65 winding complement."""
import ast
from math import gcd
from pathlib import Path
import re
import sys
from sympy import Matrix

text=Path(sys.argv[1]).read_text()
log=Path(sys.argv[2]).read_text()
def matrix(name):
    m=re.search(rf'{name}:=Matrix\(Rationals\(\),(\d+),(\d+),(\[.*?\])\);',text,re.S)
    assert m
    return Matrix(int(m[1]),int(m[2]),ast.literal_eval(m[3]))
S=matrix('CuspidalBasis');W=matrix('WindingBasis')
T=matrix('Hecke2');B=matrix('AdaptedBasis')
assert S.shape==(5,8) and W.shape==(4,8) and S.rank()==5 and W.rank()==4
assert B.det()!=0 and (B[:4,:]*S).row_join(W).shape==(4,16)
assert B[:4,:]*S==W
M=B*T*B.inv()
assert M[:4,4:]==Matrix.zeros(4,1) and M[4,4]==-1
rows=re.findall(r'OCTIC65_WINDING order=(\d+) conductor=(\d+) relative_dim=(\d+) winding_dim=(\d+) coefficient_degree=(\d+)',log)
assert len(rows)==12
weighted=0;deficits=[]
for order,cond,dim,wind,degree in rows:
    order,cond,dim,wind,degree=map(int,(order,cond,dim,wind,degree))
    phi=sum(gcd(i,order)==1 for i in range(1,order+1))
    assert degree==phi
    weighted+=phi*dim
    if dim!=wind:
        deficits.append((order,cond,dim-wind))
assert weighted==121 and deficits==[(1,1,1)]
assert 'FINITE_IMAGE_T2_PLUS1_AND_SIXTH_POWER_AND_DIAMONDS_CERTIFIED' in log
print({'level':65,'weighted_dimension':121,'nontrivial_character_winding_full':True,
       'trivial_winding_dimension':4,'only_upper_positive_dimension':1,
       'Hecke2_on_complement':-1,'exact_rational_replay':True,
       'all_diamond_quotients_inherit_finite_images':True})
