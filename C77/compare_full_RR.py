#!/usr/bin/env python3
"""Compare every intrinsic RR kernel of two C77 exports, in their actual basis."""
import ast
from pathlib import Path
import re
import sys

def parse(path):
 s=Path(path).read_text();assert s.count('print "FI_JET_DATA_READY";')==1
 def matrix(pattern):
  ms=re.findall(pattern+r':=Matrix\(GF\(2\),(\d+),(\d+),(\[.*?\])\);',s,re.S);assert len(ms)==1,pattern
  nr,nc=map(int,ms[0][:2]);es=ast.literal_eval(ms[0][2]);assert len(es)==nr*nc
  return [sum(es[i*nc+j]<<j for j in range(nc)) for i in range(nr)]
 return matrix

def echelon(rows):
 p={}
 for x in rows:
  while x:
   k=x.bit_length()-1
   if k in p:x^=p[k]
   else:p[k]=x;break
 return p

def equal(a,b):
 p=echelon(a);q=echelon(b)
 if len(p)!=len(q):return False
 for x in b:
  while x:
   k=x.bit_length()-1
   if k not in p:return False
   x^=p[k]
 return True

a,b=map(parse,sys.argv[1:3]);assert a(r'\bCartierMatrix')==b(r'\bCartierMatrix')
assert equal(a(r'\bImageBasis'),b(r'\bImageBasis'))
for i in range(1,7):
 for m in range(1,9):
  assert equal(a(rf'\bkernelspaces\[<{i},{m}>\]'),b(rf'\bkernelspaces\[<{i},{m}>\]')),(i,m)
print({'level':77,'full_direct_RR_kernels':48,'transported_RR_kernels':48,'all_intrinsic_kernel_spaces_equal':True,
       'actual_Cartier_matrices_equal':True,'image_spaces_equal':True})
