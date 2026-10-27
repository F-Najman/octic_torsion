#!/usr/bin/env python3
"""Independent column-rank replay of fresh D8 local differential data."""
import ast
from collections import Counter
from hashlib import sha256
from pathlib import Path
import re

D8=Path(__file__).resolve().parents[1]
H=D8/'C45'
s=(H/'data'/'N45_completion_jets.m').read_text()
r=(H/'data'/'N45_reference_RR_geometry_D7.m').read_text()
assert ast.literal_eval(re.search(r'\bD:=(\d+);',s)[1])==8
ds=ast.literal_eval(re.search(r'\bds:=(\[.*?\]);',s)[1])
assert Counter(ds)=={1:12,2:4,4:5,6:2}
def matrices(text,pattern):
    out={}
    for m in re.finditer(pattern+r':=Matrix\(GF\(2\),(\d+),(\d+),(\[.*?\])\);',text,re.S):
        g=m.groups();k=tuple(map(int,g[:-3]));nr,nc=map(int,g[-3:-1]);a=ast.literal_eval(g[-1])
        assert len(a)==nr*nc and set(a)<={0,1} and k not in out
        rows=[sum(a[i*nc+j]<<j for j in range(nc)) for i in range(nr)]
        cols=[sum(a[i*nc+j]<<i for i in range(nr)) for j in range(nc)]
        out[k]=(rows,cols,nr,nc)
    return out
def basis(vs):
    b={}
    for v in vs:
        while v:
            k=v & -v
            if k not in b: b[k]=v;break
            v^=b[k]
    return b
def reduce(v,b):
    while v:
        k=v & -v
        if k not in b:return v
        v^=b[k]
    return 0
CM=matrices(s,r'\bCartierMatrix')[()]
U=matrices(s,r'\bImageBasis')[()]
assert CM[2:]==(41,41) and U[2:]==(39,41)
ub=basis(U[0]);cib=basis([v^(1<<i) for i,v in enumerate(CM[0])])
assert len(ub)==len(cib)==39 and all(reduce(v,ub)==0 for v in cib.values())
A=matrices(s,r'basisjets\[(\d+)\]')
J=matrices(s,r'localjets\[<(\d+),(\d+)>\]')
assert len(A)==23 and len(J)==124
for (i,m),(rows,cols,nr,nc) in J.items():
    assert (nr,nc)==(39,m*ds[i-1])
    expected=[]
    for c in A[i,][1][:nc]:
        expected.append(sum(((v & c).bit_count()&1)<<k for k,v in enumerate(U[0])))
    assert cols==expected
assert matrices(r,r'\bCartierMatrix')==matrices(s,r'\bCartierMatrix')
assert matrices(r,r'\bImageBasis')==matrices(s,r'\bImageBasis')
assert ast.literal_eval(re.search(r'\bds:=(\[.*?\]);',r)[1])==ds
K=matrices(r,r'kernelspaces\[<(\d+),(\d+)>\]')
assert len(K)==103
for (i,m),(rows,cols,nr,nc) in K.items():
    cut=A[i,][1][:m*ds[i-1]]
    assert nc==41 and len(basis(rows))+len(basis(cut))==41
    assert all((v & c).bit_count()%2==0 for v in rows for c in cut)
tested=0
def visit(i,left,columns):
    global tested
    if i==len(ds):
        if left:return
        assert len(columns)==8 and len(basis(columns))==8
        tested+=1;return
    for m in range(left//ds[i]+1):
        extra=[] if m==0 else J[i+1,m][1]
        visit(i+1,left-m*ds[i],columns+extra)
visit(0,8,[])
assert tested==148945
good=(H/'logs'/'N45_p2_good_only.log').read_text()
assert 'NONCUSP counts deg 1..8: [ 0, 0, 0, 0, 0, 0, 0, 60 ]' in good
assert re.findall(r'SMART_CASE_RESULT type=\[ 8 \] index=(\d+) killed=(\w+)',good)==[(str(i),'true') for i in range(1,6)]
assert 'SUMMARY N=45 M=1 p=2 D=8 shard=0/1 total_survivors=0' in good
assert not re.search(r'Runtime error|User error|Assertion failed',good)
print(f'OCTIC45_INDEPENDENT_AUDIT_DONE original_kernel_equalities=103 projection_matrices=124 fresh_D8_rank_tests={tested} exceptions=0 fresh_good8_cases=5 complete=true')
