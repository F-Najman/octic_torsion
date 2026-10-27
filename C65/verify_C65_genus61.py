#!/usr/bin/env python3
"""Independent GF(2) replay of the C65 quotient kernel/transport certificate."""
import ast
from collections import Counter, deque
import json
from pathlib import Path
import re
import sys

source=Path(sys.argv[1]).read_text()
def expand_range(match):
    a,b=int(match[1]),int(match[2])
    step=int(match[3]) if match[3] else 1
    assert step>0
    return repr(list(range(a,b+1,step)))
source=re.sub(r'\[\s*(-?\d+)\s*\.\.\s*(-?\d+)(?:\s+by\s+(-?\d+))?\s*\]',expand_range,source)
tate=json.loads(Path(sys.argv[2]).read_text())
def value(name):
    matches=re.findall(rf'\b{name}:=(.*?);',source,re.S)
    assert len(matches)==1,name
    return ast.literal_eval(matches[0])

matrices={}
for name,n,m,seq in re.findall(r'(\w+(?:\[<\d+,\d+>\])?):=Matrix\(GF\(2\),(\d+),(\d+),(\[.*?\])\);',source,re.S):
    n,m=int(n),int(m)
    seq=ast.literal_eval(seq)
    assert len(seq)==n*m and all(x in (0,1) for x in seq) and name not in matrices
    matrices[name]=(m,[sum(seq[i*m+j]<<j for j in range(m)) for i in range(n)])
def get(name):
    return matrices[name][1]
def rank(rows):
    piv={}
    for x in rows:
        while x:
            k=x.bit_length()-1
            if k not in piv:
                piv[k]=x
                break
            x^=piv[k]
    return len(piv)
def multiply(A,B):
    out=[]
    for row in A:
        x=0
        while row:
            bit=row&-row
            x^=B[bit.bit_length()-1]
            row^=bit
        out.append(x)
    return out
def plus(A,B):
    assert len(A)==len(B)
    return [a^b for a,b in zip(A,B)]
def power(A,n):
    out=[1<<i for i in range(len(A))]
    while n:
        if n&1:
            out=multiply(out,A)
        A=multiply(A,A)
        n>>=1
    return out
def space_equal(A,B):
    return rank(A)==rank(B)==rank(A+B)

I=[1<<i for i in range(61)]
C=get('CartierMatrix');M5=get('Diamond5Matrix');M13=get('Diamond13Matrix')
BM=get('ImageBasis')
assert len(BM)==48 and rank(BM)==48
assert power(M5,2)==power(M13,6)==I
assert multiply(M5,M13)==multiply(M13,M5)
assert multiply(M5,C)==multiply(C,M5) and multiply(M13,C)==multiply(C,M13)
joint=multiply(plus(C,I),plus(M5,I))+multiply(plus(C,I),plus(M13,I))+multiply(plus(C,I),plus(power(C,6),I))
assert space_equal(joint,BM)
assert space_equal(multiply(BM,M5),BM) and space_equal(multiply(BM,M13),BM)

ds=value('ds');widths=value('CuspWidths');p5=value('CuspPermutation5');p13=value('CuspPermutation13')
assert len(ds)==len(widths)==len(p5)==len(p13)==22
assert Counter(zip(ds,widths))=={(1,65):12,(2,13):6,(6,1):2,(6,5):2}
assert sorted(p5)==sorted(p13)==list(range(1,23))
assert all(p5[p5[i]-1]==i+1 and p5[p13[i]-1]==p13[p5[i]-1] for i in range(22))
for i in range(22):
    j=i+1
    for _ in range(6):
        j=p13[j-1]
    assert j==i+1
    assert (ds[i],widths[i])==(ds[p5[i]-1],widths[p5[i]-1])==(ds[p13[i]-1],widths[p13[i]-1])

# Independently match the actual cusp action with the complete Tate-label action.
places=tate['target_places'];tp5=tate['diamond27_permutation'];tp13=tate['diamond41_permutation']
for width in [65,13,1,5]:
    actual=[i for i,w in enumerate(widths) if w==width]
    expected=[p['index'] for p in places if p['width']==width]
    mapping={actual[0]:expected[0]};todo=deque(actual[:1])
    while todo:
        i=todo.popleft();j=mapping[i]
        assert ds[i]==places[j]['degree']
        for ap,tp in [(p5,tp5),(p13,tp13)]:
            ni,nj=ap[i]-1,tp[j]
            if ni in mapping:
                assert mapping[ni]==nj
            else:
                mapping[ni]=nj;todo.append(ni)
    assert set(mapping)==set(actual) and set(mapping.values())==set(expected)

allowed=value('AllowedCusps');mults=value('AdmissibleMultiplicities')
assert allowed==[i+1 for i,w in enumerate(widths) if w!=5]
for i in range(22):
    expected=[0] if widths[i]==5 else list(range(0,8//ds[i]+1,2 if ds[i]==2 else 1))
    assert mults[i]==expected

reps=value('SourceRepresentatives');transports=value('DeckTransports')
assert len(reps)==3 and sorted(ds[i-1] for i in reps)==[1,2,6]
assert len(transports)==20 and sorted(t[0] for t in transports)==allowed
deck={(a,b):multiply(power(M5,a),power(M13,b)) for a in range(2) for b in range(6)}
source_count=0;local_count=0
for rep in reps:
    ms=mults[rep-1][1:]
    for m in ms:
        key=f'sourcekernels[<{rep},{m}>]'
        W=get(key);A=get(f'basisjets[<{rep},{m}>]')
        assert rank(W)==len(W) and rank(W)+rank(A)==61
        assert multiply(W,A)==[0]*len(W)
        assert matrices[f'basisjets[<{rep},{m}>]'][0]==m*ds[rep-1]
        source_count+=1
assert source_count==11

for i,rep,a,b in transports:
    M=deck[a,b];Minv=deck[a,(-b)%6]
    assert multiply(M,Minv)==I
    j=rep
    for _ in range(a):
        j=p5[j-1]
    for _ in range(b):
        j=p13[j-1]
    assert j==i
    for m in mults[i-1][1:]:
        W=get(f'kernelspaces[<{i},{m}>]');A=get(f'basisjets[<{i},{m}>]');J=get(f'localjets[<{i},{m}>]')
        SW=get(f'sourcekernels[<{rep},{m}>]');SA=get(f'basisjets[<{rep},{m}>]')
        assert space_equal(W,multiply(SW,Minv))
        assert A==multiply(M,SA) and J==multiply(BM,A)
        assert multiply(W,A)==[0]*len(W) and rank(W)+rank(A)==61
        local_count+=1
assert local_count==110
direct_count=0
for name in matrices:
    if name.startswith('directtargetkernels['):
        suffix=name[len('directtargetkernels'):]
        assert space_equal(get(name),get('kernelspaces'+suffix))
        direct_count+=1
assert direct_count==3

columns={}
for i in allowed:
    for m in mults[i-1][1:]:
        name=f'localjets[<{i},{m}>]';ncols,rows=matrices[name]
        assert ncols==m*ds[i-1] and len(rows)==48
        columns[i,m]=[sum(((row>>j)&1)<<r for r,row in enumerate(rows)) for j in range(ncols)]

tested=0;bad=[]
def check(i,remaining,cols,vector):
    global tested
    if i==22:
        if remaining==0:
            assert len(cols)==8
            tested+=1
            r=rank(cols)
            if r!=8:
                bad.append((r,vector))
        return
    for m in mults[i]:
        cost=m*ds[i]
        if cost<=remaining:
            check(i+1,remaining-cost,cols+(columns[i+1,m] if m else []),vector+[m])
check(0,8,[],[])
assert tested==83949
print({'level':65,'genus':61,'image_dimension':48,'source_kernels':source_count,
       'transported_kernels':local_count,'direct_target_checks':direct_count,
       'Tate_cusp_action_match':True,'admissible_divisors_tested':tested,
       'exceptions':len(bad),'first_exceptions':bad[:3],'independent_replay':True})
assert not bad
