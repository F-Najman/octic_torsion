#!/usr/bin/env python3
"""Independent C49 hybrid geometry/rank replay. Input is ONE frozen export.

With no exports, checks only elementary geometry and the stopping rule.
With an export, requires seventeen completed rational and all six cubic jets.
Blocks from separate runs must never be merged: equal Cartier matrices alone
do not identify the underlying global differential bases.
"""
import argparse
import ast
from collections import Counter
from hashlib import sha256
from itertools import combinations_with_replacement
import json
from math import comb,gcd
from pathlib import Path
import re

ap=argparse.ArgumentParser();ap.add_argument('export',nargs='?',type=Path)
ap.add_argument('--output',type=Path);args=ap.parse_args()
def basis(values):
    result={}
    for v in values:
        while v:
            k=v&-v
            if k not in result:result[k]=v;break
            v^=result[k]
    return result
def reduced(v,b):
    while v:
        k=v&-v
        if k not in b:return v
        v^=b[k]
    return 0
def matrices(text,pattern):
    result={}
    for m in re.finditer(pattern+r':=Matrix\(GF\(2\),(\d+),(\d+),(\[.*?\])\);',text,re.S):
        values=m.groups();key=tuple(map(int,values[:-3]));nr,nc=map(int,values[-3:-1])
        data=ast.literal_eval(values[-1]);assert len(data)==nr*nc and set(data)<={0,1}
        assert key not in result
        rows=[sum(data[i*nc+j]<<j for j in range(nc)) for i in range(nr)]
        cols=[sum(data[i*nc+j]<<i for i in range(nr)) for j in range(nc)]
        result[key]=(rows,cols,nr,nc)
    return result

# The elliptic model and the norm argument need only elementary arithmetic.
a1,a2,a3,a4,a6=1,-1,0,-2,-1;x,y=2,-1
assert y*y+a1*x*y+a3*y==x**3+a2*x*x+a4*x+a6
assert 2*y+a1*x+a3==0
b2=a1*a1+4*a2;b4=2*a4+a1*a3;b6=a3*a3+4*a6
b8=a1*a1*a6+4*a2*a6-a1*a3*a4+a2*a3*a3-a4*a4
delta=-b2*b2*b8-8*b4**3-27*b6*b6+9*b2*b4*b6
assert delta==-343 and 52367710906884085342%4==2
points=[(x,y) for x in range(2) for y in range(2)
        if (y*y+x*y-x**3+x*x+2*x+1)%2==0]
assert points==[(0,1)]

# Independent compressed primitive cusp labels and exact C21 regularity.
N=49
def cusp(a,c):
    c%=N;g=gcd(c,N)
    return min((a%g,c),((-a)%g,(-c)%N))
labels={cusp(a,c) for c in range(N) for a in range(gcd(c,N)) if gcd(a,gcd(c,N))==1}
def sigma(P,u):return cusp(u*P[0],P[1])
def diamond(P,u):return cusp(pow(u,-1,N)*P[0],u*P[1])
def orbits(values,action):
    left=set(values);result=[]
    while left:
        p=min(left);o={p};q=action(p)
        while q!=p:assert q not in o;o.add(q);q=action(q)
        left-=o;result.append(o)
    return result
assert len(labels)==60
qo=orbits(labels,lambda P:sigma(P,3));fo=orbits(labels,lambda P:sigma(P,2))
assert Counter(map(len,qo))=={1:21,6:3,21:1}
assert Counter(map(len,fo))=={1:21,3:6,21:1}
rational={next(iter(o)) for o in qo if len(o)==1}
assert len(orbits(rational,lambda P:diamond(P,3)))==1
assert 4*5<21 and 6*comb(21,5)+comb(7,2)*comb(18,2)==125307
report={'level':49,'degree':8,'elliptic_model_checks':True,
 'rational_torsion_2primary_order':2,'quotient_degree':21,
 'cusp_regular_diamond_action':True,'required_rational_cusps':17,
 'required_cubic_cusps':6,'required_mixed_divisors':125307,
 'rank_certificate_complete':False}

if args.export:
    common=None;selected={};hashes=[]
    for path in [args.export]:
        s=path.read_text();hashes.append({'path':str(path),'sha256':sha256(path.read_bytes()).hexdigest()})
        assert re.search(r'\bN:=49;',s) and re.search(r'\bD:=8;',s)
        ds=ast.literal_eval(re.search(r'\bds:=(\[.*?\]);',s)[1]);assert Counter(ds)=={1:21,3:6}
        CM=matrices(s,r'\bCartierMatrix')[()];U=matrices(s,r'\bImageBasis')[()]
        assert CM[2:]==(69,69) and U[2:]==(68,69)
        ci=basis([v^(1<<i) for i,v in enumerate(CM[0])]);ub=basis(U[0])
        assert len(ci)==len(ub)==68 and all(reduced(v,ub)==0 for v in ci.values())
        identity=(ds,CM,U)
        if common is None:common=identity
        else:assert common==identity,'global differential basis differs'
        A=matrices(s,r'basisjets\[(\d+)\]');J=matrices(s,r'localjets\[<(\d+),(\d+)>\]')
        for (i,),block in A.items():
            need=5 if ds[i-1]==1 else 2
            if not all((i,m) in J for m in range(1,need+1)):continue
            for m in range(1,need+1):
                jet=J[i,m];assert jet[2:]==(68,m*ds[i-1])
                expected=[sum(((v&c).bit_count()&1)<<k for k,v in enumerate(U[0]))
                          for c in block[1][:m*ds[i-1]]]
                assert len(expected)==m*ds[i-1] and jet[1]==expected
            selected.setdefault(i,{m:J[i,m][1] for m in range(1,need+1)})
    ds=common[0];rats=sorted(i for i in selected if ds[i-1]==1)
    cubs=sorted(i for i in selected if ds[i-1]==3)
    assert len(rats)>=17 and len(cubs)==6,(len(rats),len(cubs))
    rats=rats[:17];tested=0
    def rational_columns(multiset):
        counts=Counter(multiset);return sum((selected[i][m] for i,m in counts.items()),[])
    for c in cubs:
        for rational_multiset in combinations_with_replacement(rats,5):
            cols=selected[c][1]+rational_columns(rational_multiset)
            assert len(cols)==len(basis(cols))==8
            tested+=1
    for pair in combinations_with_replacement(cubs,2):
        cubic_columns=rational_columns(pair)
        for rational_multiset in combinations_with_replacement(rats,2):
            cols=cubic_columns+rational_columns(rational_multiset)
            assert len(cols)==len(basis(cols))==8
            tested+=1
    assert tested==125307
    report.update(rank_certificate_complete=True,rank_tests=tested,exceptions=0,
                  rational_indices=rats,cubic_indices=cubs,inputs=hashes)
if args.output:args.output.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,sort_keys=True))
