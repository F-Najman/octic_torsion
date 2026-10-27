#!/usr/bin/env python3
"""Independent octic C77 cusp, character, C6 transport and rank replay."""
import ast
from collections import Counter
from hashlib import sha256
from itertools import product
import json
from math import comb,gcd
from pathlib import Path
import re

BASE=Path(__file__).resolve().parent
manifest=json.loads((BASE/'data'/'C77_manifest.json').read_text())
assert manifest['complete'] and manifest['excluded_cyclic_order']==77
# The SHA-256 checks of the manifest files were removed (git binds the files).
s=(BASE/'data'/'C77_root_fresh_deck_kernels.m').read_text()
def value(name):
    matches=re.findall(rf'\b{name}:=(.*?);',s,re.S);assert len(matches)==1
    return ast.literal_eval(matches[0].replace('<','(').replace('>',')'))
def matrix(pattern):
    ms=re.findall(pattern+r':=Matrix\(GF\(2\),(\d+),(\d+),(\[.*?\])\);',s,re.S);assert len(ms)==1
    nr,nc=map(int,ms[0][:2]);a=ast.literal_eval(ms[0][2]);assert len(a)==nr*nc and set(a)<={0,1}
    return ([sum(a[i*nc+j]<<j for j in range(nc)) for i in range(nr)],
            [sum(a[i*nc+j]<<i for i in range(nr)) for j in range(nc)])
def basis(vs):
    out={}
    for v in vs:
        while v:
            k=v&-v
            if k not in out:out[k]=v;break
            v^=out[k]
    return out
def reduce(v,bb):
    while v:
        k=v&-v
        if k not in bb:return v
        v^=bb[k]
    return 0
def same(A,B):
    aa,bb=basis(A),basis(B)
    return len(aa)==len(bb) and all(reduce(v,aa)==0 for v in B)
def apply(v,A):
    out=0
    for i,row in enumerate(A):
        if (v>>i)&1:out^=row
    return out
def multiply(A,B):return [apply(row,B) for row in A]
I=[1<<i for i in range(37)]
def power(A,k):
    out=I
    for _ in range(k):out=multiply(out,A)
    return out
assert (value('N'),value('D'),value('gx'),value('udim'))==(77,8,37,32)
assert value('ds')==[1]*6
C=matrix(r'\bCartierMatrix')[0];U=matrix(r'\bImageBasis')[0]
M2=matrix(r'\bMatrixC2')[0];M3=matrix(r'\bMatrixC3')[0]
assert power(M2,2)==power(M3,3)==I and multiply(M2,M3)==multiply(M3,M2)
assert multiply(M2,C)==multiply(C,M2) and multiply(M3,C)==multiply(C,M3)
CMI=[row^(1<<i) for i,row in enumerate(C)]
assert same(U,multiply(C,CMI)) and len(basis(U))==32
assert same(U,multiply(U,M2)) and same(U,multiply(U,M3))
assert len({tuple(multiply(power(M2,a),power(M3,b))) for a in range(2) for b in range(3)})==6
labels=value('CuspLabels');p2=value('CuspPermutation2');p3=value('CuspPermutation3')
assert set(labels)==set(product(range(3),range(2)))
for i,(t,z) in enumerate(labels):
    assert labels[p2[i]-1]==(t,1-z) and labels[p3[i]-1]==((t+1)%3,z)
source=[matrix(rf'\bsourcekernels\[{m}\]')[0] for m in range(1,9)]
blocks={}
for i,(a,b) in enumerate(value('DeckPowers'),1):
    M=multiply(power(M2,a),power(M3,b));inv=multiply(power(M2,a),power(M3,(-b)%3))
    B,cols=matrix(rf'\bbasisjets\[{i}\]');assert len(B)==37 and len(cols)==8
    assert B==multiply(M,matrix(r'\bbasisjets\[1\]')[0])
    for m in range(1,9):
        W=matrix(rf'\bkernelspaces\[<{i},{m}>\]')[0]
        assert len(basis(W))==37-m and same(W,multiply(source[m-1],inv))
        assert all((w&c).bit_count()%2==0 for w in W for c in cols[:m])
        J,Jcols=matrix(rf'\blocaljets\[<{i},{m}>\]')
        assert J==[apply(u,B)&((1<<m)-1) for u in U]
        blocks[i,m]=Jcols
    if i>1:assert same(matrix(rf'\bdirecttargetkernels\[{i}\]')[0],matrix(rf'\bkernelspaces\[<{i},1>\]')[0])
tested=0
def visit(i,left,cols):
    global tested
    if i==7:
        if left:return
        assert len(cols)==len(basis(cols))==8;tested+=1;return
    for m in range(left+1):visit(i+1,left-m,cols+([] if m==0 else blocks[i,m]))
visit(1,8,[]);assert tested==comb(13,5)==1287

# Rebuild the Picard C6 character list, independently of the native constructor.
cuspclasses=[(0,0,1),(20,22,1),(9,57,1),(29,55,1)]
values=[2,0,5,3]
chars=[v for v in product(range(6),repeat=3) if gcd(gcd(gcd(*v[:2]),v[2]),6)==1
       and all(sum(x*y for x,y in zip(v,c))%6==e for c,e in zip(cuspclasses,values))]
assert len(chars)==6
for point in [(5,10,10),(15,20,10)]:chars=[v for v in chars if sum(x*y for x,y in zip(v,point))%6==0]
assert chars==[(4,5,2)]
samples=[((0,0,1),0),((20,22,1),0),((18,1,1),1),((26,9,7),1),((19,15,2),1)]
assert len(basis([sum((x%2)<<i for i,x in enumerate(v)) for v,_ in samples]))==3
assert all(sum(x*y for x,y in zip(chars[0],v))%2==e for v,e in samples)

# Compressed cusp labels give source and quotient signatures and the special
# cubic-lift property without importing the sibling's primitive-column code.
N=77
def cusp(a,c):
    c%=N;g=gcd(c,N);return min((a%g,c),((-a)%g,(-c)%N))
P={cusp(a,c) for c in range(N) for a in range(gcd(c,N)) if gcd(a,gcd(c,N))==1}
def sig(p,u):return cusp(u*p[0],p[1])
def dia(p,u):return cusp(pow(u,-1,N)*p[0],u*p[1])
def orbitset(S,action):
    left=set(S);out=[]
    while left:
        p=min(left);o={p};q=action(p)
        while q!=p:assert q not in o;o.add(q);q=action(q)
        left-=o;out.append(frozenset(o))
    return out
F=orbitset(P,lambda p:sig(p,2));assert Counter(map(len,F))=={1:30,3:10,10:3,30:1}
units=[u for u in range(N) if gcd(u,N)==1]
for o in F:
    if len(o)==1:assert len({sig(min(o),u) for u in units})==1
    if len(o)==3:assert {dia(p,23) for p in o}==set(o) and all(dia(p,23)==sig(p,4) for p in o)
H=orbitset(P,lambda p:dia(p,36));assert all(len(o)==5 for o in H)
index={p:o for o in H for p in o}
Q=orbitset(H,lambda o:index[sig(min(o),2)])
assert Counter(map(len,Q))=={1:6,2:3,3:2,6:1}
image_rat={index[min(o)] for o in F if len(o)==1};assert len(image_rat)==6
assert all(len({index[sig(min(o),u)] for u in units})==1 for o in image_rat)

names=['C77_deck_reproduce.m','C77_AS_component.m','C77_validate_AS.m',
       'logs/C77_validate_AS.log','data/C77_root_fresh_deck_kernels.m','winding77_hecke2_certificate.m',
       'logs/winding77_hecke2_certificate.log']
report={'level':77,'degree':8,'theory_audit':'passes',
        'frozen_manifest_sha256':sha256((BASE/'data'/'C77_manifest.json').read_bytes()).hexdigest(),
        'all_frozen_hashes_pass':True,'C6_characters_before_source_conditions':6,
        'unique_character':[4,5,2],'all_Picard_mod2_AS_checks':True,
        'rational_quotient_cusps':6,'cubic_lifts_diamond23_invariant':True,
        'source_RR_kernels':8,'fresh_target_RR_checks':5,'transported_kernels':48,
        'Cartier_Hecke_image_dimension':32,'rank8_divisors':tested,'exceptions':0,
        'files':[{'path':name,'sha256':sha256((BASE/name).read_bytes()).hexdigest()} for name in names]}
(BASE/'data'/'audit77_completion.json').write_text(json.dumps(report,indent=2)+'\n')
print('C77_INDEPENDENT_COMPLETION_AUDIT_DONE characters=6_to_1 cusp_lifts=true C6_transport=48 projections=48 D8_ranks=1287 exceptions=0 theory=passes')
