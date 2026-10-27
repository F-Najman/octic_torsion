#!/usr/bin/env python3
"""Independent C77 exact C6 transport replay; geometric RR is Magma input."""
import ast
from pathlib import Path
import re
import sys
s=Path(sys.argv[1]).read_text()
assert s.count('print "FI_JET_DATA_READY";')==1

def val(name):
    ms=re.findall(rf'\b{name}:=(.*?);',s,re.S)
    assert len(ms)==1,name
    return ast.literal_eval(ms[0].replace('<','(').replace('>',')'))

def mat(pattern):
    ms=re.findall(pattern+r':=Matrix\(GF\(2\),(\d+),(\d+),(\[.*?\])\);',s,re.S)
    assert len(ms)==1,pattern
    nr,nc=map(int,ms[0][:2]);es=ast.literal_eval(ms[0][2]);assert len(es)==nr*nc and set(es)<={0,1}
    return [sum(es[i*nc+j]<<j for j in range(nc)) for i in range(nr)],nc

def echelon(rows):
    p={}
    for x in rows:
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

def equal_span(a,b):
    pa,pb=echelon(a),echelon(b)
    return len(pa)==len(pb) and all(contains(pa,x) for x in b)

def apply(x,rows):
    y=0
    while x:
        b=x&-x;y^=rows[b.bit_length()-1];x^=b
    return y

def mul(a,b):return [apply(x,b) for x in a]
def permul(p,q):return [q[j-1] for j in p]
def power(a,k):
    r=I
    for _ in range(k):r=mul(r,a)
    return r

assert (val('N'),val('D'),val('gx'),val('udim'),val('SourceCusp'))==(77,8,37,32,1)
assert val('ds')==[1]*6
I=[1<<i for i in range(37)]
M2,n=mat(r'\bMatrixC2');M3,n3=mat(r'\bMatrixC3');C,nc=mat(r'\bCartierMatrix');U,nu=mat(r'\bImageBasis')
assert n==n3==nc==nu==37 and len(M2)==len(M3)==len(C)==37
assert mul(M2,M2)==I and power(M3,3)==I and mul(M2,M3)==mul(M3,M2)
assert mul(M2,C)==mul(C,M2) and mul(M3,C)==mul(C,M3)
assert equal_span(U,mul(U,M2)) and equal_span(U,mul(U,M3))
p2,p3=val('CuspPermutation2'),val('CuspPermutation3');labels=val('CuspLabels')
assert set(labels)=={(i,j) for i in range(3) for j in range(2)}
assert sorted(p2)==sorted(p3)==list(range(1,7))
for i,(t,z) in enumerate(labels):
    assert labels[p2[i]-1]==(t,1-z)
    assert labels[p3[i]-1]==((t+1)%3,z)
assert permul(p2,p2)==list(range(1,7))
assert permul(permul(p3,p3),p3)==list(range(1,7))
assert permul(p2,p3)==permul(p3,p2)
group={tuple(I):tuple(range(1,7))};queue=[tuple(I)]
while queue:
    a=queue.pop();p=group[a]
    for b,q in [(M2,p2),(M3,p3)]:
        new=tuple(mul(a,b));nperm=tuple(permul(p,q))
        if new in group:assert group[new]==nperm
        else:group[new]=nperm;queue.append(new)
assert len(group)==6 and {p[0] for p in group.values()}==set(range(1,7))
bs,n=mat(r'\bbasisjets\[1\]');assert len(bs)==37 and n==8
powers=val('DeckPowers');assert len(powers)==6
for i,(a,b) in enumerate(powers,1):
    assert a in range(2) and b in range(3)
    M=mul(power(M2,a),power(M3,b));assert group[tuple(M)][0]==i
    target,n=mat(rf'\bbasisjets\[{i}\]');assert target==mul(M,bs)
    # M^{-1}=M2^a*M3^{-b}.
    inv=mul(power(M2,a),power(M3,(-b)%3));assert mul(M,inv)==I
    for m in range(1,9):
        W,n=mat(rf'\bkernelspaces\[<{i},{m}>\]');S,ns=mat(rf'\bsourcekernels\[{m}\]')
        assert n==ns==37 and len(W)==len(S)==37-m
        assert equal_span(W,mul(S,inv)),(i,m,'kernel transport')
    if i>1:
        Q,n=mat(rf'\bdirecttargetkernels\[{i}\]');W,n=mat(rf'\bkernelspaces\[<{i},1>\]')
        assert equal_span(Q,W),(i,'fresh target m1 check')
print({'level':77,'deck_group_order':6,'cusp_orbit_size':6,'direct_source_RR_kernels':8,
       'independent_target_RR_checks':5,'transported_kernel_equalities':48,
       'Cartier_commutes':True,'image_stable':True,'independent_C6_transport_replay':True})
