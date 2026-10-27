"""Independent integer-Hermite replay of all quartic/last18 obstructions."""
import argparse
from collections import Counter
from pathlib import Path
import runpy
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form

class Quotient:
    def __init__(self,R,hs):
        self.n=len(R[0])-1
        assert len(R)==self.n and all(r[-1]==0 for r in R+hs)
        self.H=hermite_normal_form(Matrix([r[:-1] for r in R+hs]).T)
        assert self.H.shape==(self.n,self.n)
        assert all(self.H[i,j]==0 for i in range(self.n) for j in range(i))
        self.index=abs(int(self.H.det()))
    def key(self,v):
        assert len(v)==self.n+1 and v[-1]==0
        w=list(v[:-1])
        for j in reversed(range(self.n)):
            d=int(self.H[j,j]);assert d>0;q=w[j]//d
            for i in range(j+1):w[i]-=q*int(self.H[i,j])
        return tuple(w)
    def add(self,x,y):return self.key([a+b for a,b in zip(x,y)]+[0])

def main():
    ap=argparse.ArgumentParser();ap.add_argument('run',type=Path);P=ap.parse_args().run
    d=runpy.run_path(str(P/'quartic_filter.py'));assert d['complete']
    Q=Quotient(d['relations'],d['rational_cusp_generators']);assert Q.index==4352
    ds=d['finite_cusp_degrees'];cs=d['finite_cusp_classes'];assert sorted(ds)==[1]*16+[4]*4
    zero=(0,)*Q.n;cc=[set() for _ in range(5)];cc[0]={zero}
    for w,v in zip(ds,cs):
        z=Q.key(v)
        for n in range(w,5):cc[n]|={Q.add(x,z) for x in cc[n-w]}
    classes={};seen=set();count=Counter()
    for part,indices,a,passed,known in d['tests']:
        assert not passed and not known
        sig=(tuple(part),tuple(indices));assert sig not in seen;seen.add(sig)
        assert len(part)==len(indices) and part in [[2],[3],[4],[2,2]]
        if len(part)==1:
            classes[part[0],indices[0]]=a
            x=Q.key([-z for z in a]);assert x not in cc[4-sum(part)]
        else:
            i,j=indices;assert i<=j
            expected=[x+y for x,y in zip(classes[2,i],classes[2,j])]
            assert Q.key(expected)==Q.key(a) and Q.key(a)!=zero
        count[tuple(part)]+=1
    assert count=={(2,):20,(3,):128,(4,):602,(2,2):210}
    for de,n in [(2,20),(3,128),(4,602)]:assert {i for d,i in classes if d==de}==set(range(1,n+1))
    assert {indices for part,indices in seen if part==(2,2)}=={(i,j) for i in range(1,21) for j in range(i,21)}
    print('C230_QUARTIC_INTEGER_PASS tests=960 quotient_index=4352')

    a=runpy.run_path(str(P/'quartic_allcusps.py'));assert a['complete'] and sorted(a['cusp_degrees'])==[1]*16+[4]*4
    seen=set()
    for v,h in a['tests']:
        assert len(v)==20 and all(x>=0 for x in v) and sum(x*d for x,d in zip(v,a['cusp_degrees']))==4 and h==1
        assert tuple(v) not in seen;seen.add(tuple(v))
    assert len(seen)==3880
    def cusp_vectors(j,left,v):
        if j==20:
            if left==0:yield tuple(v)
            return
        for m in range(left//a['cusp_degrees'][j]+1):yield from cusp_vectors(j+1,left-m*a['cusp_degrees'][j],v+[m])
    assert seen==set(cusp_vectors(0,4,[]))
    print('C230_QUARTIC_ALLCUSPS_COVERAGE_PASS divisors=3880 dimension=1')

    b=runpy.run_path(str(P/'three_quotients.py'));assert b['complete'] and not b['remaining']
    Q=Quotient(b['relations'],b['rational_cusp_generators']);assert Q.index==4352
    ts={}
    for idx,desc,v,passed in b['targets']:
        assert idx not in ts and desc[0]==8
        assert bool(passed)==(Q.key(v)==(0,)*Q.n);ts[idx]=(desc,v,passed)
    assert set(ts)==set(range(1,10)) and sum(x[2] for x in ts.values())==1
    assert len(b['source_maps'])==18
    for maps in b['source_maps']:
        assert [r[0] for r in maps]==[1,2,3] and len({r[1] for r in maps})==3
        assert all(passed==ts[ti][2] for _,ti,passed in maps)
        assert not all(passed for _,ti,passed in maps)
    tests=b['candidate_tests'];assert len(tests)==18 and {r[0] for r in tests}==set(range(1,19)) and {r[1] for r in tests}==set(range(1,19))
    assert all(not r[2] for r in tests)
    print('C230_FINAL_INTEGER_PASS sources=18 maps=54 targets=9 remaining=0')

if __name__=='__main__':main()
