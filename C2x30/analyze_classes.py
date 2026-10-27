"""Independent finite-group/series classification of the native256 rows."""
import argparse
from collections import Counter,deque
import json
from pathlib import Path
import runpy

def main():
    ap=argparse.ArgumentParser();ap.add_argument('run',type=Path);args=ap.parse_args()
    P=args.run;a=runpy.run_path(str(P/'classes256.py'))
    assert a['complete'] and not a['exclusion_claimed']
    R=a['relations'];n=len(R);iv=[R[i][i] for i in range(n)]
    assert all(row==[iv[i] if j==i else 0 for j in range(n+1)] for i,row in enumerate(R))
    def key(v):
        assert len(v)==n+1 and v[-1]==0
        return tuple(x%m for x,m in zip(v,iv))
    hs=a['cusp_block_classes'];ds=a['cusp_block_degrees'];signatures=a['cusp_block_signatures']
    assert ds==[1,1,4,2,4,1,2,1] and len(hs)==8
    zero=(0,)*n;H={zero};q=deque([zero])
    while q:
        x=q.popleft()
        for h in hs:
            y=tuple((x[i]+h[i])%iv[i] for i in range(n))
            if y not in H:H.add(y);q.append(y)
    assert len(H)==256
    best={};count=0
    def enumerate_cusps(j,left,v):
        nonlocal count
        if j==len(ds):
            if left:return
            count+=1;k=tuple(sum(v[t]*hs[t][i] for t in range(8))%iv[i] for i in range(n))
            b=sum((v[t]//2)*ds[t] for t in range(8));h=max(1,b)
            if k not in best or best[k][0]<h:best[k]=(h,list(v),b)
            return
        for m in range(left//ds[j]+1):enumerate_cusps(j+1,left-m*ds[j],v+[m])
    enumerate_cusps(0,8,[]);assert count==602 and len(best)==221
    cases=a['cases'];assert len(cases)==256
    assert [r[0] for r in cases]==list(range(1,257))
    assert {key(r[1]) for r in cases}==H
    records=[];exceptions=[]
    for idx,g,known,original,h in cases:
        k=key(g);assert len(original)==8
        assert tuple(sum(original[t]*hs[t][i] for t in range(8))%iv[i] for i in range(n))==k
        if known:
            assert all(x>=0 for x in original) and sum(x*d for x,d in zip(original,ds))==8
            baseh,v,b=best[k];assert baseh<=h
            if h==1:route='unique_cuspidal'
            elif h==baseh:
                route='full_pullback' if h==4 else 'pullback_with_fixed_cusps'
                assert sum((x%2)*d for x,d in zip(v,ds))==8-2*b
                if route=='full_pullback':assert b==4
                else:assert b>0 and 8-2*b>0
            else:
                assert h==2 and baseh==1;route='exceptional_pencil';exceptions.append([idx,v,h])
            rec=dict(index=idx,route=route,dimension7=h,base_dimension=baseh,base_degree=b,vector=v)
        else:
            assert k not in best and h==0;rec=dict(index=idx,route='dimension_zero',dimension7=0,vector=original)
        rec['class_vector']=g;records.append(rec)
    counts=Counter(r['route'] for r in records)
    assert counts==dict(dimension_zero=35,unique_cuspidal=116,pullback_with_fixed_cusps=96,full_pullback=4,exceptional_pencil=5)
    expected=[[0,0,0,1,1,1,0,1],[0,0,1,0,0,1,1,1],[0,0,1,1,0,0,1,0],[0,1,0,0,1,1,1,0],[0,1,1,1,0,1,0,0]]
    assert [r[1] for r in exceptions]==expected
    out=dict(classes=256,effective_rational_cusp_divisors=602,effective_cuspidal_classes=221,route_counts=dict(counts),records=records,exceptional_pencils=exceptions)
    (P/'series_analysis.json').write_text(json.dumps(out,indent=2)+'\n')
    entries=','.join('<'+str(i)+','+repr(v)+','+str(h)+'>' for i,v,h in exceptions)
    (P/'exceptional_pencil_inputs.m').write_text('C230ExpectedCuspSignatures := '+repr(signatures)+';\nC230ExpectedCuspDegrees := '+repr(ds)+';\nC230ExceptionalCases := [*'+entries+'*];\n')
    print('C230_CLASSES_INDEPENDENT_PASS classes=256 cuspaldivisors=602 routes='+str(dict(counts)))

if __name__=='__main__':main()
