#!/usr/bin/env python3
"""Independent integer-lattice replay of the octic C54 class filter.

Geometric class maps, source coverage and the semistable model are
separate exact inputs; this script checks all exported group arithmetic.
"""
import ast
from pathlib import Path
import sys
from sympy import Matrix,diag
from sympy.matrices.normalforms import hermite_normal_form

def load(path):
    out={}
    for s in ast.parse(path.read_text()).body:
        assert isinstance(s,ast.Assign) and len(s.targets)==1
        name=s.targets[0].id
        assert name not in out
        out[name]=ast.literal_eval(s.value)
    return out
def replay(base,tests,indices):
    def relation(name):
        rows=base['relations'+name]
        assert len(rows)==3 and all(len(r)==4 and r[-1]==0 for r in rows)
        R=Matrix([r[:-1] for r in rows]);assert R.det()!=0
        return R
    R2,R5=relation('2'),relation('5')
    def finite(v):
        assert len(v)==4 and v[-1]==0
        return v[:-1]
    C2=[finite(v) for v in base['cusp2']];C5=[finite(v) for v in base['cusp5']]
    assert len(C2)==len(C5)==14
    H2=hermite_normal_form(Matrix(R2.tolist()+C2).T)
    H5=hermite_normal_form(Matrix(R5.tolist()+C5).T)
    HJ=hermite_normal_form(Matrix(diag(R2,R5).tolist()+[a+b for a,b in zip(C2,C5)]).T)
    assert abs(R2.det())==472473 and abs(H2.det())==1
    assert abs(H5.det())==6518
    assert abs(R5.det())//abs(H5.det())==472473
    assert abs(R2.det()*R5.det())//abs(HJ.det())==472473
    H5inv,HJinv=H5.inv(),HJ.inv();R2inv,R5inv=R2.T.inv(),R5.T.inv()
    def member(inv,v):
        return all(z.q==1 for z in inv*Matrix(v))
    def residue(v):
        return tuple(z%1 for z in R2inv*Matrix(v))
    candidates={residue(finite(v)) for v in base['candidate2']}
    assert len(candidates)==531
    count=members=passes=0;seen=set()
    for row in tests:
        if indices:
            i,q,m,p,im2,preimage=row
            assert i not in seen;seen.add(i)
        else:
            q,m,p,im2,preimage=row
        q=finite(q);is_member=member(H5inv,q)
        assert int(is_member)==m and p in (0,1)
        passed=False
        if is_member:
            members+=1;im2=finite(im2)
            assert len(preimage)==14 and member(HJinv,im2+q)
            got5=[sum(a*c[j] for a,c in zip(preimage,C5)) for j in range(3)]
            got2=[sum(a*c[j] for a,c in zip(preimage,C2)) for j in range(3)]
            assert member(R5inv,[x-y for x,y in zip(got5,q)])
            assert member(R2inv,[x-y for x,y in zip(got2,im2)])
            passed=residue(im2) in candidates
        else:
            assert im2==preimage==[]
        assert int(passed)==p;passes+=passed;count+=1
    assert passes==0
    return {'tests':count,'members':members,'passes':passes},seen

p=Path(__file__).resolve().parent/'data'
def option(flag,default):
    # optional: --base FILE --tests FILE replay another export in the same format (a rerun of C54_lower_reproduce.m)
    return Path(sys.argv[sys.argv.index(flag)+1]) if flag in sys.argv else p/default
base=load(option('--base','C54_class_arithmetic_base.py'))
low=load(option('--tests','C54_replayed_testclasses.py'));assert low['complete'] is True
summary,seen=replay(base,low['tests'],False)
assert summary['tests']==62301
if '--lower-only' in sys.argv:
    print({'degree':8,'level':54,'lower_and_degree7':summary,
           'pure8_not_yet_replayed':True,'independent_integer_replay':True})
    raise SystemExit(0)
allseen=set();shards=[]
for s in range(4):
    data=load(p/f'C54_good8_classes_shard{s}.py')
    assert data['shard']==s and data['complete'] is True
    result,seen=replay(data,data['tests'],True)
    expected={i for i in range(1,3580) if (i-1)%4==s}
    assert seen==expected and not allseen&seen
    allseen|=seen;shards.append(result)
assert allseen==set(range(1,3580))
print({'degree':8,'level':54,'lower_and_degree7':summary,'pure8_shards':shards,
       'all_tests':summary['tests']+sum(s['tests'] for s in shards),
       'all_passing':0,'independent_integer_replay':True})
