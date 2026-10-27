#!/usr/bin/env python3
"""Exact integer-HNF replay of the two genus10 quotient obstructions."""
import ast, hashlib, json
from pathlib import Path
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form
P=Path(__file__).resolve().parent/'data'
def read(name):
    raw=(P/name).read_bytes()
    data=ast.literal_eval(raw.decode().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
    return data,hashlib.sha256(raw).hexdigest()
report={};alltests=[]
for name in ['N28_D8_forgetful_class_arithmetic.txt','N28_D8_isogeny_class_arithmetic.txt']:
    (inv,rels,cusp,ds,tests,fp),sha=read(name)
    assert inv==[4,4,24,936,0] and len(ds)==len(cusp)==15
    assert all(r[-1]==0 for r in rels+cusp)
    R=hermite_normal_form(Matrix([r[:-1] for r in rels]).T)
    B=hermite_normal_form(Matrix([r[:-1] for r in rels+cusp]).T)
    assert R.shape==B.shape==(4,4)
    assert abs(R.det())==359424 and abs(B.det())==4
    decisions=[]
    for n,j,d,e,mult,cl,flag in tests:
        assert d==8 and e*mult==d and e in [4,8] and cl[-1]==0
        w=B.inv()*Matrix(cl[:-1]);inside=all(x.q==1 for x in w)
        assert inside==flag
        decisions.append({'record':n,'map':j,'source_degree':d,'image_degree':e,'multiplicity':mult,'in_rational_image':inside})
    report[name]={'sha256':sha,'Jacobian_order':int(abs(R.det())),
                  'rational_image_order':int(abs(R.det())/abs(B.det())),
                  'image_index':int(abs(B.det())),'tests':decisions}
    alltests+=decisions
assert len(alltests)==6
assert all(x['in_rational_image'] for x in alltests if x['record']==35)
assert all(not x['in_rational_image'] for x in alltests if x['record']==46)
report['remaining_full_Aut24_orbits']=[35]
report['exclusion_proved']=False
(P/'replay_N28_D8_quotients.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2));print('N28_D8_QUOTIENT_INTEGER_HNF_REPLAY_DONE')
