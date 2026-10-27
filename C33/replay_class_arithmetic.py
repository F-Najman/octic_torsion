#!/usr/bin/env python3
"""Fresh integer-lattice replay of the exhaustive octic X1(33) certificate."""
import ast
from collections import Counter
import json
from pathlib import Path
import re
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form

HERE=Path(__file__).resolve().parent/'data'
data={}
for stmt in ast.parse((HERE/'class_arithmetic.py').read_text()).body:
    assert isinstance(stmt,ast.Assign) and len(stmt.targets)==1
    data[stmt.targets[0].id]=ast.literal_eval(stmt.value)
assert data['complete'] is True
raw=data['relations'];n=len(raw)
assert all(len(row)==n+1 and row[-1]==0 for row in raw)
R=Matrix([row[:-1] for row in raw]);detR=abs(int(R.det()))
assert detR==15*42373650
adjR=R.T.adjugate().tolist()
def finite(v):
    assert len(v)==n+1 and v[-1]==0
    return v[:-1]
def mul(A,v):
    return tuple(sum(a*x for a,x in zip(row,v)) for row in A)
def key(v):
    return tuple(int(x)%detR for x in mul(adjR,v))
qc=[finite(v) for v in data['rational_cusp_generators']]
assert len(qc)==17
L=hermite_normal_form(Matrix(R.tolist()+qc).T)
detL=abs(int(L.det()));adjL=L.adjugate().tolist()
assert detR//detL==5*42373650 and detL==3
def inH(v):
    return all(int(x)%detL==0 for x in mul(adjL,v))

Q={key(finite(v)) for v in data['rational_effective_cusps']}
S7={key(finite(v)) for v in data['known_septic_plus_cusp']}
S8={key(finite(v)) for v in data['known_octic_classes']}
assert len(Q)==62055 and len(S7)==100 and len(S8)==5
assert not(Q&S7 or Q&S8 or S7&S8)
types=Counter();passes=Counter();qpasses=Counter();candidates=set()
for degree,v,reportedH,reportedQ in data['tests']:
    v=finite(v);k=key(v);actualH=inH(v);actualQ=k in Q
    assert bool(reportedH)==actualH and bool(reportedQ)==actualQ
    types[degree]+=1;passes[degree]+=actualH;qpasses[degree]+=actualH and actualQ
    if actualH and not actualQ:
        candidates.add(k)
assert candidates==S7|S8
assert types[5]==2*270 and types[6]==20*60 and types[7]==30*10 and types[8]==25
log=(HERE.parent/'logs'/'classify_octic_33.log').read_text()
summary=re.search(r'C33_D8_ALL_CUSPS divisors=(\d+) classes=(\d+)',log)
assert summary and int(summary[1])==62595 and types[0]==int(summary[2])
assert not re.search(r'Runtime error|User error|Assertion failed',log)
assert 'unknown=0 complete=true' in log

# Independently derive every finite-field place count from the exported
# modular-symbol polynomial by Newton identities and Mobius inversion.
coeff=ast.literal_eval((HERE/'frobenius_coefficients.txt').read_text())
assert len(coeff)==43 and coeff[-1]==1
power=[];all_closed=[]
for m in range(1,9):
    trace=-m*coeff[42-m]-sum(coeff[42-j]*power[m-j-1] for j in range(1,m))
    power.append(trace)
    total=2**m+1-trace-sum(d*all_closed[d-1] for d in range(1,m) if m%d==0)
    assert total%m==0;all_closed.append(total//m)
cusp_counts=[10,5,0,0,2,0,0,0]
noncusp=[a-b for a,b in zip(all_closed,cusp_counts)]
assert noncusp==[0,0,0,0,2,20,30,25]
report={'degree':8,'level':33,'complete':True,'noncuspidal_place_counts':noncusp,
        'tests_by_degree':dict(types),'rational_image_passes':dict(passes),
        'rational_effective_passes':dict(qpasses),'candidate_classes':len(candidates),
        'known_septic_plus_cusp':len(S7),'known_octic_classes':len(S8),'unknown':0}
(HERE/'class_arithmetic_replay.json').write_text(json.dumps(report,indent=2)+'\n')
print('C33_OCTIC_INDEPENDENT_REPLAY_DONE candidates=105 septic_plus_cusp=100 octic=5 unknown=0 actual_relations=true complete=true')
