#!/usr/bin/env python3
"""Independent actual-relation replay, including repeated quartic support."""
import ast
from collections import Counter
from pathlib import Path
import json
import re
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form

HERE = Path(__file__).resolve().parent/'data'
data = {}
for st in ast.parse((HERE/'class_arithmetic.py').read_text()).body:
    assert isinstance(st,ast.Assign) and len(st.targets) == 1
    data[st.targets[0].id] = ast.literal_eval(st.value)
assert data['complete'] is True
raw = data['relations']; n = len(raw)
assert all(len(v) == n+1 and v[-1] == 0 for v in raw)
R = Matrix([v[:-1] for v in raw]); detR = abs(int(R.det()))
assert detR == 667924070400
adjR = R.T.adjugate().tolist()
def finite(v):
    assert len(v) == n+1 and v[-1] == 0
    return v[:-1]
def mul(A,v): return tuple(sum(a*x for a,x in zip(row,v)) for row in A)
def key(v): return tuple(int(x) % detR for x in mul(adjR,v))
def add(a,b): return tuple((x+y) % detR for x,y in zip(a,b))
qc = [finite(v) for v in data['rational_cusp_generators']]
degrees = data['rational_cusp_degrees']
assert len(qc) == len(degrees) == 18 and Counter(degrees) == {1:16,8:2}
L = hermite_normal_form(Matrix(R.tolist()+qc).T)
detL = abs(int(L.det())); adjL = L.adjugate().tolist()
assert detL == 128 and detR//detL == 8760*595680
def inH(v): return all(int(x) % detL == 0 for x in mul(adjL,v))
known = {key(finite(v)) for v in data['known_octic_classes']}
assert len(known) == 7
assert data['known_reduction_prime_patterns'] == [[[1,8]],[[2,4]],[[1,4],[1,4]]]
cc = [{(0,)*n}] + [set() for _ in range(4)]
for v,d in zip(qc,degrees):
    if d != 1: continue
    a = key(v)
    for j in range(1,5): cc[j].update(add(a,b) for b in cc[j-1])
counts = Counter(); passing = Counter(); singles = {}; candidates = set()
for d,i,vector,reported in data['tests']:
    v = finite(vector); actual = inH(v); k = key(v)
    assert actual == bool(reported) and i == counts[d]+1
    counts[d] += 1; passing[d] += actual; singles[d,i] = k
    if actual: candidates.update(add(k,c) for c in cc[8-d])
assert counts == {4:26,5:32,6:80,7:256,8:734}
expected_pairs = {(i,j) for i in range(1,27) for j in range(i,27)}
assert {(i,j) for i,j,_,_ in data['quartic_pair_tests']} == expected_pairs
pair_passing = 0
for i,j,vector,reported in data['quartic_pair_tests']:
    v = finite(vector); actual = inH(v); k = key(v)
    assert k == add(singles[4,i],singles[4,j]) and actual == bool(reported)
    if actual: pair_passing += 1; candidates.add(k)
assert candidates == known
coef = data['frobenius_coefficients']; assert len(coef) == 43 and coef[-1] == 1
traces = []; closed = []
for m in range(1,9):
    tr = -m*coef[42-m]-sum(coef[42-j]*traces[m-j-1] for j in range(1,m))
    traces.append(tr)
    num = 3**m+1-tr-sum(d*closed[d-1] for d in range(1,m) if m % d == 0)
    assert num % m == 0; closed.append(num//m)
assert [a-b for a,b in zip(closed,[16,0,0,0,0,0,0,2])] == [0,0,0,26,32,80,256,734]
log = (HERE.parent/'logs'/'classify_octic_34_continuation.log').read_text()
assert not re.search(r'Runtime error|User error|Assertion failed',log)
assert 'unknown_octic_classes=0 complete=true' in log
report = {'degree':8,'level':34,'complete':True,'rational_image_index':128,
          'single_counts':dict(counts),'single_members':dict(passing),
          'quartic_pair_count':351,'quartic_pair_members':pair_passing,
          'known_octic_classes':7,'unknown':0}
(HERE/'class_arithmetic_replay.json').write_text(json.dumps(report,indent=2)+'\n')
print('C34_OCTIC_INDEPENDENT_REPLAY_DONE actual_relations=true single_tests1128=true quartic_pairs351=true repeated_support=true known_octic_classes7=true unknown0=true')
