#!/usr/bin/env python3
"""Independent actual-relation HNF replay of the octic X1(35) filter."""
import ast
from collections import Counter
import json
from pathlib import Path
import re
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form

HERE = Path(__file__).resolve().parent/'data'
data = {}
for st in ast.parse((HERE/'class_arithmetic.py').read_text()).body:
    assert isinstance(st, ast.Assign) and len(st.targets) == 1
    data[st.targets[0].id] = ast.literal_eval(st.value)
assert data['complete'] is True
raw = data['relations']; n = len(raw)
assert all(len(row) == n+1 and row[-1] == 0 for row in raw)
R = Matrix([row[:-1] for row in raw]); detR = abs(int(R.det()))
assert detR == 18446099880
adjR = R.T.adjugate().tolist()
def finite(v):
    assert len(v) == n+1 and v[-1] == 0
    return v[:-1]
def mul(A, v):
    return tuple(sum(a*x for a, x in zip(row, v)) for row in A)
def key(v):
    return tuple(int(x) % detR for x in mul(adjR, v))
qc = [finite(v) for v in data['rational_cusp_generators']]
assert len(qc) == 18
L = hermite_normal_form(Matrix(R.tolist()+qc).T)
detL = abs(int(L.det())); adjL = L.adjugate().tolist()
assert detL == 13 and detR//detL == 13*109148520
def inH(v):
    return all(int(x) % detL == 0 for x in mul(adjL, v))
bad = [finite(v) for v in data['bad_cubic_cusps']]
assert len(bad) == 4 and not any(inH(v) for v in bad)
expected_pairs = {(i,j) for i in range(1,5) for j in range(i,5)}
assert {(i,j) for i,j,_,_ in data['bad_cusp_pairs']} == expected_pairs
zero_pairs = []
for i,j,reported,expected in data['bad_cusp_pairs']:
    actual = inH([a+b for a,b in zip(bad[i-1],bad[j-1])])
    assert actual == bool(reported) == bool(expected)
    if actual: zero_pairs.append((i,j))
assert len(zero_pairs) == 2 and sorted(i for pair in zero_pairs for i in pair) == [1,2,3,4]
known = {key(finite(v)) for v in data['known_octic_classes']}
counts = Counter(); passing = Counter(); candidates = set()
for degree, vector, reported in data['tests']:
    v = finite(vector); actual = inH(v)
    assert actual == bool(reported)
    counts[degree] += 1; passing[degree] += actual
    if degree == 8 and actual: candidates.add(key(v))
assert counts == {6:24,7:24,8:30} and passing[6] == passing[7] == 0
assert candidates == known
coef = data['frobenius_coefficients']; assert len(coef) == 51 and coef[-1] == 1
traces = []; closed = []
for m in range(1,9):
    trace = -m*coef[50-m]-sum(coef[50-j]*traces[m-j-1] for j in range(1,m))
    traces.append(trace)
    number = 2**m+1-trace-sum(d*closed[d-1] for d in range(1,m) if m % d == 0)
    assert number % m == 0; closed.append(number//m)
noncusp = [a-b for a,b in zip(closed,[12,0,4,3,0,0,0,0])]
assert noncusp == [0,0,0,0,0,24,24,30]
log = (HERE.parent/'logs'/'classify_octic_35.log').read_text()
assert not re.search(r'Runtime error|User error|Assertion failed', log)
assert 'cusp_classification=true complete=true' in log
report = {'degree':8,'level':35,'complete':True,'rational_image_index':13,
          'noncuspidal_place_counts':noncusp,'classes_checked':78,
          'bad_cusp_singles_checked':4,'bad_cusp_pairs_checked':10,
          'bad_cusp_zero_pairs':zero_pairs,'sextic_members':0,'septic_members':0,
          'known_octic_classes':len(known),'unknown':0}
(HERE/'class_arithmetic_replay.json').write_text(json.dumps(report,indent=2)+'\n')
print(f'C35_OCTIC_INDEPENDENT_REPLAY_DONE actual_relations=true nc_tests=78 cusp_tests=14 known_octic_classes={len(known)} unknown=0 complete=true')
