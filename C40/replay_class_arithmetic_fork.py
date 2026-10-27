#!/usr/bin/env python3
"""Independent lattice replay of the fresh octic level 40 classification.

Uses the actual exported relation matrix, with no invariant-basis assumption.
Run only after the Magma continuation has completed its export.
"""
import ast
from collections import Counter
import json
import hashlib
from pathlib import Path
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form

HERE = Path(__file__).resolve().parent/'data'
def literals(path):
    result = {}
    for statement in ast.parse(path.read_text()).body:
        assert isinstance(statement, ast.Assign) and len(statement.targets) == 1
        result[statement.targets[0].id] = ast.literal_eval(statement.value)
    return result

data = literals(HERE/'class_arithmetic_direct_RR.py')
assert not data['complete'] and data['interrupted_prefix']
prefix_counts = Counter(row[0] for row in data['tests'])
assert [prefix_counts[d] for d in range(4,8)] == [4,9,10,51]
# The six octic entries in the interrupted prefix are deliberately unused.
data['tests'] = [row for row in data['tests'] if row[0] < 8]
master_sha = hashlib.sha256((HERE/'fixed_octic_inputs.m').read_bytes()).hexdigest()
all_indices = set()
all_witnesses = []
for shard in range(4):
    part = literals(HERE/f'fork_shard{shard}_arithmetic.py')
    assert part['complete'] and part['shard'] == shard and part['master_sha'] == master_sha
    assert part['relations'] == data['relations']
    assert part['rational_cusp_generators'] == data['rational_cusp_generators']
    expected = list(range(shard+1,97,4))
    assert part['expected_indices'] == expected
    assert [row[0] for row in part['tests']] == expected
    witnesses = literals(HERE.parent/'logs'/f'fork_shard{shard}_identities.log.witnesses.py')['witnesses']
    assert len(witnesses) == len(part['tests']) == 24
    for row,witness in zip(part['tests'],witnesses):
        i,vector,member = row
        assert i not in all_indices and vector == witness[2]
        all_indices.add(i)
        data['tests'].append([8,i,vector,member])
        all_witnesses.append((i,witness))
assert all_indices == set(range(1,97))
data['tests'].sort(key=lambda row:(row[0],row[1]))
pairs = literals(HERE/'quartic_pair_classes_fork.py')
assert pairs['complete']
data['quartic_pair_tests'] = pairs['quartic_pair_tests']
geometry = literals(HERE/'factorbase_geometry.py')
proof_state = literals(HERE/'diagnostics/proven_relation_state.py')
assert proof_state['proven'] and proof_state['order'] == 31579088486400
assert geometry['relations'] == data['relations'] == proof_state['abstract_relations']
assert geometry['factorbase_classes'] == proof_state['factorbase_abstract_classes']
assert geometry['factorbase_degrees'] == list(reversed(proof_state['factorbase_degrees']))
assert len(geometry['factorbase_classes']) == len(geometry['factorbase_ideal_bases']) == 616
assert max(geometry['factorbase_degrees']) == 7

raw = data['relations']
n = len(raw)
assert all(len(row) == n + 1 and row[-1] == 0 for row in raw)
R = Matrix([row[:-1] for row in raw])
assert abs(int(R.det())) == 31579088486400

def finite(vector):
    assert len(vector) == n + 1 and vector[-1] == 0
    return tuple(vector[:-1])

def reducer(matrix):
    H = hermite_normal_form(matrix)
    assert H.rows == H.cols == n
    columns = [tuple(int(H[i,j]) for i in range(n)) for j in range(n)]
    assert all(columns[j][j] > 0 for j in range(n))
    assert all(columns[j][i] == 0 for j in range(n) for i in range(j+1,n))
    def reduce(vector):
        v = list(vector)
        for j in range(n-1,-1,-1):
            quotient = v[j] // columns[j][j]
            if quotient:
                for i in range(j+1):
                    v[i] -= quotient * columns[j][i]
        return tuple(v)
    return reduce, abs(int(H.det()))

key, order = reducer(R.T)
qgens = [finite(v) for v in data['rational_cusp_generators']]
hkey, index = reducer(Matrix(R.tolist() + [list(v) for v in qgens]).T)
assert index == 4592 and order // index == 6876979200
zero = (0,) * n
def inH(vector):
    return hkey(vector) == zero
def add(a,b):
    return key(tuple(x+y for x,y in zip(a,b)))

def effective_classes(weights, generators):
    result = [set() for _ in range(9)]
    result[0].add(zero)
    for weight, generator in zip(weights,generators):
        for degree in range(weight,9):
            result[degree].update(add(a,generator) for a in result[degree-weight])
    return result

cw = data['finite_cusp_degrees']
qw = data['rational_cusp_degrees']
assert Counter(cw) == {1:12,2:6,4:6}
assert Counter(qw) == {1:12,2:2,4:6,8:1}
cc = effective_classes(cw,[finite(v) for v in data['finite_cusp_classes']])
qq = effective_classes(qw,qgens)
assert [len(s) for s in cc] == data['finite_effective_class_counts']
assert [len(s) for s in qq] == data['rational_effective_class_counts']
assert all(inH(v) for v in qq[8])
cuspal_members = {v for v in cc[8] if inH(v)}
unknown_cuspal = cuspal_members - qq[8]
assert len(cuspal_members) == data['finite_cuspal_degree8_members']
assert len(unknown_cuspal) == data['unknown_cuspal_classes']

known = {key(finite(v)) for v in data['known_octic_classes']}
assert len(known) == 4 and all(inH(v) for v in known)
assert data['known_reduction_prime_patterns'] == [[[1,4],[1,4]]]
counts = Counter()
candidates = set()
for degree, index_in_degree, vector, reported in data['tests']:
    assert 4 <= degree <= 8
    counts[degree] += 1
    assert index_in_degree == counts[degree]
    vector = finite(vector)
    admitted = {add(vector,z) for z in cc[8-degree] if inH(add(vector,z))}
    assert len(admitted) == reported
    candidates.update(admitted)
assert [counts[d] for d in range(4,9)] == data['representative_counts']
quartics = [finite(v) for v in data['quartic_classes']]
first_indices = data['quartic_representative_indices']
assert len(quartics) == 32 and len(first_indices) == counts[4] == 4
expected_pairs = {(i,j) for i in first_indices for j in range(1,33)}
observed_pairs = set()
for i,j,vector,reported in data['quartic_pair_tests']:
    assert (i,j) not in observed_pairs
    observed_pairs.add((i,j))
    expected = add(quartics[i-1],quartics[j-1])
    assert key(finite(vector)) == expected
    assert bool(reported) == inH(expected)
    if reported:
        candidates.add(expected)
assert observed_pairs == expected_pairs and len(observed_pairs) == 128
unknown = candidates - known - qq[8]
# Each retained principal-function witness is also replayed in the proven
# abstract relation lattice; the geometric both-order replay is separate.
fb = geometry['factorbase_classes']
fd = geometry['factorbase_degrees']
assert all(len(v) == n+1 and v[-1] == d for v,d in zip(fb,fd))
def sparse_vector(rows):
    result = [0]*(n+1)
    used = set()
    degree = 0
    for j,e in rows:
        assert 1 <= j <= 616 and j not in used and e > 0
        used.add(j)
        degree += e*fd[j-1]
        for k in range(n+1):
            result[k] += e*fb[j-1][k]
    return result,degree
for i,w in all_witnesses:
    trial,degree,vector,erows,rrows,arows,function = w
    assert 1 <= trial <= 1000 and degree == 0
    ev,ed = sparse_vector(erows)
    rv,rd = sparse_vector(rrows)
    av,ad = sparse_vector(arows)
    assert ed == rd == 25 and ad == 1
    actual = [a-b+degree*c for a,b,c in zip(ev,rv,av)]
    assert actual[-1] == 0
    assert key(finite(actual)) == key(finite(vector))
    assert len(function) == len(geometry['field_equation'])-1 == 15
    assert any(any(numerator) for numerator,denominator in function)
    for numerator,denominator in function:
        assert denominator and any(denominator)
        assert all(c in (0,1,2) for c in numerator+denominator)

coef = data['frobenius_coefficients']
assert len(coef) == 51 and coef[-1] == 1
traces,closed = [],[]
for m in range(1,9):
    trace = -m*coef[50-m] - sum(coef[50-j]*traces[m-j-1] for j in range(1,m))
    traces.append(trace)
    number = 3**m+1-trace-sum(d*closed[d-1] for d in range(1,m) if m%d == 0)
    assert number % m == 0
    closed.append(number//m)
noncuspidal = [number-Counter(cw)[d] for d,number in enumerate(closed,1)]
assert noncuspidal == [0,0,0,32,72,80,408,708]
report = {'degree':8,'level':40,'rational_image_index':index,
          'noncuspidal_place_counts':noncuspidal,
          'representative_counts':data['representative_counts'],
          'quartic_pairs':len(observed_pairs),'known_octic_classes':len(known),
          'unknown_cuspal_classes':len(unknown_cuspal),
          'unknown_noncuspidal_classes':len(unknown),'complete':not unknown and not unknown_cuspal,
          'master_sha256':master_sha,'disjoint_shards':4,'exact_class_witnesses':len(all_witnesses)}
lift_file = HERE/'residual_short_cusp_lifts.txt'
if unknown or unknown_cuspal:
    assert lift_file.exists()
    lift = ast.literal_eval(lift_file.read_text().replace('[*','[').replace('*]',']'))
    assert lift[0:2] == [40,8] and len(lift[2]) == len(qgens) == 21
    assert lift[3] == qw and len(lift[5]) == len(lift[6])
    lifted_keys = set()
    for coefficients,vector in zip(lift[5],lift[6]):
        assert len(coefficients) == len(qgens)
        actual = tuple(sum(coefficients[j]*qgens[j][k] for j in range(len(qgens))) for k in range(n))
        assert key(actual) == key(finite(vector))
        lifted_keys.add(key(actual))
    assert lifted_keys == unknown | unknown_cuspal
    assert len(lifted_keys) == len(lift[5])
    rr_file = HERE/'residual_RR_p7.txt'
    if rr_file.exists():
        rr_data = {}
        for st in ast.parse(rr_file.read_text()).body:
            assert isinstance(st,ast.Assign) and len(st.targets)==1
            rr_data[st.targets[0].id] = ast.literal_eval(st.value)
        assert rr_data['prime'] == 7 and rr_data['complete']
        assert [i for i,h0 in rr_data['tests']] == list(range(1,len(lift[5])+1))
        remaining = sum(h0>0 for i,h0 in rr_data['tests'])
        assert remaining == rr_data['remaining']
        report['local_filter_complete'] = report['complete']
        report['second_prime_RR_total'] = len(lift[5])
        report['second_prime_RR_remaining'] = remaining
        report['complete'] = remaining == 0
(HERE/'class_arithmetic_fork_replay.json').write_text(json.dumps(report,indent=2)+'\n')
print('C40_OCTIC_FORK_INDEPENDENT_REPLAY_DONE actual_relations=true '
      f'quartic_pairs=128 unknown_cuspal={len(unknown_cuspal)} '
      f'unknown_noncuspidal={len(unknown)} complete={report["complete"]}')
