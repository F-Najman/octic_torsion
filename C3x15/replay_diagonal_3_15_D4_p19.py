#!/usr/bin/env python3
"""Independent integer-HNF replay of the relative-quartic four-prime refinement."""
import ast
from collections import Counter, defaultdict
import hashlib
import json
from math import comb, prod
from pathlib import Path
from sympy import Matrix, ilcm
from sympy.matrices.normalforms import hermite_normal_form

HERE = Path(__file__).resolve().parent / 'data'


def read(name):
    text = (HERE / name).read_text()
    return ast.literal_eval(text.replace("[*", "[").replace("*]", "]")), hashlib.sha256(text.encode()).hexdigest()


def basis(rows, n):
    H = hermite_normal_form(Matrix(rows).T)
    assert H.shape == (n, n) and H.det()
    return H


def residue_function(H):
    Hi = H.inv()
    den = int(ilcm(*[x.q for x in Hi]))
    scaled = [[int(den * x) for x in row] for row in Hi.tolist()]
    assert Matrix(scaled) == den * Hi
    return lambda v: tuple(sum(a * b for a, b in zip(row, v)) % den for row in scaled)


def finite(inv, rel, *vectors):
    assert inv[-1] == 0 and all(n > 0 for n in inv[:-1])
    n = len(inv)
    assert all(len(v) == n and v[-1] == 0 for v in rel)
    for vv in vectors:
        assert all(len(v) == n and v[-1] == 0 for v in vv)
    R = [v[:-1] for v in rel]
    H = basis(R, n - 1)
    order = abs(int(H.det()))
    assert order == prod(inv[:-1])
    return R, order, residue_function(H)


def graph(invA, RA, cuspA, invB, RB, cuspB):
    rA, orderA, keyA = finite(invA, RA, cuspA)
    rB, orderB, keyB = finite(invB, RB, cuspB)
    nA, nB = len(invA) - 1, len(invB) - 1
    a = [v[:-1] for v in cuspA]
    b = [v[:-1] for v in cuspB]
    hA = orderA // abs(int(basis(rA + a, nA).det()))
    hB = orderB // abs(int(basis(rB + b, nB).det()))
    rows = [r + [0] * nB for r in rA] + [[0] * nA + r for r in rB]
    rows += [x + y for x, y in zip(a, b)]
    H = basis(rows, nA + nB)
    hdiag = orderA * orderB // abs(int(H.det()))
    assert hA == hB == hdiag == 1990656
    return keyA, keyB, residue_function(H), orderA, orderB


a, sha2 = read("diagonal_3_15_D4_arithmetic.txt")
assert len(a) == 16
(sigs, ds, inv4, R4, inv7, R7, cusp4, cusp7, invH, RH,
 embed4, embed7, states4, states7, recorded_pairs, counts) = a
assert len(sigs) == len(set(map(tuple, sigs))) == 20
assert Counter(ds) == {1: 16, 4: 4}
assert len(cusp4) == len(cusp7) == 20
assert len(states4) == 5102 and len(states7) == 6264
assert counts == [24 * comb(17, 2), 144 * 16, 396, 24 * 25 // 2]
k4, k7, diag, ord4, ord7 = graph(inv4, R4, cusp4, inv7, R7, cusp7)
assert ord4 == 5971968 and ord7 == 149299200
assert prod(invH) == 1990656 and abs(int(basis(RH, len(invH)).det())) == 1990656
assert len(embed4) == len(embed7) == len(invH)


def linear(coeffs, rows):
    return [sum(c * row[j] for c, row in zip(coeffs, rows)) for j in range(len(rows[0]))]


for r in RH:
    assert not any(k4(linear(r, embed4)[:-1]))
    assert not any(k7(linear(r, embed7)[:-1]))
for x, y in zip(embed4, embed7):
    assert not any(diag(x[:-1] + y[:-1]))

lookup = defaultdict(list)
for i, row in enumerate(states4, 1):
    typ, idx, cl, hcl, pred = row
    assert len(cl) == len(inv4) and len(pred) == len(inv7) and cl[-1] == pred[-1] == 0
    assert (typ, idx) == ((0, i) if i <= 5000 else (4, i - 5000))
    v4, v7 = linear(hcl, embed4), linear(hcl, embed7)
    assert not any(k4([x - 3 * y for x, y in zip(v4, cl)][:-1]))
    assert not any(k7([x - y for x, y in zip(v7, pred)][:-1]))
    assert not any(diag([3 * z for z in cl[:-1]] + pred[:-1]))
    lookup[k7(pred[:-1])].append(i)
assert len({k4(row[2][:-1]) for row in states4}) == 5102
assert len(lookup) == 1730
assert Counter(row[0] for row in states7) == {2: counts[0], 3: counts[1], 4: counts[2], 22: counts[3]}
pairs = []
for j, row in enumerate(states7, 1):
    assert len(row[3]) == len(inv7) and row[3][-1] == 0
    for i in lookup.get(k7([3 * z for z in row[3][:-1]]), []):
        pairs.append([i, j])
assert sorted(pairs) == sorted(recorded_pairs)
assert len(pairs) == 60
assert len({i for i, j in pairs}) == 60 and len({j for i, j in pairs}) == 12
assert all(states4[i - 1][0] == 0 and states7[j - 1][0] == 22 for i, j in pairs)
pred7 = sorted({tuple(states4[i - 1][4]) for i, j in pairs})
assert len(pred7) == 12
report = {"degree": 8, "group": [3, 15], "multiplier": 3,
          "two_prime_sha256": sha2, "cusp_group_order": 1990656,
          "two_prime_states": [5102, 6264], "two_prime_matches": 60,
          "two_prime_predictions": 12, "all_nonempty_F4_support_eliminated": True,
          "independent_HNF_graph_and_transport_replay": True}

p13file = HERE / "diagonal_3_15_D4_p13_arithmetic.txt"
if p13file.exists():
    b, sha3 = read(p13file.name)
    assert len(b) == 16
    (sig13, ds13, inv7b, R7b, inv13, R13, cusp7b, cusp13,
     old_preds, new_preds, states13, matches13, counts13, nc13, reps13, labels13) = b
    assert (sig13, ds13, inv7b, R7b, cusp7b) == (sigs, ds, inv7, R7, cusp7)
    assert sorted(map(tuple, old_preds)) == pred7 and len(new_preds) == 12
    kk7, k13, diag713, _, ord13 = graph(inv7, R7, cusp7, inv13, R13, cusp13)
    assert ord13 == 12230590464
    for x, y in zip(old_preds, new_preds):
        assert not any(diag713(x[:-1] + y[:-1]))
    targets = {k13(v[:-1]): i + 1 for i, v in enumerate(new_preds)}
    assert len(targets) == 12
    assert nc13 == [0, 84, 608, 7050]
    assert reps13[0] == 0 and reps13[1] == 3 and reps13[3] == 155
    assert all(len(labels13[d]) == nc13[d] and set(labels13[d]) == set(range(1, reps13[d] + 1))
               for d in range(1, 4))
    assert counts13 == [reps13[1] * comb(17, 2), reps13[2] * 16, reps13[3], reps13[1] * nc13[1]]
    assert len(states13) == sum(counts13)
    found = []
    for j, row in enumerate(states13, 1):
        assert len(row[3]) == len(inv13) and row[3][-1] == 0
        key = k13([3 * z for z in row[3][:-1]])
        if key in targets:
            found.append([j, targets[key]])
    assert sorted(found) == sorted(matches13)
    report.update({"third_prime_sha256": sha3, "third_prime_states": len(states13),
                   "third_prime_matches": len(found), "full_three_prime_arithmetic_replayed": True})
else:
    report["full_three_prime_arithmetic_replayed"] = False

p19file = HERE / "diagonal_3_15_D4_p19_arithmetic.txt"
b, sha4 = read(p19file.name)
assert len(b) == 16
(sig19, ds19, inv7b, R7b, inv19, R19, cusp7b, cusp19,
 old_preds, new_preds, states19, matches19, counts19, nc19, reps19, labels19) = b
assert (sig19, ds19, inv7b, R7b, cusp7b) == (sigs, ds, inv7, R7, cusp7)
assert sorted(map(tuple, old_preds)) == pred7 and len(new_preds) == 12
kk7, k19, diag719, _, ord19 = graph(inv7, R7, cusp7, inv19, R19, cusp19)
assert ord19 == 247669456896
for x, y in zip(old_preds, new_preds):
    assert not any(diag719(x[:-1] + y[:-1]))
targets = {k19(v[:-1]): i + 1 for i, v in enumerate(new_preds)}
assert len(targets) == 12
assert nc19[0] == reps19[0] == 0
assert all(len(labels19[d]) == nc19[d] and set(labels19[d]) == set(range(1, reps19[d] + 1))
           for d in range(1, 4))
# All136 rational pairs AND8 quadratic cusps occur as degree-two complements.
assert counts19 == [reps19[1] * (comb(17, 2) + 8), reps19[2] * 16,
                    reps19[3], reps19[1] * nc19[1]]
assert len(states19) == sum(counts19)
found = []
for j, row in enumerate(states19, 1):
    assert len(row[3]) == len(inv19) and row[3][-1] == 0
    key = k19([3 * z for z in row[3][:-1]])
    if key in targets:
        found.append([j, targets[key]])
assert sorted(found) == sorted(matches19)
report.update({"fourth_prime_sha256": sha4, "fourth_prime_states": len(states19),
               "fourth_prime_NC_counts": nc19, "fourth_prime_matches": len(found),
               "full_four_prime_arithmetic_replayed": True})

(HERE / "replay_diagonal_3_15_D4_p19.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
print("M3N15_DIAGONAL19_INDEPENDENT_HNF_REPLAY_DONE")
