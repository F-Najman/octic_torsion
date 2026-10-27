#!/usr/bin/env python3
"""Exact actual-relation HNF replay, including every new cubic-cusp coset."""
import ast
from collections import Counter
import hashlib
import importlib.util
import json
from math import prod
from pathlib import Path
import re
from sympy import Matrix, ilcm
from sympy.matrices.normalforms import hermite_normal_form

HERE = Path(__file__).resolve().parent / "data"
fn = HERE / "N26_D8_class_arithmetic.txt"
text = fn.read_text()
s = re.sub(r"\btrue\b", "True", text.replace("[*", "[").replace("*]", "]"))
s = re.sub(r"\bfalse\b", "False", s)
inv, relations, generators, cubic, records, fp, nc, reps = ast.literal_eval(s)
assert inv == [2, 14, 70, 1330, 3990, 11970, 0]
assert len(generators) == 18 and len(cubic) == 6 and len(records) == 239
assert nc == [0, 0, 0, 0, 36, 122, 288, 936] and reps == [6, 26, 48, 159]
assert all(len(v) == len(inv) and v[-1] == 0 for v in relations + generators + cubic + [r[2] for r in records])


def basis(rows):
    H = hermite_normal_form(Matrix([v[:-1] for v in rows]).T)
    assert H.shape == (6, 6) and H.det()
    return H


RG, RH = basis(relations), basis(relations + generators)
g_order, index = abs(int(RG.det())), abs(int(RH.det()))
h_order = g_order // index
assert g_order == prod(inv[:-1]) == 124501550040000
assert h_order == 2 * 14 * 266 * 3990 * 11970 == 355718714400
assert index == 350
Hi = RH.inv()
den = int(ilcm(*[a.q for a in Hi]))
M = [[int(den * a) for a in row] for row in Hi.tolist()]


def inside(v):
    return all(sum(a * b for a, b in zip(row, v[:-1])) % den == 0 for row in M)


assert all(inside(v) for v in generators)
tests, surviving = 0, []
seen = Counter()
for d, i, cl, accepted, recorded_matches in records:
    assert d in range(5, 9)
    seen[d] += 1
    assert i == seen[d]
    matches = [0] if inside(cl) else []
    tests += 1
    if d == 5:
        for j, c in enumerate(cubic, 1):
            if inside([a + b for a, b in zip(cl, c)]):
                matches.append(j)
            tests += 1
    assert matches == recorded_matches and bool(matches) == accepted
    if matches:
        surviving.append([d, i, matches])
assert seen == dict(zip(range(5, 9), reps))
assert tests == 275 and surviving == []

# Independently exhaust closed-point degrees using the archived Frobenius polynomial.
assert len(fp) == 51 and fp[-1] == 1 and sum(fp) == g_order
traces = []
closed = []
for n in range(1, 9):
    tr = -n * fp[50 - n]
    tr -= sum(fp[50 - j] * traces[n - j - 1] for j in range(1, n))
    traces.append(tr)
    points = 3**n + 1 - tr
    points -= sum(d * closed[d - 1] for d in range(1, n) if n % d == 0)
    assert points % n == 0
    closed.append(points // n)
assert closed == [18, 0, 6, 0, 36, 122, 288, 936]
assert closed == [v + (18 if d == 1 else 6 if d == 3 else 0) for d, v in enumerate(nc, 1)]

# The abstract Tate cusp enumeration is independent of the finite-curve class calculation.
spec = importlib.util.spec_from_file_location("tate_cusps", HERE.parents[1] / "common" / "cusp_conditions_relative.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
cusps = mod.compute(2, 26, 3, 8)
assert Counter(cusps["closed_cusp_degrees"]) == {1: 18, 3: 6}
assert cusps["rational_residue_cusps_not_F_rational"] == 0

report = {"degree": 8, "excluded_groups": [[2, 26], [1, 52]],
          "actual_relation_HNF_replay": True, "finite_Jacobian_order": g_order,
          "rational_image_order": h_order, "rational_image_index": index,
          "place_counts": nc, "diamond_representatives": reps,
          "representatives_tested": len(records), "cubic_cusp_cosets": 6,
          "individual_HNF_membership_tests": tests, "survivors": len(surviving),
          "Frobenius_closed_point_counts": closed, "Tate_rational_cusps_verified": 18,
          "arithmetic_sha256": hashlib.sha256(fn.read_bytes()).hexdigest(),
          "boundary": "Finite class maps and geometric hash are reproducible Magma stages; exact group membership and zeta/Tate arithmetic independently replayed."}
(HERE / "replay_N26_D8.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
print("N26_D8_INDEPENDENT_HNF_REPLAY_DONE")
