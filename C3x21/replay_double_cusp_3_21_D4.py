#!/usr/bin/env python3
"""Independent exhaustive GF4 rank replay of the doubled-cusp RR test."""
import ast
from collections import Counter
from hashlib import sha256
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
mode = sys.argv[1] if len(sys.argv) > 1 else ""
assert mode in ("", "dt")
suffix = "_dt" if mode else ""
path = HERE / "data" / ("double_cusp_3_21_D4" + suffix + "_arithmetic.txt")
text = path.read_text()
N, D, g, ds, keys, constraints, recorded_bad, tested, total = ast.literal_eval(
    text.replace("[*", "[").replace("*]", "]"))
assert (N, D, g) == (21, 4, 25)
assert Counter(ds) == {1: 24, 3: 8}
assert len(keys) == len(constraints) == 104 and len(set(map(tuple, keys))) == 104
MULT = [[0, 0, 0, 0], [0, 1, 2, 3], [0, 2, 3, 1], [0, 3, 1, 2]]
INV = [0, 1, 3, 2]


def rank(A):
    A = [list(row) for row in A]
    nr, nc = len(A), len(A[0]) if A else 0
    r = 0
    for j in range(nc):
        pivot = next((i for i in range(r, nr) if A[i][j]), None)
        if pivot is None:
            continue
        A[r], A[pivot] = A[pivot], A[r]
        A[r] = [MULT[INV[A[r][j]]][x] for x in A[r]]
        for i in range(r + 1, nr):
            a = A[i][j]
            if a:
                A[i] = [x ^ MULT[a][y] for x, y in zip(A[i], A[r])]
        r += 1
        if r == min(nr, nc):
            break
    return r


local = {}
for (i, m), B in zip(keys, constraints):
    assert 1 <= i <= len(ds) and 1 <= m <= D // ds[i - 1]
    assert len(B) == g and all(len(row) == 2 * m * ds[i - 1] for row in B)
    assert all(0 <= x < 4 for row in B for x in row)
    local[i, m] = B


def weighted(weights, n):
    if not weights:
        if n == 0:
            yield []
        return
    for m in range(n // weights[0] + 1):
        for tail in weighted(weights[1:], n - m * weights[0]):
            yield [m] + tail


bad, dims, count = [], Counter(), 0
for v in weighted(ds, D):
    A = [[] for _ in range(g)]
    for i, m in enumerate(v, 1):
        if m:
            A = [a + b for a, b in zip(A, local[i, m])]
    assert all(len(row) == 2 * D for row in A)
    ell = 2 * D + 1 - rank(A)
    assert ell >= 1
    dims[ell] += 1
    if ell != 1:
        bad.append(v)
    count += 1
assert count == tested == total == 17742
assert bad == recorded_bad
report = {"degree": 8, "group": [3, 21], "relative_degree": 4,
          "finite_field": 4, "genus": g, "local_constraint_matrices": len(local),
          "cusp_divisors_tested": count, "RR_dimensions": dict(dims),
          "exceptions": len(bad), "independent_GF4_exhaustive_replay": True,
          "method": "direct residue-lift recurrence" if mode else "Hensel completion",
          "arithmetic_sha256": sha256(text.encode()).hexdigest()}
(HERE / "data" / ("replay_double_cusp_3_21_D4" + suffix + ".json")).write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
print("M3N21_D4_DOUBLE_CUSP_INDEPENDENT_REPLAY_DONE")
