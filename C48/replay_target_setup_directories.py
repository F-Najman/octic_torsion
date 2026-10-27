"""Independent integer-lattice/Newton replay of fresh C48 target setup."""
import hashlib
import json
import runpy
from pathlib import Path
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form
from sympy.polys.domains import ZZ

HERE = Path(__file__).resolve().parent / 'data'
import sys
DATA = HERE / (sys.argv[1] if len(sys.argv)>1 else ".")
src = DATA / "target_2_24_arithmetic.py"
a = runpy.run_path(str(src))
assert a["proven"] is True

def closed(coeffs, p=5, n=8):
    aa = coeffs[::-1]
    assert aa[0] == 1
    powers, counts = [0], [0]
    for m in range(1, n + 1):
        powers.append(-m * aa[m] - sum(aa[j] * powers[m-j] for j in range(1, m)))
        pts = p**m + 1 - powers[m]
        rest = pts - sum(d * counts[d] for d in range(1, m) if m % d == 0)
        assert rest % m == 0
        counts.append(rest // m)
    return counts[1:]

for name in ("source", "target"):
    assert closed(a[f"{name}_frobenius_coefficients"]) == a[f"{name}_closed_counts"]

R = Matrix(a["relations"])
n = R.rows
assert R.cols == n + 1 and n == 8
assert R[:, n] == Matrix.zeros(n, 1)
R = R[:, :n]
assert R.det() == 4246732800 == sum(a["target_frobenius_coefficients"])
assert all(v[-1] == 0 and len(v) == n + 1 for v in a["rational_cusp_generators"] + a["finite_cusp_classes"])
M = R.T.row_join(Matrix([v[:n] for v in a["rational_cusp_generators"]]).T)
H = hermite_normal_form(M)
index = abs(int(H.det()))
assert index == 4608 and int(R.det()) // index == 921600
S = smith_normal_form(H, domain=ZZ)
assert [abs(int(S[i, i])) for i in range(n) if abs(S[i, i]) != 1] == [2, 2, 24, 48]

def reduce(v):
    v = list(v)
    for i in reversed(range(n)):
        q = v[i] // int(H[i, i])
        for j in range(i + 1):
            v[j] -= q * int(H[j, i])
    return tuple(v)

def add(v, w):
    return reduce([x + y for x, y in zip(v, w)])

zero = (0,) * n
cc = [{zero}] + [set() for _ in range(8)]
for w, v in zip(a["finite_cusp_degrees"], a["finite_cusp_classes"]):
    v = reduce(v[:n])
    for j in range(w, 9):
        cc[j].update(add(z, v) for z in list(cc[j-w]))
counts = list(map(len, cc))
assert counts == a["complement_counts"] == [1, 5, 23, 75, 178, 378, 661, 1065, 1487]
assert sorted(a["block_degrees"]) == [1] * 8 + [2] * 6 + [4] * 5
assert len(a["block_signatures"]) == 19 == len(set(map(tuple, a["block_signatures"])))
assert a["source_noncuspidal_counts"] == [0, 0, 16, 64, 656, 2712, 10784, 48992]
assert a["target_noncuspidal_counts"] == [0, 0, 16, 96, 672, 2712, 11184, 48760]

result = {
    "arithmetic_sha256": hashlib.sha256(src.read_bytes()).hexdigest(),
    "target_group_order": int(R.det()),
    "rational_cusp_image_order": 921600,
    "finite_quotient_index": index,
    "complement_counts": counts,
    "newton_counts_verified": True,
    "integer_lattice_verified": True,
    "setup_only": True,
    "exclusion_claimed": False,
    "complete": True,
}
(DATA / "target_setup_replay.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
