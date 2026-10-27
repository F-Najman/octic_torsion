#!/usr/bin/env python3
"""Exact index arithmetic and complete finite-field log coverage audit."""
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
raw = (HERE / "logs" / "check_local_cyclotomic.log").read_text()
assert not re.search(r"error|interrupt|failed", raw, re.I)
wanted = {(5, 15, 11, 1, d) for d in [1, 2]} | {
    (5, 20, 11, 1, d) for d in [1, 2]} | {
    (3, 27, 2, 2, d) for d in [1, 2, 3, 4]}
found = {}
for match in re.finditer(r"OCTIC_LOCAL_CYCLOTOMIC\s+(.*?)hits\s+(\d+)\s+([\d.]+)",
                         raw, re.S):
    fields = re.match(r"(\d+) (\d+) (\d+) (\d+) (\d+) (\d+) (\d+) (\d+)\s+(.*)",
                      match.group(1), re.S)
    assert fields is not None
    m, n, p, f, d, size, ncurves, ncard = map(int, fields.group(*range(1, 9)))
    key = (m, n, p, f, d)
    assert key not in found and key in wanted and size == p**(f*d)
    assert int(match.group(2)) == 0
    expected_curves = {11: 22, 121: 248, 4: 13, 16: 37, 64: 133, 256: 517}
    assert ncurves == expected_curves[size]
    expected_card = {(5, 20, 11, 1, 2): (2, [10, 10]),
                     (3, 27, 2, 2, 3): (1, [9, 9])}
    if key in expected_card:
        count, invariants = expected_card[key]
        assert ncard == count
        assert list(map(int, re.findall(r"\d+", fields.group(9)))) == invariants
    else:
        assert ncard == 0 and fields.group(9).strip() == "{}"
    found[key] = {"field_size": size, "isomorphism_classes": ncurves,
                  "order_divisible_cases": ncard, "subgroup_hits": 0}
assert set(found) == wanted
assert "OCTIC_LOCAL_CYCLOTOMIC_DONE cases=8 hits=0 complete=true" in raw

indices = []
for m, n, factors, d, factor, expected in [
        (5, 15, [3, 5], 2, 1, 480),
        (5, 20, [2, 5], 2, 1, 720),
        (3, 27, [3], 4, 2, 972)]:
    index = Fraction(m*n*n, 2)
    for ell in factors:
        index *= Fraction(ell*ell-1, ell*ell)
    assert index == expected
    bound = Fraction(325*index, 32768)
    assert bound > factor*d
    indices.append({"group": [m, n], "relative_degree": d, "index": int(index),
                    "gonality_lower_bound": str(bound),
                    "required_threshold": factor*d})
result = {"status": "passed", "cases": [dict(group=list(k[:2]), p=k[2],
              residue_base_degree=k[3], relative_degree=k[4], **v)
              for k, v in sorted(found.items())], "index_checks": indices,
          "log_sha256": hashlib.sha256(raw.encode()).hexdigest(),
          "verification_boundary": "Exact index arithmetic and all finite-field "
              "case coverage are independently checked. Finite elliptic group "
              "enumeration remains the reproducible Magma stage."}
(HERE / "data" / "cyclotomic_local_replay.json").write_text(json.dumps(result, indent=2)+"\n")
print(json.dumps(result, sort_keys=True))
