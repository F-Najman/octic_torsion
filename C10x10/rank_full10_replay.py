#!/usr/bin/env python3
"""Replay exact winding dimensions and the H-fixed character arithmetic."""
from pathlib import Path
import hashlib
import json
from math import gcd
import re

HERE = Path(__file__).resolve().parent
rows = []
for line in (HERE / "logs" / "winding100_octic.log").read_text().splitlines():
    m = re.match(r"WINDING N=100 order=(\d+) conductor=(\d+) relative_dimension=(\d+) winding_dimension=(\d+) coefficient_degree=(\d+)", line)
    if m:
        rows.append(tuple(map(int, m.groups())))
assert len(rows) == 6 and all(r[2] == r[3] for r in rows)
assert sum(r[2] * r[4] for r in rows) == 231
U = {a for a in range(1, 100) if gcd(a, 100) == 1}
H = {pow(11, i, 100) for i in range(10)}
assert len(H) == 10 and len(U) == 40 and 99 not in H
powers = {pow(3, k, 100): k for k in range(20)}
assert len(powers) == 20 and 99 not in powers and powers[89] == 18
assert U == set(powers) | {(-x) % 100 for x in powers}
fixed = [k for k in range(20) if 18 * k % 20 == 0]
assert fixed == [0, 10]
orders = [20 // gcd(20, k) for k in fixed]
assert orders == [1, 2]
dim = sum(r[2] * r[4] for r in rows if r[0] in orders)
assert dim == 13
report = {"degree": 8, "curve": "X(10) over Q(zeta5)",
          "claim": "rank zero; no exclusion claim",
          "full_winding_absolute_dimension": 231, "H_order": 10,
          "H_fixed_even_character_orders": orders, "H_fixed_dimension": dim,
          "local_twist_hypothesis": "H-fixed nebentypes have conductor1 or5, hence tame determinant at5; conductor exponent<=2 forces tame wild-inertia action, and every character-mod5 twist stays within level100.",
          "files": {n: hashlib.sha256(p.read_bytes()).hexdigest()
                    for n, p in [("winding_probe_100.m", HERE.parent / "common" / "winding_probe.m"),
                                 ("winding100_octic.log", HERE / "logs" / "winding100_octic.log")]}}
(HERE / "data" / "rank_full10_replay.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
print("FULL10_RANK_CHARACTER_ARITHMETIC_REPLAY_DONE")
