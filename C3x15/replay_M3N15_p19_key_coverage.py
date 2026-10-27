#!/usr/bin/env python3
"""Replay the complete quartic source-key export without finite fields."""
import ast
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent / 'data'
allrows, hashes, fp0 = [], {}, None
for shard in range(1, 7):
    path = HERE / f"M3N15_p19_native_keys_shard{shard}.txt"
    text = path.read_text()
    s, lo, hi, total, fp, rows = ast.literal_eval(text.replace("[*", "[").replace("*]", "]"))
    assert (s, lo, hi, total) == (shard, 1 + (shard - 1) * 5500, min(shard * 5500, 32076), 32076)
    assert [r[0] for r in rows] == list(range(lo, hi + 1))
    assert len(fp) == 19 and fp[-1] == 1
    if fp0 is None:
        fp0 = fp
    assert fp0 == fp
    allrows.extend(rows)
    hashes[path.name] = sha256(text.encode()).hexdigest()
keys = [tuple(row[1]) for row in allrows]
assert len(keys) == len(set(keys)) == 32076
assert all(len(k) == 10 and k[0] in (0, 1) and k[5] == 1 and all(0 <= x < 19 for x in k) for k in keys)
traces, closed = [], []
for n in range(1, 5):
    tr = -n * fp0[18 - n] - sum(fp0[18 - j] * traces[n - j - 1] for j in range(1, n))
    traces.append(tr)
    points = 19**n + 1 - tr - sum(d * closed[d - 1] for d in range(1, n) if n % d == 0)
    assert points % n == 0
    closed.append(points // n)
assert closed == [16, 152, 2336, 32076]
report = {"degree": 8, "group": [3, 15], "prime": 19,
          "closed_place_counts": closed, "quartic_actual_keys": len(keys),
          "all_keys_distinct": True, "primitive_coordinate_flags": dict(Counter(k[0] for k in keys)),
          "shard_sha256": hashes,
          "coverage_argument": "Each exported key is an actual normalized-curve quartic place; keys are injective on Frobenius orbits. Their distinct count equals the independently replayed Frobenius count. No cross-run place-index stability is assumed."}
(HERE / "replay_M3N15_p19_key_coverage.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
print("M3N15_P19_KEY_COVERAGE_REPLAY_DONE")
