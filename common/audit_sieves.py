#!/usr/bin/env python3
"""Strict degree-eight partition coverage and independent place-count audit."""
from pathlib import Path
import hashlib
import json
import re
import sys

ROOT = Path(__file__).resolve().parents[1]


def partitions(n, bound=None):
    if not n:
        yield ()
        return
    for first in range(min(n, bound or n), 0, -1):
        for rest in partitions(n-first, first):
            yield (first,) + rest


def audit(n, p, log_name):
    path = ROOT / f"C{n}" / "logs" / log_name
    source = path.read_text()
    assert not re.search(
        r"Runtime error|User error|System error|Assertion failed|syntax error|"
        r"failed memory request|^ERROR", source, re.I | re.M
    )
    assert f"SUMMARY N={n} M=1 p={p} D=8 shard=0/1 total_survivors=0" in source
    assert source.rstrip().endswith("SIEVE_DONE")
    done = {}
    for match in re.finditer(
        r"TYPE \[([^\]]+)\]: no divisors|"
        r"RESULT type=\[([^\]]+)\] total=(\d+) survivors=(\d+)", source
    ):
        part = tuple(map(int, (match[1] or match[2]).replace(",", " ").split()))
        assert part not in done
        done[part] = (int(match[3] or 0), int(match[4] or 0))
    expected = {a for degree in range(1, 9) for a in partitions(degree)}
    assert len(expected) == 66
    assert set(done) == expected and not any(v[1] for v in done.values())
    match = re.search(r"NONCUSP counts deg 1\.\.8: \[([^\]]+)\]", source)
    counts = list(map(int, match[1].replace(",", " ").split()))
    independent = (ROOT / "common" / "logs" / "independent_counts.log").read_text()
    match = re.search(
        rf"DEGREE8_FAST_COUNT m=1 n={n} p={p} .*?"
        r"closed_place_counts=\[([^\]]+)\]", independent
    )
    assert match
    assert counts == list(map(int, match[1].replace(",", " ").split()))
    result = {
        "degree": 8, "level": n, "prime": p, "partition_types": 66,
        "tested_divisors": sum(v[0] for v in done.values()),
        "nonempty_types": [
            {"type": list(a), "total": done[a][0]}
            for a in sorted(done) if done[a][0]
        ],
        "closed_place_counts": counts, "survivors": 0,
        "independent_counts_agree": True, "complete": True,
        "log": log_name,
        "log_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }
    (ROOT / f"C{n}" / "data" / f"N{n}_audit.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))
    print("DEGREE8_SIEVE_COVERAGE_AUDIT_DONE")


if __name__ == "__main__":
    n = int(sys.argv[1])
    p = int(sys.argv[2])
    audit(n, p, sys.argv[3])
