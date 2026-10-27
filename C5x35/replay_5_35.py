#!/usr/bin/env python3
"""Independent exact elementary arithmetic for the octic C5 x C35 proof."""
from fractions import Fraction
from math import gcd, lcm, isqrt
from pathlib import Path
import hashlib
import json
import re

HERE = Path(__file__).resolve().parent
WINDING = HERE.parent / "large_levels" / "logs" / "winding175_octic.log"


def phi(n):
    return sum(gcd(n, k) == 1 for k in range(1, n + 1))


def char_order(a, b):
    return lcm(20 // gcd(20, a), 6 // gcd(6, b))


def conductor(a, b):
    f25 = 1 if a == 0 else (5 if a % 5 == 0 else 25)
    return f25 * (1 if b == 0 else 7)


assert pow(2, 20, 25) == 1
assert all(pow(2, k, 25) != 1 for k in range(1, 20))
assert pow(3, 6, 7) == 1
assert all(pow(3, k, 7) != 1 for k in range(1, 6))
# epsilon(2 mod25)=exp(2pi i a/20), epsilon(3 mod7)=exp(2pi i b/6).
even_chars = {(a, b) for a in range(20) for b in range(6) if (a + b) % 2 == 0}
orbit_scalars = [k for k in range(60) if gcd(k, 60) == 1]
unseen = set(even_chars)
orbits = []
while unseen:
    a, b = min(unseen)
    orbit = {((a * k) % 20, (b * k) % 6) for k in orbit_scalars}
    assert orbit <= even_chars
    unseen -= orbit
    orbits.append(orbit)
assert len(orbits) == 12 and sum(map(len, orbits)) == 60

pattern = re.compile(
    r"WINDING N=175 order=(\d+) conductor=(\d+) relative_dimension=(\d+) "
    r"winding_dimension=(\d+) coefficient_degree=(\d+)"
)
rows = [tuple(map(int, m.groups())) for m in pattern.finditer(WINDING.read_text())]
assert len(rows) == 12
assert "WINDING_CHARACTER_PROBE_DONE" in WINDING.read_text()
assert len({r[0] for r in rows}) == 12
by_order = {r[0]: r for r in rows}
for orbit in orbits:
    a, b = min(orbit)
    order = char_order(a, b)
    assert len(orbit) == phi(order)
    row = by_order[order]
    assert row[1] == conductor(a, b) and row[4] == len(orbit)
    assert {(char_order(x, y), conductor(x, y)) for x, y in orbit} == {(order, row[1])}

index_175 = Fraction(175**2, 2) * Fraction(24, 25) * Fraction(48, 49)
cusp_sum = sum(phi(d) * phi(175 // d) for d in range(1, 176) if 175 % d == 0)
genus_175 = 1 + index_175 / 12 - Fraction(cusp_sum, 4)
assert index_175 == 14400 and genus_175 == 1033
assert sum(dim * coeff for _, _, dim, _, coeff in rows) == genus_175
bad_orders = [order for order, _, dim, wind, _ in rows if dim != wind]
assert set(bad_orders) == {1, 5}

H = {pow(36, k, 175) for k in range(5)}
assert H == {1, 36, 71, 106, 141}
pm_H = H | {(-h) % 175 for h in H}
assert pow(24, 3, 175) == 174
assert 24 not in pm_H and pow(24, 2, 175) not in pm_H
assert 24 % 35 not in {1, 34} and pow(24, 2, 35) not in {1, 34}

# 36 == 2^16 mod25 and 36 == 1 mod7, hence H-fixed iff5|a.
assert pow(2, 16, 25) == 36 % 25 and 36 % 7 == 1
hfixed = {(a, b) for a, b in even_chars if 16 * a % 20 == 0}
assert len(hfixed) == 12 and all(35 % conductor(a, b) == 0 for a, b in hfixed)
# 24 == -1 mod25 and 24 == 3 mod7.
delta_value = lambda a, b: (Fraction(a, 2) + Fraction(b, 6)) % 1
moving = {(a, b) for a, b in hfixed if delta_value(a, b)}
assert len(moving) == 8
assert {char_order(a, b) for a, b in moving} == {3, 6, 12}
assert all(by_order[char_order(a, b)][2] == by_order[char_order(a, b)][3]
           for a, b in moving)
# Twisting by a character psi mod5 changes nebentype by psi^2:
# a becomes a+10r mod20 for r=0,1,2,3.
for a, b in moving:
    for r in range(4):
        aa = (a + 10 * r) % 20
        assert (aa, b) in hfixed
        assert delta_value(aa, b) == delta_value(a, b)
        assert by_order[char_order(aa, b)][2] == by_order[char_order(aa, b)][3]

assert next(k for k in range(1, 5) if pow(2, k, 5) == 1) == 4
local_intervals = []
for q in (16, 256):
    bound = isqrt(4 * q)
    interval = (q + 1 - bound, q + 1 + bound)
    assert not any(n % 175 == 0 for n in range(interval[0], interval[1] + 1))
    local_intervals.append([q, *interval])
index_mixed = index_175 / len(H)
gonality_lower = Fraction(325, 32768) * index_mixed
assert index_mixed == 2880 and gonality_lower == Fraction(14625, 512)
assert gonality_lower > 8
hfixed_dimension = sum(by_order[char_order(a, b)][2] for a, b in hfixed)
moving_dimension = sum(by_order[char_order(a, b)][2] for a, b in moving)
assert hfixed_dimension == 169 and moving_dimension == 112

report = {
    "degree": 8, "excluded_group": [5, 35],
    "winding_log_sha256": hashlib.sha256(WINDING.read_bytes()).hexdigest(),
    "even_characters": len(even_chars), "galois_character_orbits": len(orbits),
    "genus_X1_175": int(genus_175), "bad_winding_orders": sorted(bad_orders),
    "H": sorted(H), "H_fixed_character_orders": sorted({char_order(*x) for x in hfixed}),
    "delta_moving_character_orders": sorted({char_order(*x) for x in moving}),
    "H_fixed_dimension": hfixed_dimension, "delta_moving_dimension": moving_dimension,
    "PSL_index_mixed": int(index_mixed), "gonality_lower": str(gonality_lower),
    "Hasse_intervals": local_intervals,
    "quartic_nebentype_twists_checked": len(moving) * 4,
    "boundary": "Elementary arithmetic and recorded winding exhaustion replayed; modular symbols remain the exact Magma stage."
}
(HERE / "data" / "replay_5_35.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
print("FIVE_THIRTY_FIVE_ARITHMETIC_REPLAY_DONE")
