#!/usr/bin/env python3
"""Exact fixed-pairing Tate cusp Galois conditions over Q(zeta_m)."""
from collections import Counter
from math import gcd
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent / "data"


def compute(m, n, p, degree):
    assert n % m == 0 and gcd(p, n) == 1
    residue_degree = 1
    while pow(p, residue_degree, m) != 1:
        residue_degree += 1
    labels = {(a, b, c, d) for a in range(n) for b in range(n)
              if gcd(gcd(a, b), n) == 1 for c in range(m) for d in range(m)
              if (c*b-d*a) % m == 1}

    def orbit(x):
        a, b, c, d = x
        return frozenset(((s*(a+k*b)) % n, s*b % n,
                          s*(c+k*d) % m, s*d % m)
                         for s in (-1, 1) for k in range(n))

    classes = sorted(set(map(orbit, labels)), key=lambda s: sorted(s))
    index = {x: i for i, cs in enumerate(classes) for x in cs}
    gal = [t for t in range(1, n) if gcd(t, n) == 1 and t % m == 1 % m]

    def perm(t):
        ans = []
        for cs in classes:
            images = {index[(t*a % n, b, t*c % m, d)] for a, b, c, d in cs}
            assert len(images) == 1
            ans.append(images.pop())
        assert len(set(ans)) == len(classes)
        return ans

    perms = {t: perm(t) for t in gal}
    frob = perms[pow(p, residue_degree, n)]
    left = set(range(len(classes)))
    cycles = []
    while left:
        i = min(left)
        cy = []
        while i not in cy:
            cy.append(i)
            left.remove(i)
            i = frob[i]
        assert i == cy[0]
        cycles.append(cy)
    allowed = [[t for t in gal if all(perms[t][cy[0]] in cy for cy in cycles
                                      if len(cy) <= k)]
               for k in range(1, degree)]
    return {"m": m, "n": n, "p": p, "residue_degree": residue_degree,
            "relative_divisor_degree": degree, "geometric_cusps": len(classes),
            "closed_cusp_degrees": sorted(map(len, cycles)),
            "allowed_residues_by_complement_degree": allowed,
            "rational_residue_cusps_not_F_rational": sum(
                len(cy) == 1 and any(perms[t][cy[0]] != cy[0] for t in gal)
                for cy in cycles)}


if __name__ == "__main__":
    cases = [(3, 15, 7), (3, 18, 7), (3, 21, 13), (3, 24, 7),
             (3, 27, 7), (4, 16, 5), (4, 20, 3), (4, 24, 5),
             (5, 10, 11), (5, 15, 11), (5, 20, 11), (5, 25, 11),
             (8, 8, 3), (10, 10, 3)]
    phi = {3: 2, 4: 2, 5: 4, 8: 4, 10: 4}
    records = [compute(m, n, p, 8//phi[m]) for m, n, p in cases]
    (HERE / "cusp_conditions_relative.json").write_text(
        json.dumps(records, indent=2) + "\n")
    for row in records:
        print(row["m"], row["n"], row["p"],
              dict(Counter(row["closed_cusp_degrees"])),
              row["allowed_residues_by_complement_degree"],
              "nonF_rational_residue_cusps", row["rational_residue_cusps_not_F_rational"])
    print("OCTIC_RELATIVE_CUSP_CONDITIONS_DONE complete=true")
