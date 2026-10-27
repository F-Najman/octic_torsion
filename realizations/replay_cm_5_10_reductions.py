#!/usr/bin/env python3
"""Independent exhaustive finite-field upper bounds for the octic CM model."""
import json

rows = []
for p, z in [(41, 2), (61, 8)]:
    assert (z**8-z**6+z**4-z*z+1) % p == 0
    i = pow(z, 5, p)
    assert (i*i+1) % p == 0 and (11-2*i) % p != 0
    count = 1 + sum((y*y+(1-i)*x*y-i*y-x*x*x+i*x*x) % p == 0
                    for x in range(p) for y in range(p))
    assert count == 50
    rows.append({"prime": p, "zeta20": z, "i": i, "order": count})
print(json.dumps({"status": "passed", "reductions": rows,
                  "full_torsion_upper_order": 50}, sort_keys=True))
