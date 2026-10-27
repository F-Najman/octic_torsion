# cands8.py -- independent enumeration of the minimal obstructions for degree 8.
# Groups C_m x C_n (m | n), phi(m) | 8 (Weil pairing), all primes dividing n at most 23 (Khawaja, S(8)).
# A minimal obstruction is a group outside G all of whose proper subgroups lie in G.
# Subgroups of C_m x C_n are exactly the C_a x C_b with a | b, a | m, b | n (enumerated directly).
import json, sys
from sympy import primefactors, totient, divisors

S8 = {2, 3, 5, 7, 11, 13, 17, 19, 23}
G = set()
for n in list(range(1, 29)) + [30, 32, 36]: G.add((1, n))
for k in range(1, 13): G.add((2, 2 * k))
for k in range(1, 5): G.add((3, 3 * k))
for k in range(1, 4): G.add((4, 4 * k))
G |= {(5, 5), (6, 6)}
assert len(G) == 52, len(G)          # Phi^infty(8), Najman-Varivoda
G |= {(1, 33), (1, 34), (1, 35), (1, 39), (1, 40), (1, 42), (5, 10), (3, 15)}
assert len(G) == 60

def subgroups(m, n):
    return {(a, b) for b in divisors(n) for a in divisors(m) if b % a == 0}

# G must be closed under subgroups (sanity check of the statement)
for (m, n) in G:
    assert subgroups(m, n) <= G, (m, n)

def okprimes(n): return all(p in S8 for p in primefactors(n))
ms = [m for m in range(1, 200) if 8 % totient(m) == 0]
print("admissible m:", ms)
NMAX = 3000                                                   # generous; the bound 42*23 = 966 is re-derived below
obs = []
for m in ms:
    for n in range(m, NMAX + 1, m):
        if not okprimes(n) or (m, n) in G: continue
        if all(s in G for s in subgroups(m, n) - {(m, n)}):
            obs.append((m, n))
cyc = sorted(n for m, n in obs if m == 1)
noncyc = sorted((m, n) for m, n in obs if m > 1)
print("minimal obstructions:", len(obs), "cyclic:", len(cyc), "noncyclic:", len(noncyc))
print("cyclic n:", cyc)
print("noncyclic (m,n):", noncyc)
print("max cyclic n:", max(cyc), " max noncyclic n:", max(n for m, n in noncyc))
# compare with the Codex ledger (read-only)
L = json.load(open("common/data/ledger.json"))
ex = L["exclusions"]
print("ledger exclusions:", type(ex).__name__, len(ex))
def parse(x):
    if isinstance(x, dict) and "group" in x and isinstance(x["group"], list): return tuple(int(t) for t in x["group"])
    # ledger entries: try to extract (m,n) from a dict or a string
    if isinstance(x, dict):
        for k in ("group", "name", "id", "label"):
            if k in x: return parse(x[k])
        if "m" in x and "n" in x: return (int(x["m"]), int(x["n"]))
        return None
    s = str(x).replace(" ", "").replace("C", "").replace("x", ",").replace("×", ",").replace("(", "").replace(")", "")
    parts = [p for p in s.split(",") if p]
    if len(parts) == 1: return (1, int(parts[0]))
    if len(parts) == 2: return (int(parts[0]), int(parts[1]))
    return None
keys = list(ex.keys()) if isinstance(ex, dict) else ex
led = set()
for k in keys:
    g = parse(k)
    if g is None and isinstance(ex, dict): g = parse(ex[k])
    if g is None: print("unparsed ledger entry:", str(k)[:100]); continue
    led.add(g)
mine = set(obs)
print("ledger parsed:", len(led))
print("in mine, not in ledger:", sorted(mine - led))
print("in ledger, not in mine:", sorted(led - mine))
print("CANDS8_DONE", "AGREE" if mine == led and len(mine) == 80 else "DISAGREE")
