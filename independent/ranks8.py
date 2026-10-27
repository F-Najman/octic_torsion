# ranks8.py -- independent LMFDB check of the analytic-rank / character / dimension inputs of the
# degree-8 draft.  Data: mf_newforms (weight 2) on devmirror.lmfdb.xyz.  Character values are taken
# from the LMFDB field char_values (generators and exponents) and cross-checked against an independent
# implementation of Conrey characters.  Genera of X_H(N) are recomputed with Sage (GammaH), independent
# of LMFDB, and compared with the LMFDB dimension sums (which certifies completeness of the fetched data).
# Output: one line per claim, CONFIRMED / CONTRADICTED / NOT CHECKABLE.   Closing marker RANKS8_DONE.
import json, sys, itertools
from fractions import Fraction as Fr
from math import gcd
from sage.all import GammaH, Gamma1, Gamma0, euler_phi, divisors as sdiv, kronecker_symbol
import psycopg2

def divisors(n): return [d for d in range(1, n + 1) if n % d == 0]
def sigma0(n): return len(divisors(n))

# ---------------------------------------------------------------- fetch
TOP = [57, 63, 65, 69, 75, 77, 81, 85, 91, 95, 99, 105, 115, 119, 121, 125, 133, 143, 161, 169, 175, 187, 209,
       221, 247, 253, 289, 299, 323, 361, 391, 437, 529,
       52, 56, 60, 64, 72, 84, 54, 80, 96, 100, 45, 48, 36, 28, 38, 44, 46, 49, 50, 51, 55, 33, 34, 35, 39, 40,
       21, 23, 25, 27, 30, 15, 42, 26, 32, 24, 20, 18, 14]
levels = sorted({d for n in TOP for d in divisors(n) if d >= 11})
CACHE = "independent/data/lmfdb_newforms_cache.json"
if "--offline" in sys.argv:
    rows = json.load(open(CACHE))
else:
    conn = psycopg2.connect(host="devmirror.lmfdb.xyz", dbname="lmfdb", user="lmfdb", password="lmfdb", connect_timeout=60)
    cur = conn.cursor()
    cur.execute("""select label, level, dim, char_orbit_label, conrey_index, char_order, char_conductor, char_values,
                   analytic_rank, analytic_rank_proved, traces[1:8], char_parity
                   from mf_newforms where weight=2 and level = ANY(%s)""", (levels,))
    rows = [dict(label=r[0], level=r[1], dim=r[2], orbit=r[3], conrey=r[4], order=r[5], cond=r[6], cv=r[7],
                 rank=r[8], proved=r[9], traces=[int(x) for x in r[10]], parity=r[11]) for r in cur.fetchall()]
    json.dump(rows, open(CACHE, "w"))
print("fetched %d weight-2 newforms at %d levels" % (len(rows), len(levels)))
unk = [r["label"] for r in rows if r["rank"] is None]
print("forms with unknown analytic rank:", unk[:20], len(unk))
notproved_pos = [r["label"] for r in rows if r["rank"] and not r["proved"]]
print("positive-rank forms whose rank is not flagged proved:", notproved_pos)

# ---------------------------------------------------------------- characters
def chi_lmfdb(r, a):
    """chi(a) = exp(2 pi i t) from LMFDB char_values [N, order, gens, exps]; returns t in [0,1)."""
    N, order, gens, exps = r["cv"]
    a %= N
    assert gcd(a, N) == 1
    ords = []
    for g in gens:
        k, y = 1, g % N
        while y != 1: y = y * g % N; k += 1
        ords.append(k)
    for es in itertools.product(*[range(o) for o in ords]):
        y = 1
        for g, e in zip(gens, es): y = y * pow(g, e, N) % N
        if y == a:
            return sum(Fr(e * v, order) for e, v in zip(es, exps)) % 1
    raise ValueError("no discrete log")

def factor(n):
    out, p = {}, 2
    while p * p <= n:
        while n % p == 0: out[p] = out.get(p, 0) + 1; n //= p
        p += 1
    if n > 1: out[n] = out.get(n, 0) + 1
    return out
def primroot(q):
    phi = int(euler_phi(q)); fac = list(factor(phi))
    for g in range(2, q):
        if gcd(g, q) == 1 and all(pow(g, phi // r, q) != 1 for r in fac): return g
def dlog(g, x, q):
    y, k = 1, 0
    while y != x % q: y = y * g % q; k += 1
    return k
def conrey(N, c, a):
    """independent Conrey character chi_N(c, a), t in [0,1)."""
    t = Fr(0)
    for p, e in factor(N).items():
        q = p ** e
        if p != 2:
            g = primroot(q); t += Fr(dlog(g, c, q) * dlog(g, a, q), q - q // p)
        elif e == 2:
            t += Fr((1 - (1 if c % 4 == 1 else -1)) * (1 - (1 if a % 4 == 1 else -1)), 8)
        elif e >= 3:
            ec = 1 if c % 4 == 1 else -1; ea = 1 if a % 4 == 1 else -1
            ac = dlog(5, (ec * c) % q, q); aa = dlog(5, (ea * a) % q, q)
            t += Fr((1 - ec) * (1 - ea), 8) + Fr(ac * aa, 2 ** (e - 2))
    return t % 1

# cross-check the two character evaluations on every fetched form
mism = 0
for r in rows:
    N = r["level"]
    for a in range(1, min(N, 60)):
        if gcd(a, N) == 1 and chi_lmfdb(r, a) != conrey(N, r["conrey"], a): mism += 1
print("character cross-check (LMFDB char_values vs independent Conrey formula):", "OK" if mism == 0 else "MISMATCHES %d" % mism)

def chi(r, a): return chi_lmfdb(r, a % r["level"])
def forms(n): return [r for r in rows if n % r["level"] == 0]
def pos(n): return [r for r in forms(n) if r["rank"] and r["rank"] > 0]
def show(rs): return ", ".join("%s(o%d,c%d,r%d%s)" % (r["label"], r["order"], r["cond"], r["rank"], "" if r["proved"] else ",UNPROVED") for r in rs)


def orbit_dims(n, pred=lambda r: True):
    """group forms of level | n by the kernel of their character mod n (= Galois orbit of characters mod n);
    returns {kernel: [order, conductor, relative dim incl. oldforms, positive-rank relative dim]}"""
    units = [a for a in range(1, n) if gcd(a, n) == 1]
    out = {}
    for r in forms(n):
        if not pred(r): continue
        ker = frozenset(a for a in units if chi(r, a) == 0)
        e = out.setdefault(ker, [r["order"], r["cond"], Fr(0), Fr(0)])
        w = Fr(r["dim"] * sigma0(n // r["level"]), int(euler_phi(r["order"])))
        e[2] += w
        if r["rank"]: e[3] += w
    return out

results = []
def report(tag, ok, claim, evidence):
    status = {True: "CONFIRMED", False: "CONTRADICTED", None: "NOT CHECKABLE"}[ok]
    results.append((tag, status))
    print("[%s] %s: %s\n      evidence: %s" % (status, tag, claim, evidence))

# ---------------------------------------------------------------- dimension / genus completeness
def dimsum(n, pred=lambda r: True):
    return sum(r["dim"] * sigma0(n // r["level"]) for r in forms(n) if pred(r))
bad = []
for n in sorted(set(TOP)):
    g = Gamma1(n).genus()
    if dimsum(n) != g: bad.append((n, dimsum(n), g))
report("data-completeness", not bad, "for every level n used, sum over newforms of level M|n of dim*sigma0(n/M) = g(X1(n))",
       "all %d levels agree with Sage genus" % len(set(TOP)) if not bad else str(bad))
for n, g in [(54, 52), (55, 81), (57, 85), (65, 121), (69, 133), (75, 145), (77, 181), (81, 190), (85, 225), (91, 265),
             (95, 289), (99, 281), (100, 231), (105, 289), (175, 1033)]:
    report("genus X1(%d)" % n, Gamma1(n).genus() == g and dimsum(n) == g, "g(X1(%d)) = %d (weighted winding dimension)" % (n, g),
           "Sage %d, LMFDB dim sum %d" % (Gamma1(n).genus(), dimsum(n)))

def trivial_on(r, H): return all(chi(r, h) == 0 for h in H)
def quot_genus(n, H): return dimsum(n, lambda r: trivial_on(r, H))
for (n, H, g, txt) in [(55, [16], 17, "X1(55)/<16>"), (57, [7], 25, "X1(57)/<7>"), (57, [25], 9, "X1(57)/C9 (C9=<25>)"),
                       (65, [51], 61, "X1(65)/<51>"), (65, [4, 51], 21, "X1(65)/<4,51>"), (77, [36], 37, "X1(77)/<36>"),
                       (63, [8], 49, "X1(63)/<8>"), (75, [26], 73, "X1(75)/<26>"), (39, [16], 9, "X1(39)/<16>"),
                       (100, [11], 13, "X(10)=X1(100)/<11>"),
                       (69, [a for a in range(1, 69) if gcd(a, 69) == 1 and kronecker_symbol(69, a) == 1], 13, "X1(69)/ker eps")]:
    sg = GammaH(n, H).genus(); lg = quot_genus(n, H)
    report("genus " + txt, sg == g and lg == g, "genus %d" % g, "Sage GammaH genus %d, LMFDB invariant-dimension %d" % (sg, lg))
report("genus X0(69)", Gamma0(69).genus() == 7, "g(X0(69)) = 7", str(Gamma0(69).genus()))
for (m, n, g) in [(2, 32, 37), (2, 36, 41), (2, 42, 61), (3, 18, 10), (3, 21, 25), (3, 24, 25), (4, 16, 13), (4, 20, 25),
                  (4, 24, 33), (6, 12, 9), (2, 24, 13), (2, 18, 7), (10, 10, 13)]:
    H = [a for a in range(1, m * n) if gcd(a, m * n) == 1 and (a % n == 1)]
    sg = GammaH(m * n, H).genus()
    report("genus X1(%d,%d)" % (m, n), sg == g, "g = %d" % g, "Sage genus of X_H(%d), H={a = 1 mod %d}: %d" % (m * n, n, sg))
for (n, g) in [(36, 17), (45, 41), (48, 37), (51, 65), (25, 12)]:
    report("genus X1(%d)" % n, Gamma1(n).genus() == g, "g = %d" % g, str(Gamma1(n).genus()))

# ---------------------------------------------------------------- rank zero of J1(n) over Q
for n in [38, 44, 45, 46, 49, 50, 51, 54, 55, 69, 75, 81, 33, 34, 35, 39, 40, 48, 28, 36, 21, 23, 25, 27, 30, 15]:
    p = pos(n)
    report("rank J1(%d)(Q)=0" % n, not p, "no newform of level dividing %d has positive analytic rank" % n,
           "positive-rank forms: " + (show(p) if p else "none"))
# X1(2,N) = X_H(2N), H = {+-1, +-(N+1)} over Q
for N in [14, 18, 24, 26, 28, 30, 32, 36, 42]:
    L = 2 * N; p = pos(L); badf = [r for r in p if chi(r, N + 1) == 0]
    report("rank J1(2,%d)(Q)=0" % N, not badf, "every positive-rank form of level | %d has chi(%d) != 1" % (L, N + 1),
           "positive-rank forms at levels | %d: %s; values chi(%d): %s" % (L, show(p) or "none", N + 1, [str(chi(r, N + 1)) for r in p]))
# X1(m,n) over Q(zeta_m) = X_H(mn), H = {a = +-1 mod n}; the set of newforms with chi trivial on H at level | mn
# is stable under twisting by characters mod m (Shimura: level | lcm(mn, m^2, m*cond(chi)) = mn), so rank zero over
# Q(zeta_m) <=> no positive-rank form of level | mn with chi trivial on H.
for (m, n) in [(3, 15), (3, 18), (3, 21), (3, 24), (3, 27), (4, 16), (4, 20), (4, 24), (5, 15), (5, 20), (6, 12), (10, 10)]:
    L = m * n; H = [a for a in range(1, L) if gcd(a, L) == 1 and a % n == 1]
    p = pos(L); badf = [r for r in p if trivial_on(r, H)]
    report("rank J1(%d,%d)(Q(zeta_%d))=0" % (m, n, m), not badf,
           "no positive-rank form of level | %d with character trivial on H={a=1 mod %d}" % (L, n),
           "positive-rank forms at levels | %d: %s; trivial on H: %s" % (L, show(p) or "none", show(badf) or "none"))

# ---------------------------------------------------------------- kernel claims (<a>-1) J1(n)(Q) finite
def kernel_claim(n, a, extra=""):
    p = pos(n); badf = [r for r in p if chi(r, a) != 0]
    report("(<%d>-1)J1(%d)(Q) finite" % (a, n), not badf, "every positive-rank form of level | %d has chi(%d)=1 %s" % (n, a, extra),
           "positive-rank: %s; chi(%d) values: %s" % (show(p) or "none", a, [str(chi(r, a)) for r in p]))
    return p
TABLE = [(115, 3, 44, 7, 6336, 88), (119, 3, 48, 7, 6912, 96), (121, 2, 55, 8, 7260, 110), (125, 3, 50, 7, 7500, 50),
         (133, 64, 3, 8, 8640, 36), (143, 67, 12, 8, 10080, 120), (161, 3, 66, 7, 12672, 132), (169, 3, 39, 7, 14196, 52),
         (187, 122, 16, 8, 17280, 160), (209, 34, 18, 8, 21600, 180), (221, 3, 48, 7, 24192, 96), (247, 144, 3, 8, 30240, 72),
         (253, 3, 55, 7, 31680, 220), (289, 3, 136, 7, 41616, 136), (299, 3, 33, 7, 44352, 264), (323, 101, 18, 8, 51840, 288),
         (361, 3, 171, 7, 64980, 114), (391, 3, 176, 7, 76032, 352), (437, 3, 198, 7, 95040, 396), (529, 3, 253, 7, 139656, 506)]
def Ind(n):
    I = Fr(n * n, 2)
    for p in factor(n): I *= Fr(p * p - 1, p * p)
    return int(I)
def ordpm(a, n):
    k, y = 1, a % n
    while y not in (1, n - 1): y = y * a % n; k += 1
    return k
for (n, a, o, k3, I, dcm) in TABLE:
    kernel_claim(n, a)
    G = Fr(325 * Ind(n), 2 ** 15)
    k3ok = (k3 == 7) == (a % n in (3, pow(3, -1, n)))
    report("Prop3.34 row %d" % n, ordpm(a, n) == o and Ind(n) == I and 8 * k3 < G and n % 3 != 0 and 8 % o != 0 and k3ok,
           "ord<a> in (Z/n)^*/+-1 = %d, I_n = %d, 8k3 = %d < G_n, 3 does not divide n, ord does not divide 8; k3=7 iff a=3^{+-1} (d_CM=%d not checked here)" % (o, I, 8 * k3, dcm),
           "ord=%d, I=%d, G_n=%.2f, k3 rule %s" % (ordpm(a, n), Ind(n), float(G), "consistent" if k3ok else "inconsistent"))
for n in [115, 161, 253, 299, 391, 437, 529]:
    p = pos(n); badf = [r for r in p if r["order"] != 1]
    report("only trivial-character positive rank at %d" % n, not badf, "full winding on every nontrivial-character component",
           "positive-rank: " + (show(p) or "none"))
# 85, 91, 95, 99, 105, 175
p = kernel_claim(85, 2)
nt = sorted({r["orbit"] for r in p if r["order"] != 1})
report("85: only nontrivial positive-rank orbit is 85.j (order 4)", nt == ["j"] and all(r["order"] == 4 for r in p if r["orbit"] == "j"),
       "", "nontrivial-character positive-rank orbits: %s" % nt)
p = kernel_claim(91, 3)
nt = sorted({(r["orbit"], r["order"]) for r in p if r["order"] != 1})
report("91: nontrivial positive-rank orbits are 91.g, 91.u (orders 3, 6)", nt == [("g", 3), ("u", 6)], "", str(nt))
p = kernel_claim(95, 11)
nt = sorted({(r["orbit"], r["order"]) for r in p if r["order"] != 1})
report("95: only nontrivial positive-rank orbit is 95.l (order 12)", nt == [("l", 12)], "", str(nt))
p = pos(99); report("99: positive rank only for trivial character", all(r["order"] == 1 for r in p), "", show(p))
p = kernel_claim(105, 22)
nt = sorted({(r["orbit"], r["order"], r["cond"]) for r in p if r["order"] != 1})
report("105: unique positive-rank (deficient) nontrivial orbit has order 6, conductor 21", len({x[0] for x in nt}) == 1 and all(x[1] == 6 and x[2] == 21 for x in nt), "", str(nt))
report("105: 22 = 1 mod 21", 22 % 21 == 1, "", "")
p0 = [r for r in pos(105) if r["order"] == 1]
report("J0(105)(Q) rank 0", not p0 and Gamma0(105).genus() == 13 and dimsum(105, lambda r: r["order"] == 1) == 13,
       "trivial-character part has dimension 13 and no positive-rank form", "trivial-char positive-rank: %s; dim %d" % (show(p0) or "none", dimsum(105, lambda r: r["order"] == 1)))
p = kernel_claim(175, 24)
report("175: positive-rank characters have order 1 or 5", all(r["order"] in (1, 5) for r in pos(175)), "", show(pos(175)))
# Lemma 3.77 table: characters mod 175 trivial on H = {1,36,71,106,141}
H175 = [1, 36, 71, 106, 141]
od = orbit_dims(175, lambda r: trivial_on(r, H175))
mine = sorted((v[0], v[1], str(v[2]), str(v[2] - v[3])) for v in od.values())
claimed = sorted([(1, 1, "15", "13"), (2, 5, "14", "14"), (3, 7, "14", "14"), (4, 35, "14", "14"), (6, 35, "14", "14"), (12, 35, "14", "14")])
report("Lemma 3.77 table (175, chars trivial on H)", mine == claimed,
       "orbits (order,cond, dim S2(175,eps), dim Wind): (1,1,15,13),(2,5,14,14),(3,7,14,14),(4,35,14,14),(6,35,14,14),(12,35,14,14)",
       "LMFDB (order, cond, relative dim incl. oldforms, relative dim minus positive-rank part): %s" % mine)
report("Lemma 3.77: 24^3 = -1 mod 175 and 24,24^2 not in +-H", pow(24, 3, 175) == 174 and all(x % 175 not in [h % 175 for h in H175] + [(-h) % 175 for h in H175] for x in (24, 576)), "", "")

# 57, 65, 77 trivial-character analysis
def triv_analysis(n, ap_claims, dimtriv):
    p = pos(n); nt = [r for r in p if r["order"] != 1]; tr = [r for r in p if r["order"] == 1]
    posdim = sum(r["dim"] * sigma0(n // r["level"]) for r in tr)
    ok = not nt and posdim == 1 and Gamma0(n).genus() == dimtriv and dimsum(n, lambda r: r["order"] == 1) == dimtriv
    aps = {}
    if len(tr) == 1:
        for pr, v in ap_claims.items(): aps[pr] = tr[0]["traces"][pr - 1]
        ok = ok and all(aps[pr] == v for pr, v in ap_claims.items())
    report("%d: positive rank only trivial char, 1-dim, dim S2(G0(%d))=%d, a_p=%s" % (n, n, dimtriv, ap_claims), ok, "",
           "positive-rank: %s; positive-rank trivial dimension %d; a_p of that form: %s; genus X0=%d" % (show(p), posdim, aps, Gamma0(n).genus()))
triv_analysis(57, {2: -2, 3: -1}, 5)
triv_analysis(65, {2: -1}, 5)
triv_analysis(77, {2: 0}, 7)
# 63: quotient by <8>
p = pos(63); badf = [r for r in p if chi(r, 8) == 0]
report("J(X1(63)/<8>)(Q) rank 0", not badf, "every positive-rank form of level | 63 has chi(8) != 1", "positive-rank: %s" % show(p))
od = orbit_dims(63, lambda r: chi(r, 8) == 0)
lst = sorted((v[0], v[1], str(v[2])) for v in od.values())
report("63: five character orbits with chi(8)=1, weighted dims sum to 49", len(od) == 5 and sum(v[2] * int(euler_phi(v[0])) for v in od.values()) == 49,
       "", "orbits (order, cond, rel dim): %s" % lst)
# 69: epsilon = (69/.)
report("69: eps(5)=1", kronecker_symbol(69, 5) == 1, "", "(69/5) = %d" % kronecker_symbol(69, 5))
# 81 dims 4 + 18*8 + 6*6 + 2*3
report("81: dimension decomposition 4+18*8+6*6+2*3=190", 4 + 18 * 8 + 6 * 6 + 2 * 3 == 190 == dimsum(81), "", "dimsum %d" % dimsum(81))
od = orbit_dims(81)
lst = sorted((int(euler_phi(v[0])), str(v[2])) for v in od.values() if v[2])
report("81: weighted decomposition 4 + 18*8 + 6*6 + 2*3 by character orbits", lst == sorted([(1, "4"), (18, "8"), (6, "6"), (2, "3")]), "",
       "(phi(order), relative dim) per orbit: %s" % lst)
# 100: Lemma 3.48 table
od = orbit_dims(100)
lst = sorted((v[0], v[1], str(v[2])) for v in od.values() if v[2])
report("Lemma 3.48 table (level 100, six even-character orbits)", sorted((v[0], str(v[2])) for v in od.values() if v[2]) == sorted([(1, "7"), (20, "13"), (10, "12"), (5, "12"), (4, "9"), (2, "6")]),
       "orders 1,20,10,5,4,2 with relative dims 7,13,12,12,9,6", "LMFDB orbits with nonzero S2 (order, cond, rel dim): %s" % lst)
p = pos(100)
report("Lemma 3.48: every newform of level | 100 has nonzero central L-value", not p, "", "positive-rank: " + (show(p) or "none"))
# 125 (Lemma 3.75): S2(Gamma1(25)) and the order 5/10 spaces at 125 trivial on <26>
f25 = [r for r in forms(25)]
report("Lemma 3.75: S2(Gamma1(25)) = chars of order 5 (rel dim 1) and 10 (rel dim 2), coefficient fields of degree 4",
       sorted((r["order"], r["dim"] // int(euler_phi(r["order"])), r["dim"]) for r in f25 if r["level"] == 25) == [(5, 1, 4), (10, 2, 8)] and dimsum(25) == 12,
       "", "level-25 newforms (order, rel dim, dim): %s" % sorted((r["label"], r["order"], r["dim"]) for r in f25))
od = orbit_dims(125, lambda r: chi(r, 26) == 0 and r["order"] in (5, 10))
lst = sorted((v[0], v[1], str(v[2])) for v in od.values())
report("Lemma 3.75: dims 7 and 8 of the order-5/order-10 spaces at 125 trivial on <26>", sorted((v[0], str(v[2])) for v in od.values()) == [(5, "7"), (10, "8")], "",
       "LMFDB orbits (order, cond, rel dim incl. oldforms): %s" % lst)
p = [r for r in pos(125) if r["order"] in (5, 10)]
report("Lemma 3.75: rank J1(25)(Q(zeta5)) = 0 (no positive-rank form at level | 125 with char of order 5 or 10)", not p, "",
       "positive-rank forms at level | 125: %s" % show(pos(125)))
# 55, 54
for n in [54, 55]:
    report("J1(%d)(Q) finite" % n, not pos(n), "", show(pos(n)) or "none")

# ---------------------------------------------------------------- gonality index arithmetic
for (n, num, den) in [(54, 325 * 972, 32768), (57, 325 * 1440, 32768), (77, 14625, 512), (105, 2925, 64), (81, 236925, 8192),
                      (85, 8775, 256), (91, 20475, 512), (99, 43875, 1024), (95, 43875, 1024), (69, 10725, 512), (175, 73125, 512)]:
    report("index bound X1(%d)" % n, Fr(325 * Ind(n), 32768) == Fr(num, den), "325*I/32768 = %s" % Fr(num, den),
           "I=%d, value %s" % (Ind(n), Fr(325 * Ind(n), 32768)))
for (m, n, I) in [(5, 15, 480), (5, 20, 720), (3, 27, 972), (3, 24, 576), (4, 20, 576), (5, 35, 2880), (2, 42, 1152), (10, 10, 360), (4, 24, 768)]:
    L = m * n; H = [a for a in range(1, L) if gcd(a, L) == 1 and a % n in (1, n - 1)]
    sI = int(GammaH(L, H).index())          # -1 in Gamma_H, so this is the PSL2 index
    report("index X1(%d,%d)" % (m, n), sI == I == m * Ind(n), "PSL2-index %d" % I, "Sage %d, m*I(X1(n)) = %d" % (sI, m * Ind(n)))

# ---------------------------------------------------------------- J1(5,15), J1(5,20) over F = Q(zeta5) (coordinator request)
# X1(5,5k)_F = X_H(25k)_F with H = {a = +-1 mod 5k}.  Res_{F/Q} Jac has l-adic representation the sum of rho_f (x) psi over
# newforms f (with old multiplicities) of level M | 25k with chi_f trivial on H, and the four characters psi mod 5.
# Shimura (Prop. 3.64): f (x) psi has level dividing lcm(M, 25, 5*cond(chi_f)) and character chi_f psi^2.
units5 = [1, 2, 3, 4]
def psi_val(t_index, a):  # characters mod 5: psi_t(2) = i^t  (2 generates (Z/5)^*)
    k, y = 0, 1
    while y != a % 5: y = y * 2 % 5; k += 1
    return Fr(t_index * k, 4) % 1
for (m, n) in [(5, 15), (5, 20)]:
    L = m * n; H = [a for a in range(1, L) if gcd(a, L) == 1 and a % n in (1, n - 1)]
    S = [r for r in forms(L) if trivial_on(r, H)]
    def lcm(a, b): return a * b // gcd(a, b)
    bounds = sorted({lcm(lcm(r["level"], 25), 5 * r["cond"]) for r in S})
    # character of the twist is trivial on H: chi(h) * psi(h)^2 with h = +-1 mod 5
    chars_ok = all(all((chi(r, h) + 2 * psi_val(t, h % 5)) % 1 == 0 for h in H) for r in S for t in range(4))
    posS = [r for r in S if r["rank"]]
    allpos = pos(L)
    report("rank J1(%d,%d)(Q(zeta5)) = 0 via twists" % (m, n), all(L % b == 0 for b in bounds) and chars_ok and not posS and not allpos,
           "every twist f(x)psi (f in J_H(%d), psi mod 5) is a newform of level | %d with character trivial on H, and all of these have analytic rank 0" % (L, L),
           "%d newforms f with chi_f trivial on H (labels %s); twist-level bounds lcm(M,25,5cond) = %s (all divide %d); twisted characters trivial on H: %s; "
           "positive-rank forms among them: %s; positive-rank forms at ANY level | %d: %s; all ranks proved: %s"
           % (len(S), sorted(r["label"] for r in S), bounds, L, chars_ok, show(posS) or "none", L, show(allpos) or "none", all(r["proved"] for r in forms(L))))

# not a claim of the draft (it uses J1(25)(Q(zeta5)) instead), recorded for information:
p = [r for r in pos(125) if chi(r, 26) == 0]
print("[INFO] X1(5,25)=X_H(125), H={a=+-1 mod 25}: positive-rank forms trivial on H: %s -> J1(5,25)(Q(zeta5)) has positive rank" % show(p))
print("\nSUMMARY:", {s: sum(1 for _, x in results if x == s) for s in ("CONFIRMED", "CONTRADICTED", "NOT CHECKABLE")})
for t, s in results:
    if s != "CONFIRMED": print("  ", s, t)
print("RANKS8_DONE")
