# placecounts_pari.py -- second, independent implementation (PARI/GP via cypari2) of the closed-place counts,
# for residue characteristic p >= 5.  Written because Magma 2.30 point counting over F_{5^7} hangs for some j
# (the Magma run of placecounts8.m with q=5 never finished degree 7).
# Every elliptic curve over F_{q} (q = p^e), all j, all twists (j=0, 1728 handled explicitly; mass formula checked),
# exact group structure by ellgroup; counts of level structures as in placecounts8.m:
#   X1(n): #{P exact order n}/#Aut;  X1(m,n), m>1 (fixed pairing; p^e = 1 mod m): m*#{Q exact order n}/#Aut if E[m] rational.
# Closed places by Moebius inversion.  Extra modes:
#   --x27   : X1(27)/X1(54) analysis over F_5 (k = number of rational nonzero 2-torsion points over the residue field,
#             exact degrees via Tate normal form, diamond orbits) for Prop. 3.13;
#   --orbits: orbits of X1(48) under diamonds and of X1(2,36), X1(2,42) under the automorphism groups of Prop. 3.59/3.62.
# Usage: python3 placecounts_pari.py p D spec1 spec2 ... [--x27] [--orbits]     spec = m:n
import sys, itertools, time
from math import gcd
from cypari2 import Pari
pari = Pari(); pari.allocatemem(4 * 10**9)
pari.set_real_precision(38)

args = [a for a in sys.argv[1:] if not a.startswith("--")]
p = int(args[0]); D = int(args[1]); specs = [tuple(int(t) for t in s.split(":")) for s in args[2:]]
X27 = "--x27" in sys.argv; ORB = "--orbits" in sys.argv
assert p >= 5

def mobius(n):
    r, m, d = 1, n, 2
    while d * d <= m:
        if m % d == 0:
            m //= d
            if m % d == 0: return 0
            r = -r
        d += 1
    return -r if m > 1 else r
def divisors(n): return [d for d in range(1, n + 1) if n % d == 0]
def exact(a, b, n): return sum(mobius(n // e) * gcd(a, e) * gcd(b, e) for e in divisors(n))
def contrib(sp, a, b, aut):
    m, n = sp
    if m == 1: c = exact(a, b, n)
    else:
        if a % m: return 0
        c = m * exact(a, b, n)
    assert c % aut == 0, (sp, a, b, aut)
    return c // aut

# GP helpers -------------------------------------------------------------------------------------------
pari('tatemap(E,P) = my(x0=P[1],y0=P[2],A=E.a4,s,A1,A2,A3); s=(3*x0^2+A)/(2*y0); A1=2*s; A3=2*y0; A2=3*x0-s^2; [-A2^3/A3^2, 1-A1*A2/A3, x0, y0, s, A3/A2]')
pari('tate(E,P) = my(t=tatemap(E,P)); [t[1],t[2]]')
pari('coordsPQ(E,P,Q) = my(t=tatemap(E,Q),X,Y); X=P[1]-t[3]; Y=P[2]-t[4]-t[5]*X; [t[1],t[2],X/t[6]^2,Y/t[6]^3]')
pari('fdeg(v,e,p) = fordiv(e,f, if(vector(#v,i,v[i]^(p^f))==v, return(f))); e')
pari('enc(x,e) = Vecrev(x.pol, e)')
pari('ckey(v,e,p) = my(L=vector(e,i, concat(vector(#v,k, enc(v[k]^(p^(i-1)),e))))); vecsort(L)[1]')
# PARI's ellgroup(E,,1) generators need not have orders equal to the invariants; build a true basis
pari('basisE(E) = my(G=ellgroup(E,0,1),cyc=G[2],gens=G[3],g1,g2,d2,m); if(#cyc<=1, return([cyc,gens])); '
     'g1=gens[1]; g2=gens[2]; if(ellorder(E,g1)!=cyc[1], my(t=g1); g1=g2; g2=t); if(ellorder(E,g1)!=cyc[1], error("no element of exponent order")); '
     'd2=cyc[2]; m=elllog(E,ellmul(E,g2,d2),g1,cyc[1]); if(m%d2, error("bad dlog")); g2=elladd(E,g2,ellmul(E,g1,-m/d2)); '
     'if(ellorder(E,g2)!=d2, error("bad complement")); [cyc,[g1,g2]]')
pari('ptsorder(E,N) = my(B=basisE(E),cyc=B[1],gens=B[2],h=vector(#cyc,i, ellmul(E,gens[i],cyc[i]/gcd(cyc[i],N))),g=vector(#cyc,i,gcd(cyc[i],N)),L=List()); '
     'if(#cyc==1, for(u=0,g[1]-1, my(P=ellmul(E,h[1],u)); if(P!=[0] && ellorder(E,P)==N, listput(L,P))), '
     ' for(u=0,g[1]-1, for(v=0,g[2]-1, my(P=elladd(E,ellmul(E,h[1],u),ellmul(E,h[2],v))); if(P!=[0] && ellorder(E,P)==N, listput(L,P))))); Vec(L)')
fdeg = pari('fdeg'); ckey = pari('ckey'); tate = pari('tate'); coordsPQ = pari('coordsPQ'); ptsorder = pari('ptsorder')

def curves(q, g, j, prim):
    """all twists over F_q with j-invariant j: list of (ellinit, #Aut)."""
    if j == 0:
        k = gcd(6, q - 1)
        return [(pari.ellinit([0, 0, 0, 0, prim**i], g), k) for i in range(k)]
    if j == 1728:
        k = gcd(4, q - 1)
        return [(pari.ellinit([0, 0, 0, prim**i, 0], g), k) for i in range(k)]
    A = 3 * j * (1728 - j); B = 2 * j * (1728 - j)**2
    return [(pari.ellinit([A, B], g), 2), (pari.ellinit([A * prim**2, B * prim**3], g), 2)]

pts = {sp: [0] * (D + 1) for sp in specs}
A27 = {}; ORB27 = {}
orbspecs = [(1, 48, "diamond"), (2, 36, "aut2"), (2, 42, "aut2_42")] if ORB else []
ORBS = {(sp, e): [] for sp in orbspecs for e in range(1, D + 1)}
t0 = time.time()
for e in range(1, D + 1):
    q = p**e
    g = pari('ffgen(%d^%d,\'a)' % (p, e)); prim = pari.ffprimroot(g)
    basis = [g**i for i in range(e)]
    seen27 = set(); orbit27 = []
    seenO = {sp: set() for sp in orbspecs}
    massbad = 0
    for cf in itertools.product(range(p), repeat=e):
        j = sum(c * b for c, b in zip(cf, basis)) if any(cf) else 0 * g
        cl = curves(q, g, j, prim)
        if sum(pari(1) / aut for E, aut in cl) != 1: massbad += 1
        jdeg = None
        for E, aut in cl:
            N = int(pari.ellcard(E))
            rel = [sp for sp in specs if N % (sp[0] * sp[1]) == 0]
            need27 = X27 and N % 27 == 0
            needO = [sp for sp in orbspecs if N % (sp[0] * sp[1]) == 0]
            if not rel and not need27 and not needO: continue
            cyc = [int(x) for x in pari.ellgroup(E)]
            b_ = cyc[0]; a_ = cyc[1] if len(cyc) == 2 else 1
            assert a_ * b_ == N
            for sp in rel: pts[sp][e] += contrib(sp, a_, b_, aut)
            if need27 or needO:
                if jdeg is None: jdeg = int(fdeg([j], e, p))
                # orbit computations use one j per Frobenius orbit: is j the minimal conjugate?
                isrep = (ckey([j], e, p) == pari('enc')(j, e)) if (needO or e in (2, 3, 4, 5, 6)) else True
            if need27:
                k = gcd(a_, 2) * gcd(b_, 2) - 1
                ex = exact(a_, b_, 27)
                if ex:
                    A27.setdefault((e, k), 0)
                    if jdeg == e and e not in (2, 3, 4, 5, 6):
                        A27[(e, k)] += ex // aut
                    else:
                        L = ptsorder(E, 27); assert len(L) == ex
                        good = 0
                        for P in L:
                            bc = tate(E, P)
                            if int(fdeg(bc, e, p)) != e: continue
                            good += 1
                            if isrep:
                                kk = tuple(int(z) for z in ckey(bc, e, p))
                                if kk not in seen27:
                                    orb = {tuple(int(z) for z in ckey(tate(E, pari.ellmul(E, P, d)), e, p)) for d in (1, 2, 4, 5, 7, 8, 10, 11, 13)}
                                    seen27 |= orb; orbit27.append((k, len(orb)))
                        assert good % aut == 0
                        A27[(e, k)] += good // aut
            for sp in needO:
                if not isrep: continue
                m, n, grp = sp
                if m == 2 and b_ % 2: continue
                if m == 2 and a_ % 2: continue
                Qs = ptsorder(E, n)
                Ps = [None] if m == 1 else [P for P in ptsorder(E, 2)]
                for Q in Qs:
                    for P in Ps:
                        if m == 1:
                            v = tate(E, Q); imgs = lambda: [tate(E, pari.ellmul(E, Q, d)) for d in range(1, n // 2 + 1) if gcd(d, n) == 1]
                        else:
                            if P == pari.ellmul(E, Q, n // 2): continue
                            v = coordsPQ(E, P, Q)
                            def imgs(P=P, Q=Q):
                                out = []
                                if grp == "aut2":
                                    for bb in (0, 1):
                                        for cc in (0, 1):
                                            for d in range(1, n // 2 + 1):
                                                if gcd(d, n) != 1: continue
                                                P2 = pari.elladd(E, P, pari.ellmul(E, Q, (n // 2) * bb)); Q2 = pari.elladd(E, pari.ellmul(E, P, cc), pari.ellmul(E, Q, d))
                                                out.append(coordsPQ(E, P2, Q2))
                                else:   # X1(2,42): (S_i, S_j + 2uQ), S1=P, S2=21Q, S3=P+21Q
                                    S = [P, pari.ellmul(E, Q, 21), pari.elladd(E, P, pari.ellmul(E, Q, 21))]
                                    for i in range(3):
                                        for jx in range(3):
                                            if i == jx: continue
                                            for u in (1, 2, 4, 5, 8, 10):
                                                P2 = S[i]; Q2 = pari.elladd(E, S[jx], pari.ellmul(E, Q, 2 * u))
                                                assert int(pari.ellorder(E, P2)) == 2 and int(pari.ellorder(E, Q2)) == 42 and P2 != pari.ellmul(E, Q2, 21)
                                                out.append(coordsPQ(E, P2, Q2))
                                return out
                        if int(fdeg(v, e, p)) != e: continue
                        kk = tuple(int(z) for z in ckey(v, e, p))
                        if kk in seenO[sp]: continue
                        orb = {tuple(int(z) for z in ckey(w, e, p)) for w in imgs()}
                        seenO[sp] |= orb; ORBS[(sp, e)].append(len(orb))
    if X27:
        ORB27[e] = orbit27
    print("FIELD p=%d e=%d done, mass errors %d, %.1f s" % (p, e, massbad, time.time() - t0), flush=True)
    for sp in orbspecs:
        if ORBS[(sp, e)]:
            from collections import Counter
            print("ORBP X1%s over F_%d degree %d, group %s: %d closed points, %d orbits, sizes %s" % (
                "(%d)" % sp[1] if sp[0] == 1 else "(%d,%d)" % (sp[0], sp[1]), p, e, sp[2], len(seenO[sp]), len(ORBS[(sp, e)]), dict(Counter(ORBS[(sp, e)]))), flush=True)
for sp in specs:
    P_ = pts[sp][1:]
    places = []
    for d in range(1, D + 1):
        s = sum(mobius(d // e) * P_[e - 1] for e in divisors(d)); assert s % d == 0; places.append(s // d)
    name = "X1(%d)" % sp[1] if sp[0] == 1 else "X1(%d,%d)" % sp
    print("POINTS %s p=%d: %s" % (name, p, P_))
    print("PLACES %s p=%d: %s" % (name, p, places))
if X27:
    from collections import Counter
    for e in range(1, D + 1):
        print("X27 e=%d closed points of X1(27) of exact degree e with k nonzero rational 2-torsion points: k=0:%s k=1:%s k=3:%s (k>=1: %s)" % (
            e, A27.get((e, 0), 0) / e, A27.get((e, 1), 0) / e, A27.get((e, 3), 0) / e, (A27.get((e, 1), 0) + A27.get((e, 3), 0)) / e))
        if ORB27.get(e):
            for k in (0, 1, 3):
                Lk = [s for kk, s in ORB27[e] if kk == k]
                if Lk: print("X27ORB e=%d k=%d: %d diamond orbits, sizes %s" % (e, k, len(Lk), dict(Counter(Lk))))
            print("X27ORB e=%d k>=1: %d diamond orbits" % (e, len([1 for kk, s in ORB27[e] if kk >= 1])))
print("PLACECOUNTS_PARI_DONE")
