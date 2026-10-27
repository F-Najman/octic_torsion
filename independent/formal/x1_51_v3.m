// x1_51_v2.m: X1(51) over F2 (Prop. 3.43), two-step tower model of cartier_51b_d8.m (which computed
// x1_51_v3.m: independent check of the X1(51) Cartier certificate (Prop. 3.43), d = 8, with a
// third model: G1 = X1(51)/<16> written as a SIMPLE extension of F2(x1) (x1 = Tate parameter r of
// X1(17); generator x3 = x(P3), minimal polynomial of degree 28 = norm of psi_3 from X1(17)),
// then F = G1(y(P3)) by the Weierstrass equation (degree 2).  Kernel spaces by Riemann-Roch
// (DifferentialSpace(m*P)), written to logs/x1_51_v3_kernels.txt as they are computed.
// Run from this directory:  magma -n -b MemGB:=48 x1_51_v3.m
SetColumns(0);
if not assigned MemGB then MemGB := 48; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
load "independent/formal/formal_common.m";
pf := "independent/formal/logs/x1_51_v3.progress";
kf := "independent/formal/data/x1_51_v3_kernels.txt";
Write(pf, "start" : Overwrite := true);
Write(kf, "kernel spaces of X1(51)/F2 (model v3), coordinates w.r.t. DifferentialSpace(0) basis" : Overwrite := true);
d := 8;
Rep := func< x, n | [x : i in [1..n]] >;
t0 := Cputime();
F2 := GF(2);
A2<r,s> := PolynomialRing(F2, 2);
Q2 := FieldOfFractions(A2);
bb := Q2!(r*s*(r-1)); cc := Q2!(s*(r-1));
EE := EllipticCurve([1-cc, -bb, -bb, 0, 0]);
PP := EE![0,0];
P8 := 8*PP; P9 := 9*PP;
rel := Numerator(P8[1]/P8[3] - P9[1]/P9[3]);
fac := Factorization(rel);
Say(pf, Sprintf("factors of x(8P)-x(9P): %o", [<TotalDegree(f[1]), f[2]> : f in fac]));
found := false;
for fpair in fac do
  f := fpair[1];
  if Degree(f, s) eq 0 or Degree(f, r) eq 0 then continue; end if;
  P1<x1> := FunctionField(F2);
  RS<S1> := PolynomialRing(P1);
  gg := Evaluate(f, [x1, S1]);
  if not IsIrreducible(gg) then continue; end if;
  Ftry := FunctionField(gg/LeadingCoefficient(gg));
  gen := Genus(Ftry);
  Say(pf, Sprintf("factor of bidegree (%o,%o): genus %o", Degree(f,r), Degree(f,s), gen));
  if gen eq 5 then FX := Ftry; found := true; break; end if;
end for;
assert found;
rr := FX!x1; ss := FX.1;
b := rr*ss*(rr-1); c := ss*(rr-1);
E0 := EllipticCurve([1-c, -b, -b, 0, 0]);
assert 17*E0![0,0] eq E0!0 and E0![0,0] ne E0!0;
psi := DivisionPolynomial(E0, 3); psi := psi/LeadingCoefficient(psi);
Say(pf, Sprintf("psi_3 degree %o irreducible %o", Degree(psi), IsIrreducible(psi)));
// X1(51)/<16> = F2(x1)(x3): minimal polynomial of x3 over F2(x1) = norm of psi_3
R3<Sv,Xw> := PolynomialRing(P1, 2);
Mbiv := &+[ &+[ R3 | Eltseq(Coefficient(psi, k))[l+1] * Sv^l : l in [0..Degree(FX)-1] ] * Xw^k : k in [0..Degree(psi)] ];
ggS := Evaluate(DefiningPolynomial(FX), Sv);
N := Resultant(ggS, Mbiv, Sv);
SZ<Z> := PolynomialRing(P1);
Nx := Evaluate(N, [0, Z]); Nx := Nx/LeadingCoefficient(Nx);
Say(pf, Sprintf("minimal polynomial of x(P3) over F2(x1): degree %o, max x1-degree %o, irreducible %o",
    Degree(Nx), Max([Max(Degree(Numerator(cf)), Degree(Denominator(cf))) : cf in Eltseq(Nx)]), IsIrreducible(Nx)));
assert Degree(Nx) eq 28 and IsIrreducible(Nx);
G1<x3> := FunctionField(Nx);
SG<T> := PolynomialRing(G1);
G := GCD(Evaluate(ggS, [T, G1!0]), Evaluate(Mbiv, [T, G1!x3]));
assert Degree(G) eq 1;
s0 := -Coefficient(G, 0)/Coefficient(G, 1);
iota := hom< FX -> G1 | s0 >;
assert Evaluate(Polynomial([iota(cf) : cf in Eltseq(psi)]), x3) eq 0;
Say(pf, Sprintf("genus of X1(51)/<16> model: %o (expected 33)  (%o s)", Genus(G1), Cputime(t0)));
A := [iota(a) : a in aInvariants(E0)];
SY<Y> := PolynomialRing(G1);
q := Y^2 + (A[1]*x3 + A[3])*Y - (x3^3 + A[2]*x3^2 + A[4]*x3 + A[5]);
assert IsIrreducible(q);
F<y3> := FunctionField(q);
g := Genus(F);
Say(pf, Sprintf("GENUS %o (claimed 65)  constant field %o  (%o s)", g, ConstantField(F), Cputime(t0)));
EF := EllipticCurve([F!a : a in A]);
P3 := EF![F!x3, y3]; Q := EF![0,0] + P3;
assert 3*P3 eq EF!0 and P3 ne EF!0 and 51*Q eq EF!0 and 17*Q ne EF!0 and 3*Q ne EF!0;
Say(pf, "point of exact order 51 on the generic fibre: true");

jG := G1!iota(FX!jInvariant(E0));
cuspsG := Poles(jG);
Say(pf, Sprintf("cusps of X1(51)/<16>: degrees %o  (%o s)", Sort([Degree(Q0) : Q0 in cuspsG]), Cputime(t0)));
cusps := []; jv := [];
for Q0 in cuspsG do
  for P in DedupPlaces(Decomposition(F, Q0)) do
    Append(~cusps, P); Append(~jv, RamificationIndex(P) * Valuation(jG, Q0));
  end for;
end for;
assert #DedupPlaces(cusps) eq #cusps;
degs := [Degree(P) : P in cusps];
Say(pf, Sprintf("CUSPS %o closed cusps, sorted degrees %o (claimed 1^16 2^8 8^4), j-degree %o (expected 1152)  (%o s)",
       #cusps, Sort(degs), &+[-jv[i]*degs[i] : i in [1..#cusps]], Cputime(t0)));

V, h := DifferentialSpace(DivisorGroup(F)!0);
assert Dimension(V) eq g;
w1 := h(V.1); fx := F.1 + 1;
assert IsZero(Cartier(Differential(fx)));
assert Cartier(fx^2 * w1) eq fx * Cartier(w1);
U, MC := CartierU(V, h, g);
Say(pf, Sprintf("DIMU dim U = %o (claimed 62); rank C = %o  (%o s)", Dimension(U), Rank(MC), Cputime(t0)));

small := [i : i in [1..#cusps] | degs[i] le d];
W := AssociativeArray();
for i in small do
  for m in [1..(d div degs[i])] do
    t1 := Cputime();
    W[<i,m>] := KernelCoords(cusps[i], m, h, F2, g);
    assert Dimension(W[<i,m>]) ge g - m*degs[i];
    if m gt 1 then assert W[<i,m>] subset W[<i,m-1>]; end if;
    Write(kf, Sprintf("SPACE %o %o %o", i, m, Dimension(W[<i,m>])));
    for v in Basis(W[<i,m>]) do Write(kf, &cat[Sprint(x) : x in Eltseq(v)]); end for;
    Say(pf, Sprintf("KERNEL cusp %o deg %o m %o codim %o [%o s, total %o s]", i, degs[i], m, g - Dimension(W[<i,m>]), Cputime(t1), Cputime(t0)));
  end for;
end for;
WU := AssociativeArray();
for k in Keys(W) do WU[k] := U meet W[k]; end for;
Say(pf, Sprintf("kernel spaces done: %o  (%o s)", #Keys(W), Cputime(t0)));
expected_count := CountDivisors([degs[i] : i in small], d);
cnt := 0; bad := 0; minr := 10^6; maxr := -1;
EnumRanks(1, d, U, small, degs, WU, Dimension(U), d, [], ~cnt, ~bad, ~minr, ~maxr, pf);
Say(pf, Sprintf("RESULT X1(51)/F2: genus %o, dim U %o, cusp degrees %o, divisors enumerated %o (generating function %o, claimed 1080616), rank min %o max %o, failures %o",
       g, Dimension(U), Sort(degs), cnt, expected_count, minr, maxr, bad));
Say(pf, Sprintf("X1_51_V3_DONE (%o s)", Cputime(t0)));
quit;
