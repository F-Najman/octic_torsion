// Independent check for Prop. 3.41 of the Codex phi8 proof: X = X1(49) over F2,
// U = (C-1) H^0(X_F2, Omega^1) has dimension 68, and restriction from U has rank 8 on
// every effective cuspidal degree-8 divisor supported on rational cusps and the six cubic
// cusps that contains a cubic cusp.  The proof uses a chosen set A of 17 of the 21
// rational cusps (125307 divisors); here ALL 21 rational cusps are allowed (323631
// divisors), which covers every choice of A, and the count for the first 17 is reported.
// Model (own construction): Tate normal form E7 over F2(t), adjoin x(Q) for Q with
// 7Q = +-(0,0), i.e. a root of the numerator of x([7]Q) (degree 49).
// Run from this directory:  magma -n -b x1_49_cartier.m
SetColumns(0);
if not assigned MemGB then MemGB := 24; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
load "independent/formal/formal_common.m";
pf := "independent/formal/logs/x1_49_cartier.progress";
Write(pf, "start" : Overwrite := true);
d := 8;
t0 := Cputime();
F2 := GF(2);
P1<t> := FunctionField(F2);
b := t^3-t^2; c := t^2-t;
E7 := EllipticCurve([1-c, -b, -b, 0, 0]);
assert 7*E7![0,0] eq E7!0 and E7![0,0] ne E7!0;
phi := XMultNumerator(E7, 7);
phi := phi/LeadingCoefficient(phi);
Say(pf, Sprintf("numerator of x([7]Q): degree %o, max t-degree of coefficients %o", Degree(phi),
    Max([Max(Degree(Numerator(cf)), Degree(Denominator(cf))) : cf in Eltseq(phi)])));
Say(pf, Sprintf("irreducible over F2(t): %o  (%o s)", IsIrreducible(phi), Cputime(t0)));
F<xq> := FunctionField(phi);
g := Genus(F);
Say(pf, Sprintf("GENUS %o (expected 69)  constant field %o  (%o s)", g, ConstantField(F), Cputime(t0)));
jF := F!jInvariant(E7);
cusps := DedupPlaces(Poles(jF));
degs := [Degree(P) : P in cusps];
Say(pf, Sprintf("CUSPS %o closed cusps, sorted degrees %o (claimed 1^21 3^6 21), j-degree %o (expected 1176)  (%o s)",
       #cusps, Sort(degs), &+[-Valuation(jF,P)*Degree(P) : P in cusps], Cputime(t0)));

V, h := DifferentialSpace(DivisorGroup(F)!0);
assert Dimension(V) eq g;
w1 := h(V.1); fx := F.1 + 1;
assert IsZero(Cartier(Differential(fx)));
assert Cartier(fx^2 * w1) eq fx * Cartier(w1);
U, MC := CartierU(V, h, g);
Say(pf, Sprintf("DIMU dim U = %o (claimed 68); rank C = %o  (%o s)", Dimension(U), Rank(MC), Cputime(t0)));

rat := [i : i in [1..#cusps] | degs[i] eq 1];
cub := [i : i in [1..#cusps] | degs[i] eq 3];
assert #rat eq 21 and #cub eq 6;
mmax := AssociativeArray();
for i in cub do mmax[i] := 2; end for;
for i in rat do mmax[i] := 5; end for;    // a divisor containing a cubic cusp has rational part of degree <= 5
W := AllKernels(cusps, cub cat rat, mmax, h, F2, g, t0, pf);
WU := AssociativeArray();
for k in Keys(W) do WU[k] := U meet W[k]; end for;
Say(pf, Sprintf("kernel spaces done: %o  (%o s)", #Keys(W), Cputime(t0)));

// cubic parts of degree 3 or 6
cubparts := [ [<a,1>] : a in cub ] cat [ [<a,2>] : a in cub ] cat [ [<cub[x],1>, <cub[y],1>] : x, y in [1..6] | x lt y ];
assert #cubparts eq 27;
procedure RunFamily(ratidx, cubparts, cusps, degs, U, WU, d, pf, ~cnt, ~bad, ~minr, ~maxr)
  for cp in cubparts do
    Kc := U; dc := 0;
    for pr in cp do Kc meet:= WU[pr]; dc +:= pr[2]*degs[pr[1]]; end for;
    EnumRanks(1, d - dc, Kc, ratidx, degs, WU, Dimension(U), d, cp, ~cnt, ~bad, ~minr, ~maxr, pf);
  end for;
end procedure;
cnt := 0; bad := 0; minr := 10^6; maxr := -1;
RunFamily(rat, cubparts, cusps, degs, U, WU, d, pf, ~cnt, ~bad, ~minr, ~maxr);
Say(pf, Sprintf("all 21 rational cusps: divisors %o (expected 6*C(25,5)+21*C(22,2) = %o), rank min %o max %o, failures %o",
       cnt, 6*Binomial(25,5) + 21*Binomial(22,2), minr, maxr, bad));
cnt17 := 0; bad17 := 0; minr17 := 10^6; maxr17 := -1;
RunFamily(rat[1..17], cubparts, cusps, degs, U, WU, d, pf, ~cnt17, ~bad17, ~minr17, ~maxr17);
Say(pf, Sprintf("first 17 rational cusps: divisors %o (claimed 125307), rank min %o max %o, failures %o",
       cnt17, minr17, maxr17, bad17));
Say(pf, Sprintf("RESULT X1(49)/F2: genus %o, cusp degrees %o, dim U %o, divisors with a cubic cusp on 21 rational + 6 cubic cusps %o, failures %o",
       g, Sort(degs), Dimension(U), cnt, bad));

// generic fibre check (last, may be slow): a y-coordinate in F and the order of Q
a1,a2,a3,a4,a6 := Explode([F!a : a in aInvariants(E7)]);
SY<Y> := PolynomialRing(F);
rts := Roots(Y^2 + (a1*xq + a3)*Y - (xq^3 + a2*xq^2 + a4*xq + a6));
EF := BaseChange(E7, F);
ok := #rts eq 2;
for rt in rts do
  Q := EF![xq, rt[1]];
  ok := ok and (7*Q eq EF![0,0] or 7*Q eq -EF![0,0]) and 49*Q eq EF!0 and 7*Q ne EF!0;
end for;
Say(pf, Sprintf("generic fibre: both y-roots give points of exact order 49 with 7Q = +-(0,0): %o  (%o s)", ok, Cputime(t0)));
Say(pf, Sprintf("X1_49_CARTIER_DONE (%o s)", Cputime(t0)));
quit;
