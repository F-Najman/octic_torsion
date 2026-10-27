// cartier_51d_d8.m: same as cartier_d8.m, Case "51" (X1(51), Prop. 3.43), with a SIMPLE (one-step) model:
// F = F2(x1)(y(P3)) where now x1 = s (the Tate parameter of degree 5 on X1(17), instead of r of
// degree 7 used in cartier_51c_d8.m), so y(P3) has degree 40 over F2(x1).
// Independent check of the degree-8 Cartier / formal-immersion certificates at 2
// (Codex phi8 proof, Lemma 3.4):  U = (C-1) H^0(Y_F2, Omega^1); for every effective
// cuspidal divisor B = sum m_P P of degree 8 over F_2, dim U - dim(U meet Omega(-B)) = 8.
//
//   Case "45": Y = X1(45)          (Prop. 3.39)  E9 Tate normal form + point of order 5
//   Case "63": Y = X1(63)/<8>      (Prop. 3.28)  E9 Tate normal form + x(P7)
//   Case "51": Y = X1(51)          (Prop. 3.43)  X1(17) from x(8P)=x(9P), + point of order 3
//
// Kernel spaces are Riemann-Roch spaces DifferentialSpace(m*P) (no power series).
// Run from this directory:  magma -n -b cartier_51d_d8.m
SetColumns(0);
Case := "51d";
if not assigned MemGB then MemGB := 20; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
load "independent/formal/formal_common.m";
pf := "independent/formal/logs/cartier_" cat Case cat "_d8.progress";
Write(pf, "start" : Overwrite := true);
d := 8;
Rep := func< x, n | [x : i in [1..n]] >;
t0 := Cputime();
F2 := GF(2);
Say(pf, Sprintf("CASE %o  d = %o  Magma %o", Case, d, GetVersion()));

// Case "51b": X1(51) over F2 as a TWO-step tower (the three-step tower X1(17) -> x(P3) -> y(P3)
// of cartier_d8.m Case:=51 did not finish computing its genus in 4.7 hours):
// FX = X1(17) (genus 5, from x(8P) = x(9P)), then F = FX(y(P3)), where y(P3) is a root of
// m(Y) = Res_X( psi_3(X), Y^2 + (a1 X + a3) Y - (X^3 + a2 X^2 + a4 X + a6) ), of degree 8.
// If m is irreducible, FX(y(P3)) = FX(x(P3), y(P3)) is the function field of X1(51).
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
  gg := Evaluate(f, [S1, x1]);     // r -> generator, s -> x1
  if not IsIrreducible(gg) then continue; end if;
  Ftry := FunctionField(gg/LeadingCoefficient(gg));
  gen := Genus(Ftry);
  Say(pf, Sprintf("factor of bidegree (%o,%o): genus %o", Degree(f,r), Degree(f,s), gen));
  if gen eq 5 then FX := Ftry; found := true; break; end if;
end for;
assert found;
rr := FX.1; ss := FX!x1;
b := rr*ss*(rr-1); c := ss*(rr-1);
E0 := EllipticCurve([1-c, -b, -b, 0, 0]);
assert 17*E0![0,0] eq E0!0 and E0![0,0] ne E0!0;
psi := DivisionPolynomial(E0, 3); psi := psi/LeadingCoefficient(psi);
Say(pf, Sprintf("psi_3 degree %o irreducible %o", Degree(psi), IsIrreducible(psi)));
a1,a2,a3,a4,a6 := Explode(aInvariants(E0));
R2<X,Yv> := PolynomialRing(FX, 2);
res := Resultant(Evaluate(psi, X), Yv^2 + (a1*X + a3)*Yv - (X^3 + a2*X^2 + a4*X + a6), X);
SY<Y> := PolynomialRing(FX);
m := Evaluate(res, [0, Y]); m := m/LeadingCoefficient(m);
Say(pf, Sprintf("minimal polynomial of y(P3) over X1(17): degree %o, irreducible %o", Degree(m), IsIrreducible(m)));
assert Degree(m) eq 8 and IsIrreducible(m);
// norm of m from FX[Y] down to F2(x1)[Y]: resultant in S with the X1(17) equation gg
R3<Sv,Yw> := PolynomialRing(P1, 2);
Mbiv := &+[ &+[ R3 | Eltseq(Coefficient(m, k))[l+1] * Sv^l : l in [0..Degree(FX)-1] ] * Yw^k : k in [0..Degree(m)] ];
ggS := Evaluate(DefiningPolynomial(FX), Sv);
N := Resultant(ggS, Mbiv, Sv);
SZ<Z> := PolynomialRing(P1);
mabs := Evaluate(N, [0, Z]); mabs := mabs/LeadingCoefficient(mabs);
Say(pf, Sprintf("minimal polynomial of y(P3) over F2(x1): degree %o, max x1-degree of coefficients %o",
    Degree(mabs), Max([Max(Degree(Numerator(cf)), Degree(Denominator(cf))) : cf in Eltseq(mabs)])));
irr := IsIrreducible(mabs);
Say(pf, Sprintf("irreducible over F2(x1): %o  (%o s)", irr, Cputime(t0)));
assert Degree(mabs) eq 8*Degree(FX) and irr;
F<y3> := FunctionField(mabs);
// recover the X1(17) coordinate ss in F as the common root of gg(S) and M(S, y3)
SF<T> := PolynomialRing(F);
G := GCD(Evaluate(ggS, [T, F!0]), Evaluate(Mbiv, [T, F!y3]));
Say(pf, Sprintf("gcd of gg(S) and M(S,y3) over F has degree %o", Degree(G)));
assert Degree(G) eq 1;
s0 := -Coefficient(G, 0)/Coefficient(G, 1);
iota := hom< FX -> F | s0 >;
assert Evaluate(Polynomial([iota(cf) : cf in Eltseq(m)]), y3) eq 0;
jF := iota(FX!jInvariant(E0));
expect := <65, 1152, Rep(1,16) cat Rep(2,8) cat Rep(8,4), 62, 1080616>;

g := Genus(F);
Say(pf, Sprintf("GENUS %o (claimed %o)  constant field %o  (%o s)", g, expect[1], ConstantField(F), Cputime(t0)));
cusps := DedupPlaces(Poles(jF));
degs := [Degree(P) : P in cusps];
jdeg := &+[-Valuation(jF, P)*Degree(P) : P in cusps];
Say(pf, Sprintf("CUSPS %o closed cusps, sorted degrees %o (claimed %o), j-degree %o (expected %o)  (%o s)",
       #cusps, Sort(degs), expect[3], jdeg, expect[2], Cputime(t0)));
Say(pf, Sprintf("j-valuations at cusps: %o", [Valuation(jF, P) : P in cusps]));

V, h := DifferentialSpace(DivisorGroup(F)!0);
assert Dimension(V) eq g;
// sanity checks of the Cartier operator: kills exact differentials, C(f^2 w) = f C(w)
w1 := h(V.1); fx := F.1 + 1;
assert IsZero(Cartier(Differential(fx)));
assert Cartier(fx^2 * w1) eq fx * Cartier(w1);
U, MC := CartierU(V, h, g);
Say(pf, Sprintf("DIMU dim U = %o (claimed %o); rank C = %o  (%o s)", Dimension(U), expect[4], Rank(MC), Cputime(t0)));

small := [i : i in [1..#cusps] | degs[i] le d];
mmax := AssociativeArray();
for i in small do mmax[i] := d div degs[i]; end for;
Say(pf, Sprintf("cusps of degree <= %o: %o; kernel spaces to compute: %o", d, #small, &+[mmax[i] : i in small]));
W := AllKernels(cusps, small, mmax, h, F2, g, t0, pf);
WU := AssociativeArray();
for k in Keys(W) do WU[k] := U meet W[k]; end for;
Say(pf, Sprintf("kernel spaces done: %o  (%o s)", #Keys(W), Cputime(t0)));

expected_count := CountDivisors([degs[i] : i in small], d);
cnt := 0; bad := 0; minr := 10^6; maxr := -1;
EnumRanks(1, d, U, small, degs, WU, Dimension(U), d, [], ~cnt, ~bad, ~minr, ~maxr, pf);
Say(pf, Sprintf("RESULT case %o: genus %o, dim U %o, cusp degrees %o, divisors enumerated %o (generating function %o, claimed %o), rank min %o max %o, failures %o",
       Case, g, Dimension(U), Sort(degs), cnt, expected_count, expect[5], minr, maxr, bad));
// generic fibre check (last): x(P3) in F, and P17 + P3 of exact order 51
SX<Xr> := PolynomialRing(F);
xs := [rt[1] : rt in Roots(Polynomial([iota(cf) : cf in Eltseq(psi)]))];
EF := EllipticCurve([iota(a) : a in aInvariants(E0)]);
okpt := false;
for x0 in xs do
  if IsOnCurve(EF, [x0, y3, 1]) then
    P3 := EF![x0, y3]; Q := EF![0,0] + P3;
    okpt := 3*P3 eq EF!0 and P3 ne EF!0 and 51*Q eq EF!0 and 17*Q ne EF!0 and 3*Q ne EF!0;
  end if;
end for;
Say(pf, Sprintf("generic fibre: (x(P3), y3) is a point of order 3 and (0,0)+P3 has exact order 51: %o  (%o s)", okpt, Cputime(t0)));
Say(pf, Sprintf("CARTIER_D8_DONE case %o (%o s)", Case, Cputime(t0)));
quit;
