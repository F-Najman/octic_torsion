// x1_51_v2.m: X1(51) over F2 (Prop. 3.43), two-step tower model of cartier_51b_d8.m (which computed
// genus 65 in 5768 s and then ran out of its 20 GB cap inside Poles(j)).  Changes: cusps are found
// by decomposing the cusps of X1(17); kernel spaces from local expansions (KernelsByExpansion,
// validated against DifferentialSpace on X1(45) and X1(63)/<8>), saved to logs/x1_51_v2_kernels.txt.
// Independent check of the degree-8 Cartier / formal-immersion certificates at 2
// (Codex phi8 proof, Lemma 3.4):  U = (C-1) H^0(Y_F2, Omega^1); for every effective
// cuspidal divisor B = sum m_P P of degree 8 over F_2, dim U - dim(U meet Omega(-B)) = 8.
//
//   Case "45": Y = X1(45)          (Prop. 3.39)  E9 Tate normal form + point of order 5
//   Case "63": Y = X1(63)/<8>      (Prop. 3.28)  E9 Tate normal form + x(P7)
//   Case "51": Y = X1(51)          (Prop. 3.43)  X1(17) from x(8P)=x(9P), + point of order 3
//
// Kernel spaces are Riemann-Roch spaces DifferentialSpace(m*P) (no power series).
// Run from this directory:  magma -n -b MemGB:=48 x1_51_v2.m
SetColumns(0);
Case := "51v2";
if not assigned MemGB then MemGB := 20; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
load "independent/formal/formal_common.m";
pf := "independent/formal/logs/x1_51_v2.progress";
kf := "independent/formal/data/x1_51_v2_kernels.txt";
Write(kf, "kernel spaces of X1(51)/F2, coordinates w.r.t. DifferentialSpace(0) basis" : Overwrite := true);
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
a1,a2,a3,a4,a6 := Explode(aInvariants(E0));
R2<X,Yv> := PolynomialRing(FX, 2);
res := Resultant(Evaluate(psi, X), Yv^2 + (a1*X + a3)*Yv - (X^3 + a2*X^2 + a4*X + a6), X);
SY<Y> := PolynomialRing(FX);
m := Evaluate(res, [0, Y]); m := m/LeadingCoefficient(m);
Say(pf, Sprintf("minimal polynomial of y(P3) over X1(17): degree %o, irreducible %o", Degree(m), IsIrreducible(m)));
assert Degree(m) eq 8 and IsIrreducible(m);
F<y3> := FunctionField(m);
jF := F!jInvariant(E0);
expect := <65, 1152, Rep(1,16) cat Rep(2,8) cat Rep(8,4), 62, 1080616>;

g := Genus(F);
Say(pf, Sprintf("GENUS %o (claimed %o)  constant field %o  (%o s)", g, expect[1], ConstantField(F), Cputime(t0)));
jX := FX!jInvariant(E0);
cuspsX := Poles(jX);
Say(pf, Sprintf("cusps of X1(17) over F2: degrees %o, j-degree %o", [Degree(Q) : Q in cuspsX], &+[-Valuation(jX,Q)*Degree(Q) : Q in cuspsX]));
cusps := []; jv := [];
for Q0 in cuspsX do
  for P in DedupPlaces(Decomposition(F, Q0)) do
    Append(~cusps, P); Append(~jv, RamificationIndex(P) * Valuation(jX, Q0));
  end for;
end for;
assert #DedupPlaces(cusps) eq #cusps;
degs := [Degree(P) : P in cusps];
jdeg := &+[-jv[i]*degs[i] : i in [1..#cusps]];
Say(pf, Sprintf("CUSPS %o closed cusps, sorted degrees %o (claimed %o), j-degree %o (expected %o)  (%o s)",
       #cusps, Sort(degs), expect[3], jdeg, expect[2], Cputime(t0)));
Say(pf, Sprintf("j-valuations at cusps (e(P/Q) v_Q(j)): %o", jv));

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
W := AllKernelsByExpansion(cusps, small, mmax, h, F2, g, t0, pf, kf);
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
xs := [rt[1] : rt in Roots(Polynomial([F!cf : cf in Eltseq(psi)]))];
EF := BaseChange(E0, F);
okpt := false;
for x0 in xs do
  if IsOnCurve(EF, [x0, y3, 1]) then
    P3 := EF![x0, y3]; Q := EF![0,0] + P3;
    okpt := 3*P3 eq EF!0 and P3 ne EF!0 and 51*Q eq EF!0 and 17*Q ne EF!0 and 3*Q ne EF!0;
  end if;
end for;
Say(pf, Sprintf("generic fibre: (x(P3), y3) is a point of order 3 and (0,0)+P3 has exact order 51: %o  (%o s)", okpt, Cputime(t0)));
// direct Riemann-Roch cross-check of two kernel spaces on this curve
r1 := [i : i in small | degs[i] eq 1][1]; q1 := [i : i in small | degs[i] eq 2][1];
for pr in [<r1, 1>, <q1, 1>] do
  t1 := Cputime();
  WD := KernelCoords(cusps[pr[1]], pr[2], h, F2, g);
  Say(pf, Sprintf("CROSSCHECK DifferentialSpace(%o*P_%o) (deg %o) equals expansion kernel: %o  [%o s]",
      pr[2], pr[1], degs[pr[1]], WD eq W[pr], Cputime(t1)));
end for;
Say(pf, Sprintf("X1_51_V2_DONE (%o s)", Cputime(t0)));
quit;
