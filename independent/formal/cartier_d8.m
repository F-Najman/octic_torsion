// Independent check of the degree-8 Cartier / formal-immersion certificates at 2
// (Codex phi8 proof, Lemma 3.4):  U = (C-1) H^0(Y_F2, Omega^1); for every effective
// cuspidal divisor B = sum m_P P of degree 8 over F_2, dim U - dim(U meet Omega(-B)) = 8.
//
//   Case "45": Y = X1(45)          (Prop. 3.39)  E9 Tate normal form + point of order 5
//   Case "63": Y = X1(63)/<8>      (Prop. 3.28)  E9 Tate normal form + x(P7)
//   Case "51": Y = X1(51)          (Prop. 3.43)  X1(17) from x(8P)=x(9P), + point of order 3
//
// Kernel spaces are Riemann-Roch spaces DifferentialSpace(m*P) (no power series).
// Run from this directory:  magma -n -b Case:=45 cartier_d8.m
SetColumns(0);
if not assigned Case then Case := "45"; end if;
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

if Case in {"45", "63"} then
  P1<t> := FunctionField(F2);
  b := t^2*(t-1)*(t^2-t+1); c := t^2*(t-1);
  E0 := EllipticCurve([1-c, -b, -b, 0, 0]);
  assert 9*E0![0,0] eq E0!0 and 3*E0![0,0] ne E0!0;
  if Case eq "63" then
    psi := DivisionPolynomial(E0, 7); psi := psi/LeadingCoefficient(psi);
    Say(pf, Sprintf("psi_7 degree %o irreducible %o", Degree(psi), IsIrreducible(psi)));
    F<u> := FunctionField(psi);
    expect := <49, 864, Rep(1,9) cat Rep(2,3) cat Rep(3,9) cat [6], 45, 44046>;
  else
    psi := DivisionPolynomial(E0, 5); psi := psi/LeadingCoefficient(psi);
    Say(pf, Sprintf("psi_5 degree %o irreducible %o", Degree(psi), IsIrreducible(psi)));
    F1<x5> := FunctionField(psi);
    a1,a2,a3,a4,a6 := Explode([F1!a : a in aInvariants(E0)]);
    S<Y> := PolynomialRing(F1);
    q := Y^2 + a1*x5*Y + a3*Y - (x5^3 + a2*x5^2 + a4*x5 + a6);
    Say(pf, Sprintf("y-equation irreducible %o", IsIrreducible(q)));
    F<y5> := FunctionField(q);
    EF := BaseChange(E0, F);
    P5 := EF![F!x5, y5]; P9 := EF![0,0];
    assert P5 ne EF!0 and 5*P5 eq EF!0;
    Q := P9 + P5;
    assert 45*Q eq EF!0 and 15*Q ne EF!0 and 9*Q ne EF!0;
    Say(pf, Sprintf("point of exact order 45 on the generic fibre: true"));
    expect := <41, 864, Rep(1,12) cat Rep(2,4) cat Rep(4,5) cat Rep(6,2) cat [12], 39, 148945>;
  end if;
  jF := F!jInvariant(E0);
else
  // X1(17) over F2 from the Tate normal form b = rs(r-1), c = s(r-1), relation x(8P) = x(9P)
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
  F1<x3> := FunctionField(psi);
  a1,a2,a3,a4,a6 := Explode([F1!a : a in aInvariants(E0)]);
  S<Y> := PolynomialRing(F1);
  q := Y^2 + a1*x3*Y + a3*Y - (x3^3 + a2*x3^2 + a4*x3 + a6);
  Say(pf, Sprintf("y-equation irreducible %o", IsIrreducible(q)));
  F<y3> := FunctionField(q);
  EF := BaseChange(E0, F);
  P3 := EF![F!x3, y3]; P17 := EF![0,0];
  assert P3 ne EF!0 and 3*P3 eq EF!0;
  Q := P17 + P3;
  assert 51*Q eq EF!0 and 17*Q ne EF!0 and 3*Q ne EF!0;
  Say(pf, Sprintf("point of exact order 51 on the generic fibre: true"));
  jF := F!jInvariant(E0);
  expect := <65, 1152, Rep(1,16) cat Rep(2,8) cat Rep(8,4), 62, 1080616>;
end if;

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
Say(pf, Sprintf("CARTIER_D8_DONE case %o (%o s)", Case, Cputime(t0)));
quit;
