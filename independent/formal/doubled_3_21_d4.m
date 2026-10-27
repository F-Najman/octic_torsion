// Independent check of the C3 x C21 doubled-cuspidal certificate in degree 8
// (Codex phi8 proof, Prop. 3.66): on the fixed-pairing component X of X1(3,21) over F4,
// l(2C) = 1 for every effective cuspidal divisor C of degree 4, i.e. the space
// Omega(-2C) = {omega : div(omega) >= 2C} has codimension 8 in H^0(Omega)
// (Riemann-Roch: l(2C) = 8 + 1 - g + (g - codim) = 9 - codim).
// Model (own construction): Tate normal form over X1(7)_F4, adjoin a point P3 of order 3
// (x and y), then the x-coordinate of a second 3-torsion point (root of psi_3/(X-x(P3))).
// Run from this directory:  magma -n -b doubled_3_21_d4.m
SetColumns(0);
SetMemoryLimit(16*10^9);
load "independent/formal/formal_common.m";
pf := "independent/formal/logs/doubled_3_21_d4.progress";
Write(pf, "start" : Overwrite := true);
d := 4;
t0 := Cputime();
F4<zeta> := GF(4);
P1<t> := FunctionField(F4);
b := t^3-t^2; c := t^2-t;
E7 := EllipticCurve([1-c, -b, -b, 0, 0]);
assert 7*E7![0,0] eq E7!0 and E7![0,0] ne E7!0;
psi := DivisionPolynomial(E7, 3); psi := psi/LeadingCoefficient(psi);
Say(pf, Sprintf("psi_3 degree %o irreducible over F4(t) %o", Degree(psi), IsIrreducible(psi)));
F1<x3> := FunctionField(psi);
a1,a2,a3,a4,a6 := Explode([F1!a : a in aInvariants(E7)]);
S1<Y> := PolynomialRing(F1);
qy := Y^2 + a1*x3*Y + a3*Y - (x3^3 + a2*x3^2 + a4*x3 + a6);
Say(pf, Sprintf("y-equation irreducible %o", IsIrreducible(qy)));
F2f<y3> := FunctionField(qy);
S2<X> := PolynomialRing(F2f);
psiF := Polynomial([F2f!F1!cf : cf in Eltseq(psi)]);
cub := psiF div (X - F2f!x3);
assert psiF eq cub*(X - F2f!x3);
fc := Factorization(cub);
Say(pf, Sprintf("degrees of the factors of psi_3/(X - x(P3)) over X1(21): %o", [Degree(f[1]) : f in fc]));
assert #fc eq 1 and Degree(fc[1][1]) eq 3;
F<xq> := FunctionField(fc[1][1]);
g := Genus(F);
Say(pf, Sprintf("GENUS %o (claimed 25)  constant field %o  (%o s)", g, ConstantField(F), Cputime(t0)));
// generic fibre check: P7 = (0,0), P3, and a second 3-torsion point Q3 with x(Q3) = xq
EF := BaseChange(E7, F);
P3 := EF![F!F2f!x3, F!y3];
assert 3*P3 eq EF!0 and P3 ne EF!0;
assert Evaluate(Polynomial([F!F2f!F1!cf : cf in Eltseq(psi)]), xq) eq 0 and xq ne F!F2f!x3;
Say(pf, "generic fibre: (0,0) of order 7, P3 of order 3, x(Q3) another 3-torsion abscissa: true");

jF := F!jInvariant(E7);
cusps := DedupPlaces(Poles(jF));
degs := [Degree(P) : P in cusps];
Say(pf, Sprintf("CUSPS %o closed cusps over F4, sorted degrees %o (claimed 1^24 3^8), total %o, j-degree %o  (%o s)",
       #cusps, Sort(degs), &+degs, &+[-Valuation(jF,P)*Degree(P) : P in cusps], Cputime(t0)));

V, h := DifferentialSpace(DivisorGroup(F)!0);
assert Dimension(V) eq g;
small := [i : i in [1..#cusps] | degs[i] le d];
mlist := AssociativeArray();
for i in small do mlist[i] := [2*m : m in [1..(d div degs[i])]]; end for;   // 2m*P for m*deg P <= d
Say(pf, Sprintf("cusps of degree <= %o: %o", d, #small));
W := AllKernelsList(cusps, small, mlist, h, F4, g, t0, pf);
// doubled spaces: W2[<i,m>] = Omega(-2m P_i)
W2 := AssociativeArray();
for i in small do
  for m in [1..(d div degs[i])] do W2[<i,m>] := W[<i,2*m>]; end for;
end for;
Say(pf, Sprintf("doubled kernel spaces: %o  (%o s)", #Keys(W2), Cputime(t0)));

expected_count := CountDivisors([degs[i] : i in small], d);
cnt := 0; bad := 0; minr := 10^6; maxr := -1;
EnumRanks(1, d, VectorSpace(F4, g), small, degs, W2, g, 2*d, [], ~cnt, ~bad, ~minr, ~maxr, pf);
Say(pf, Sprintf("RESULT X1(3,21)/F4: genus %o, cusp degrees %o, degree-%o divisors %o (generating function %o, claimed 17742), codim of Omega(-2C) min %o max %o (need %o), failures %o",
       g, Sort(degs), d, cnt, expected_count, minr, maxr, 2*d, bad));
Say(pf, Sprintf("DOUBLED_3_21_D4_DONE (%o s)", Cputime(t0)));
quit;
