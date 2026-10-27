// Independent check for Prop. 3.76 of the Codex phi8 proof: on X = X1(25) over F16
// (residue field of the prime above 2 in Q(zeta5)), for all 171 effective degree-2
// divisors C supported on the 18 rational cusps, h^0(K - C) = 10, i.e. h^0(C) = 1.
// Model (own construction): Tate normal form E5 (b = c = t) over F16(t), adjoin
// x(Q) for Q with 5Q = +-(0,0), i.e. a root of the numerator of x([5]Q).
// Run from this directory:  magma -n -b x1_25_F16.m
SetColumns(0);
SetMemoryLimit(12*10^9);
load "independent/formal/formal_common.m";
pf := "independent/formal/logs/x1_25_F16.progress";
Write(pf, "start" : Overwrite := true);
t0 := Cputime();
F16<z> := GF(16);
P1<t> := FunctionField(F16);
E5 := EllipticCurve([1-t, -t, -t, 0, 0]);
assert 5*E5![0,0] eq E5!0 and E5![0,0] ne E5!0;
phi := XMultNumerator(E5, 5);
phi := phi/LeadingCoefficient(phi);
Say(pf, Sprintf("numerator of x([5]Q): degree %o, irreducible over F16(t) %o", Degree(phi), IsIrreducible(phi)));
F<xq> := FunctionField(phi);
g := Genus(F);
Say(pf, Sprintf("GENUS %o (expected 12)  constant field %o  (%o s)", g, ConstantField(F), Cputime(t0)));

// generic point of order 25 with 5Q = +-(0,0)
a1,a2,a3,a4,a6 := Explode([F!a : a in aInvariants(E5)]);
SY<Y> := PolynomialRing(F);
rts := Roots(Y^2 + (a1*xq + a3)*Y - (xq^3 + a2*xq^2 + a4*xq + a6));
EF := BaseChange(E5, F);
ok := #rts eq 2;
for rt in rts do
  Q := EF![xq, rt[1]];
  ok := ok and (5*Q eq EF![0,0] or 5*Q eq -EF![0,0]) and 25*Q eq EF!0 and 5*Q ne EF!0;
end for;
Say(pf, Sprintf("generic fibre: both y-roots in F give points Q of exact order 25 with 5Q = +-(0,0): %o", ok));

jF := F!jInvariant(E5);
cusps := DedupPlaces(Poles(jF));
degs := [Degree(P) : P in cusps];
Say(pf, Sprintf("CUSPS %o closed cusps over F16, sorted degrees %o (claimed 1^18 5^2), j-degree %o (expected 300)",
       #cusps, Sort(degs), &+[-Valuation(jF,P)*Degree(P) : P in cusps]));
pl1 := Places(F, 1);
Say(pf, Sprintf("number of F16-rational places of X1(25): %o (claimed 18), all cusps: %o", #pl1, &and[P in cusps : P in pl1]));

V, h := DifferentialSpace(DivisorGroup(F)!0);
assert Dimension(V) eq g;
rat := [i : i in [1..#cusps] | degs[i] eq 1];
mmax := AssociativeArray();
for i in rat do mmax[i] := 2; end for;
W := AllKernels(cusps, rat, mmax, h, F16, g, t0, pf);
cnt := 0; bad := 0; minr := 10^6; maxr := -1;
EnumRanks(1, 2, VectorSpace(F16, g), rat, degs, W, g, 2, [], ~cnt, ~bad, ~minr, ~maxr, pf);
Say(pf, Sprintf("via kernel spaces: divisors %o (expected 171), codim of Omega(-C) min %o max %o (need 2, i.e. h^0(K-C) = %o), failures %o",
       cnt, minr, maxr, g-2, bad));

// direct Riemann-Roch on every divisor, as a cross-check
bad2 := 0; n2 := 0;
for a in [1..#rat] do
  for b in [a..#rat] do
    C := cusps[rat[a]] + cusps[rat[b]];
    n2 +:= 1;
    hKC := Dimension(DifferentialSpace(C));
    hC := Dimension(RiemannRochSpace(C));
    if hKC ne g - 2 or hC ne 1 then
      bad2 +:= 1; Say(pf, Sprintf("FAILURE direct %o %o: h0(K-C) %o h0(C) %o", rat[a], rat[b], hKC, hC));
    end if;
  end for;
end for;
Say(pf, Sprintf("RESULT X1(25)/F16: genus %o, cusp degrees %o, rational places %o; divisors %o, h0(K-C) = %o and h0(C) = 1 for all (direct RR failures %o, kernel-space failures %o)",
       g, Sort(degs), #pl1, n2, g-2, bad2, bad));
Say(pf, Sprintf("X1_25_F16_DONE (%o s)", Cputime(t0)));
quit;
