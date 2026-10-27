// Validation of KernelsByExpansion (formal_common.m) against Riemann-Roch spaces:
// on X1(45) (Case 45) and X1(63)/<8> (Case 63) over F2 -- the models of cartier_d8.m --
// compare, for every cusp P of degree <= 8 and every m <= 8/deg P, the coefficient space
// of DifferentialSpace(m*P) with the space obtained from local expansions (Expand).
// Run from this directory:  magma -n -b Case:=45 validate_residue_kernels.m
SetColumns(0);
SetMemoryLimit(16*10^9);
if not assigned Case then Case := "45"; end if;
load "independent/formal/formal_common.m";
pf := "independent/formal/logs/validate_residue_kernels_" cat Case cat ".progress";
Write(pf, "start" : Overwrite := true);
d := 8; t0 := Cputime(); F2 := GF(2);
P1<t> := FunctionField(F2);
b := t^2*(t-1)*(t^2-t+1); c := t^2*(t-1);
E0 := EllipticCurve([1-c, -b, -b, 0, 0]);
if Case eq "63" then
  psi := DivisionPolynomial(E0, 7); psi := psi/LeadingCoefficient(psi);
  F<u> := FunctionField(psi);
else
  psi := DivisionPolynomial(E0, 5); psi := psi/LeadingCoefficient(psi);
  F1<x5> := FunctionField(psi);
  a1,a2,a3,a4,a6 := Explode([F1!a : a in aInvariants(E0)]);
  S<Y> := PolynomialRing(F1);
  F<y5> := FunctionField(Y^2 + a1*x5*Y + a3*Y - (x5^3 + a2*x5^2 + a4*x5 + a6));
end if;
g := Genus(F);
cusps := DedupPlaces(Poles(F!jInvariant(E0)));
degs := [Degree(P) : P in cusps];
Say(pf, Sprintf("case %o genus %o cusp degrees %o", Case, g, Sort(degs)));
V, h := DifferentialSpace(DivisorGroup(F)!0);
small := [i : i in [1..#cusps] | degs[i] le d];
mmax := AssociativeArray(); for i in small do mmax[i] := d div degs[i]; end for;
t1 := Cputime();
WR := AllKernelsByExpansion(cusps, small, mmax, h, F2, g, t0, pf, "");
Say(pf, Sprintf("expansion kernels: %o spaces in %o s", #Keys(WR), Cputime(t1)));
t1 := Cputime();
WD := AllKernels(cusps, small, mmax, h, F2, g, t0, "");
Say(pf, Sprintf("Riemann-Roch kernels: %o spaces in %o s", #Keys(WD), Cputime(t1)));
nagree := #[k : k in Keys(WD) | WD[k] eq WR[k]];
Say(pf, Sprintf("RESULT case %o: expansion and Riemann-Roch kernel spaces agree in %o of %o cases", Case, nagree, #Keys(WD)));
Say(pf, Sprintf("VALIDATE_DONE case %o (%o s)", Case, Cputime(t0)));
quit;
