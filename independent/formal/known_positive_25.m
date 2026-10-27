// Known-positive test of the rank-enumeration pipeline (EnumRanks + kernel spaces):
// X1(25) over F16 maps to X0(25) = P^1 with degree 10, so the fibre over a cusp of X0(25)
// whose preimages are rational cusps is a reduced cuspidal divisor D of degree 10 with
// h^0(D) >= 2, i.e. codim of Omega(-D) in H^0(Omega) is deg D + 1 - h^0(D) <= 9 < 10.
// We enumerate all reduced degree-10 divisors on the 18 rational cusps (C(18,10) = 43758)
// and expect the pipeline to report failures; each reported failure is re-checked by a
// direct RiemannRochSpace computation.  Both kernel methods (Riemann-Roch and expansions) are used.
// Run from this directory:  magma -n -b known_positive_25.m
SetColumns(0);
SetMemoryLimit(12*10^9);
load "independent/formal/formal_common.m";
pf := "independent/formal/logs/known_positive_25.progress";
Write(pf, "start" : Overwrite := true);
t0 := Cputime();
F16<z> := GF(16);
P1<t> := FunctionField(F16);
E5 := EllipticCurve([1-t, -t, -t, 0, 0]);
phi := XMultNumerator(E5, 5); phi := phi/LeadingCoefficient(phi);
F<xq> := FunctionField(phi);
g := Genus(F);
cusps := DedupPlaces(Poles(F!jInvariant(E5)));
degs := [Degree(P) : P in cusps];
rat := [i : i in [1..#cusps] | degs[i] eq 1];
Say(pf, Sprintf("genus %o, %o rational cusps", g, #rat));
V, h := DifferentialSpace(DivisorGroup(F)!0);
mmax := AssociativeArray(); for i in rat do mmax[i] := 1; end for;
WD := AllKernels(cusps, rat, mmax, h, F16, g, t0, "");
WE := AllKernelsByExpansion(cusps, rat, mmax, h, F16, g, t0, "", "");
Say(pf, Sprintf("kernel spaces agree: %o", &and[WD[k] eq WE[k] : k in Keys(WD)]));
// reduced divisors only: give each cusp "degree" 1 and allow multiplicity 1 by using weight 1 and
// restricting the enumeration to 0/1 multiplicities (rem div 1 could exceed 1, so use a wrapper)
procedure Enum01(j, rem, Kc, idx, W, g, target, cur, ~cnt, ~bad, ~fails)
  if rem eq 0 then
    cnt +:= 1;
    if g - Dimension(Kc) ne target then bad +:= 1; Append(~fails, cur); end if;
    return;
  end if;
  if j gt #idx or #idx - j + 1 lt rem then return; end if;
  Enum01(j+1, rem, Kc, idx, W, g, target, cur, ~cnt, ~bad, ~fails);
  Enum01(j+1, rem-1, Kc meet W[<idx[j],1>], idx, W, g, target, Append(cur, idx[j]), ~cnt, ~bad, ~fails);
end procedure;
for pair in [<"Riemann-Roch", WD>, <"expansion", WE>] do
  cnt := 0; bad := 0; fails := [];
  Enum01(1, 10, VectorSpace(F16, g), rat, pair[2], g, 10, [], ~cnt, ~bad, ~fails);
  Say(pf, Sprintf("%o kernels: reduced degree-10 divisors %o, failures (codim < 10) %o", pair[1], cnt, bad));
  for D in fails[1..Min(#fails, 5)] do
    Dv := &+[cusps[i] : i in D];
    Say(pf, Sprintf("  failure %o: direct h0(D) = %o, codim %o", D, Dimension(RiemannRochSpace(Dv)), g - Dimension(DifferentialSpace(Dv))));
  end for;
end for;
Say(pf, Sprintf("KNOWN_POSITIVE_DONE (%o s)", Cputime(t0)));
quit;
