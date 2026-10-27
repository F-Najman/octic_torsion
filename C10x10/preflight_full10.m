/* Exploratory exact model diagnostics; no point-classification claim. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
if not assigned p then p:=11;elif Type(p) eq MonStgElt then p:=StringToInteger(p);end if;
assert IsPrime(p) and p mod 5 eq 1 and GCD(p,30) eq 1;
F:=GF(p);A<t>:=FunctionField(F);P<x>:=PolynomialRing(A);
aa:=t^20-228*t^15+494*t^10+228*t^5+1;
bb:=t^30+522*t^25-10005*t^20-10005*t^10-522*t^5+1;
den:=t^5*(t^10+11*t^5-1)^5;
assert bb^2-aa^3 eq 1728*den;
j:=-aa^3/den;E:=EllipticCurve([A!0,0,0,-27*aa,54*bb]);
assert jInvariant(E) eq j;
f:=x^3-27*aa*x+54*bb;assert IsIrreducible(f);
B<r>:=ext<A|f>;PB<z>:=PolynomialRing(B);
ff:=z^2-(-3*r^2+108*(B!aa));assert IsIrreducible(ff);
K<s>:=ext<B|ff>;
assert Genus(K) eq 13;
print "FULL10_PREFLIGHT_GENUS",p,Degree(K,A),Genus(K),Cputime(t0);
cs:=Poles(K!j);assert #cs eq 36;
assert {*Degree(c):c in cs*} eq {*1^^36*};
assert &and[Valuation(K!j,c) eq -10:c in cs];
print "FULL10_PREFLIGHT_CUSPS",p,#cs,Sort([Degree(c):c in cs]),Cputime(t0);
J:=JH(100,[11]);assert Dimension(J) eq 13;
orders:=[];
for q in [11,31,41] do
 fp:=FrobeniusPolynomial(J,q);n:=Integers()!Evaluate(fp,1);Append(~orders,n);
 print "FULL10_PREFLIGHT_JAC_ORDER",q,n,Factorization(n),Cputime(t0);
end for;
print "FULL10_PREFLIGHT_GCD",GCD(orders),Cputime(t0);
print "FULL10_PREFLIGHT_DONE complete=false",Cputime(t0);quit;
