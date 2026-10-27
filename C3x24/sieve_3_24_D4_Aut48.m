/* Fresh full relative-degree4 sieve, with all rational cusp complements. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
X:=MDX11(3,24,GF(7));assert Genus(X) eq 25;
a:=eval Read("C3x24/data/preflight_3_24_D4_Aut48_input.txt");records:=a[1];assert #records eq 14;
R:=PolynomialRing(GF(7));tests:=[* *];remaining:=[];
for n in [1..#records] do
 d,i,def,vals:=Explode(records[n]);assert d in [3,4];f:=R!def;assert IsIrreducible(f) and Degree(f) eq d;
 K<z>:=ext<GF(7)|f>;qt:=[Evaluate(R!v,z):v in vals];x:=MDPlaceCommonZeros(X`_coordinates,qt);assert Degree(x) eq d;killed:=false;
 for q in [13,19] do
  T:=HeckeOperatorFast(X,q,x);A:=T-Divisor(x)-q*Divisor(DiamondOperator(X,q,x));assert Degree(A) eq 0;
  pr:=IsPrincipal(A);Append(~tests,[*n,d,i,q,pr select 1 else 0*]);print "M3N24_D4_FULLAUT_HECKE_TEST",n,d,i,q,pr,Cputime(t0);
  if not pr then killed:=true;break;end if;
 end for;
 if not killed then Append(~remaining,n);end if;
 Write("C3x24/data/sieve_3_24_D4_Aut48_arithmetic.txt",Sprint([*tests,remaining*]):Overwrite:=true);
end for;
print "M3N24_D4_FULLAUT_HECKE_DONE tested=14 remaining",remaining,"complete=true",Cputime(t0);quit;
