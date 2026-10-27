/* Fresh third-prime test of portable survivors from the two-prime runs. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
if Type(Records) eq MonStgElt then Records:=eval Records;end if;
assert #Records gt 0 and #Seqset(Records) eq #Records;
X:=MDX11(2,32,GF(3));assert Genus(X) eq 37 and X`_E[1] eq 0 and X`_E[3] eq 0;
a:=eval Read("C2x32/data/N32_D8_fullaut_Hecke_input.txt");records:=a[1];
R:=PolynomialRing(GF(3));tests:=[* *];
for n in Records do
 d,i,def,vals:=Explode(records[n]);assert d in [6,7,8];
 f:=R!def;assert IsIrreducible(f) and Degree(f) eq d;
 K<z>:=ext<GF(3)|f>;bc:=[Evaluate(R!v,z):v in vals];
 x:=MDPlaceCommonZeros(X`_coordinates,bc);assert Degree(x) eq d;
 q:=d eq 6 select 11 else 5;A:=HeckeOperatorFast(X,q,x)-Divisor(x)-q*Divisor(DiamondOperator(X,q,x));assert Degree(A) eq 0;
 pr:=IsPrincipal(A);Append(~tests,[*n,d,i,q,pr*]);
 print "N32_D8_FULLAUT_TIMING_HECKE_CASE",n,d,i,q,pr,Cputime(t0);
end for;
Write(Output,Sprint(tests):Overwrite:=true);
print "N32_D8_FULLAUT_TIMING_HECKE_DONE",#tests,Cputime(t0);quit;
