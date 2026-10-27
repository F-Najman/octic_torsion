/* Fresh octic single-support sieve; all possible cusp complements covered. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
if Type(Shard) eq MonStgElt then Shard:=StringToInteger(Shard);end if;assert Shard in [0,1];
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
X:=MDX11(2,32,GF(3));assert Genus(X) eq 37 and X`_E[1] eq 0 and X`_E[3] eq 0;
a:=eval Read("C2x32/data/N32_D8_fullaut_Hecke_input.txt");records:=a[1];assert #records eq 50;
R:=PolynomialRing(GF(3));tested:=0;remaining:=[];
for n in [1..#records] do
 if (n-1) mod 2 ne Shard then continue;end if;
 d,i,def,vals:=Explode(records[n]);assert d in [6,7,8];f:=R!def;assert IsIrreducible(f) and Degree(f) eq d;
 K<z>:=ext<GF(3)|f>;bc:=[Evaluate(R!v,z):v in vals];
 x:=MDPlaceCommonZeros(X`_coordinates,bc);assert Degree(x) eq d;
 // Exact Tate labels: q11/19 annihilate every complement of degree2;
 // q5/11 annihilate every complement of degree<=1.
 qs:=d eq 6 select [11,19] else [5,11];killed:=false;
 for q in qs do
  A:=HeckeOperatorFast(X,q,x)-Divisor(x)-q*Divisor(DiamondOperator(X,q,x));assert Degree(A) eq 0;
  pr:=IsPrincipal(A);print "N32_D8_FULLAUT_HIGH_HECKE_CASE",Shard,n,d,i,q,pr,Cputime(t0);
  if not pr then killed:=true;break;end if;
 end for;
 tested+:=1;if not killed then Append(~remaining,n);end if;
end for;
assert tested eq 25;
print "N32_D8_FULLAUT_HIGH_HECKE_SURVIVORS",Shard,remaining;
printf "N32_D8_FULLAUT_HIGH_HECKE_DONE shard=%o tested=%o remaining=%o complete=true time=%o\n",Shard,tested,#remaining,Cputime(t0);quit;
