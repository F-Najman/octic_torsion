/* Complete D4 support coverage: singles2/3/4, and ALL2+2 pairs inclrepetition. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
Attach("common/fast_diamond_hash_m2.m");
Attach("common/fast_full_aut_hash_m4_16.m");
Attach("common/fast_full_aut_hash_m4_16_relative.m");
Attach("C4x20/fast_full_aut_hash_m4_20.m");
X:=MDX11(4,20,GF(9));C:=Curve(X);old:=eval Read("C4x20/data/preflight_4_20_D4_Aut96_arithmetic.txt");
allxs:=[* *];allrs:=[* *];tests:=[* *];remaining:=[* *];
for d in [2..4] do
 xs:=NoncuspidalPlaces(X,d);rs,labels,keys:=PlacesUpToFullAutM4Twenty(X,xs);
 assert keys eq old[d-1][2] and labels eq old[d-1][3];Append(~allxs,xs);Append(~allrs,rs);
 qs:=d eq 2 select [29] else [13,17];
 for i in [1..#rs] do
  killed:=false;
  for q in qs do
   A:=HeckeOperatorFast(X,q,rs[i])-Divisor(rs[i])-q*Divisor(DiamondOperator(X,q,rs[i]));pr:=IsPrincipal(A);
   Append(~tests,[*1,d,i,0,q,pr select 1 else 0*]);print "M4N20_D4_COMPLETE_SINGLE",d,i,q,pr,Cputime(t0);
   if not pr then killed:=true;break;end if;
  end for;
  if not killed then Append(~remaining,[*1,d,i,0*]);end if;
 end for;
end for;
cache:=AssociativeArray();xs:=allxs[1];rs:=allrs[1];assert #xs eq 96 and #rs eq 2;
for i in [1..#rs] do
 a:=Index(xs,rs[i]);assert a eq old[1][5][i];
 for j in [1..#xs] do
  killed:=false;
  for q in [13,17] do
   for k in Seqset([a,j]) do
    ky:=[q,k];
    if not IsDefined(cache,ky) then cache[ky]:=HeckeOperatorFast(X,q,xs[k])-Divisor(xs[k])-q*Divisor(DiamondOperator(X,q,xs[k]));end if;
   end for;
   A:=cache[[q,a]]+cache[[q,j]];pr:=IsPrincipal(A);
   Append(~tests,[*2,2,i,j,q,pr select 1 else 0*]);print "M4N20_D4_COMPLETE_PAIR",i,j,q,pr,Cputime(t0);
   if not pr then killed:=true;break;end if;
  end for;
  if not killed then Append(~remaining,[*2,2,i,j*]);end if;
  Write("C4x20/data/sieve_4_20_D4_Aut96_arithmetic.txt",Sprint([*tests,remaining*]):Overwrite:=true);
 end for;
end for;
print "M4N20_D4_COMPLETE_SIEVE_DONE singles=18 pairs=192 remaining",remaining,"complete=true",Cputime(t0);quit;
