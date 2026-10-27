/* Complete D4 support coverage: singles2/3/4, and ALL2+2 pairs inclrepetition. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
Attach("common/fast_diamond_hash_m2.m");
Attach("common/fast_full_aut_hash_m4_16.m");
Attach("common/fast_full_aut_hash_m4_16_relative.m");
Attach("C4x20/fast_full_aut_hash_m4_20.m");
X:=MDX11(4,20,GF(9));C:=Curve(X);old:=eval Read("C4x20/data/preflight_4_20_D4_Aut96_arithmetic.txt");
xs:=NoncuspidalPlaces(X,2);rs,labels,keys:=PlacesUpToFullAutM4Twenty(X,xs);
assert keys eq old[1][2] and labels eq old[1][3];cache:=AssociativeArray();tests:=[**];remaining:=[**];
for i in [1..#rs] do
 a:=Index(xs,rs[i]);assert a eq old[1][5][i];
 start:=i eq 1 select 83 else 1;
 for j in [start..#xs] do
  killed:=false;
  for q in [13,17] do
   for k in Seqset([a,j]) do
    ky:=[q,k];if not IsDefined(cache,ky) then cache[ky]:=HeckeOperatorFast(X,q,xs[k])-Divisor(xs[k])-q*Divisor(DiamondOperator(X,q,xs[k]));end if;
   end for;
   pr:=IsPrincipal(cache[[q,a]]+cache[[q,j]]);
   Append(~tests,[*i,j,q,pr select 1 else 0*]);print "M4N20_ROOT_PAIR_TAIL",i,j,q,pr,Cputime(t0);
   if not pr then killed:=true;break;end if;
  end for;
  if not killed then Append(~remaining,[*i,j*]);end if;
  Write("C4x20/data/M4N20_pair_tail_root_arithmetic.txt",Sprint([*tests,remaining*]):Overwrite:=true);
 end for;
end for;
print "M4N20_ROOT_PAIR_TAIL_DONE",#tests,remaining,Cputime(t0);quit;
