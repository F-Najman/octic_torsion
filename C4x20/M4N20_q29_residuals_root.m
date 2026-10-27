/* Complete D4 support coverage: singles2/3/4, and ALL2+2 pairs inclrepetition. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
Attach("common/fast_diamond_hash_m2.m");
Attach("common/fast_full_aut_hash_m4_16.m");
Attach("common/fast_full_aut_hash_m4_16_relative.m");
Attach("C4x20/fast_full_aut_hash_m4_20.m");
X:=MDX11(4,20,GF(9));C:=Curve(X);old:=eval Read("C4x20/data/preflight_4_20_D4_Aut96_arithmetic.txt");
tests:=[**];
for d in [2..4] do
 xs:=NoncuspidalPlaces(X,d);rs,labels,keys:=PlacesUpToFullAutM4Twenty(X,xs);
 assert keys eq old[d-1][2] and labels eq old[d-1][3];
 indices:=d eq 4 select [4,12] else [1..#rs];
 for i in indices do
  A:=HeckeOperatorFast(X,29,rs[i])-Divisor(rs[i])-29*Divisor(DiamondOperator(X,29,rs[i]));
  principal:=IsPrincipal(A);Append(~tests,[*d,i,0,29,principal select 1 else 0*]);
  print "M4N20_ROOT_Q29_SINGLE",d,i,principal,Cputime(t0);
 end for;
end for;
xs:=NoncuspidalPlaces(X,2);a:=old[1][5][1];j:=82;
A:=&+[HeckeOperatorFast(X,29,xs[k])-Divisor(xs[k])-29*Divisor(DiamondOperator(X,29,xs[k])):k in [a,j]];
principal:=IsPrincipal(A);Append(~tests,[*2,1,j,29,principal select 1 else 0*]);
print "M4N20_ROOT_Q29_PAIR",1,j,principal,Cputime(t0);
Write("C4x20/data/M4N20_q29_residuals_root_arithmetic.txt",Sprint(tests):Overwrite:=true);
print "M4N20_ROOT_Q29_DONE",Cputime(t0);quit;
