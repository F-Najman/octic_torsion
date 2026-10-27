/* Independent fresh full64-image export for every F9 target orbit. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
Attach("common/fast_diamond_hash_m2.m");
Attach("common/fast_full_aut_hash_m4_16.m");
Attach("common/fast_full_aut_hash_m4_16_relative.m");
a:=eval Read("C4x16/data/diagonal_4_16_D4_p5_p3_arithmetic.txt");
X:=MDX11(4,16,GF(9));F:=BaseRing(X);rows:=[* *];
for d in [3,4] do
 xs:=NoncuspidalPlaces(X,d);rs,labels,keys:=PlacesUpToFullAutM4SixteenRelative(X,xs);K:=ext<F|d>;
 assert [*d,keys,labels,[Index(xs,x):x in rs]*] eq a[13][d-2];
 orbs:=[M4SixteenRelativeOrbitKeys(X,x,K):x in rs];
 for i in [1..#rs] do assert Seqset(orbs[i]) eq {keys[j]:j in [1..#keys]|labels[j] eq i};end for;
 Append(~rows,[*d,keys,labels,orbs,[Index(xs,x):x in rs]*]);
 print "M4N16_D4_FULL64_ORBIT_EXPORT",d,#xs,#rs,Cputime(t0);
end for;
Write("C4x16/data/M4N16_D4_Aut64_coverage_arithmetic.txt",Sprint(rows):Overwrite:=true);
print "M4N16_D4_FULL64_COVERAGE_DONE",Cputime(t0);quit;
