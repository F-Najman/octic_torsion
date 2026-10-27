/* Complete D4 support coverage: singles2/3/4, and ALL2+2 pairs inclrepetition. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
Attach("common/fast_diamond_hash_m2.m");
Attach("common/fast_full_aut_hash_m4_16.m");
Attach("common/fast_full_aut_hash_m4_16_relative.m");
Attach("C4x20/fast_full_aut_hash_m4_20.m");
X:=MDX11(4,20,GF(9));C:=Curve(X);old:=eval Read("C4x20/data/preflight_4_20_D4_Aut96_arithmetic.txt");
records:=[**];
for d in [2,4] do
 xs:=NoncuspidalPlaces(X,d);rs,labels,keys:=PlacesUpToFullAutM4Twenty(X,xs);
 assert keys eq old[d-1][2] and labels eq old[d-1][3];K:=ext<GF(9)|d>;
 if d eq 4 then
  for i in [4,12] do
   a:=Index(xs,rs[i]);orbitkeys:=M4TwentyOrbitKeys(X,rs[i],K);
   stabilizer:=[j:j in [1..96]|orbitkeys[j] eq keys[a]];assert #stabilizer eq 8;
   Append(~records,[*1,d,i,a,keys[a],orbitkeys,stabilizer*]);
   print "M4N20_ROOT_QUARTIC_STABILIZER",i,#stabilizer,Cputime(t0);
  end for;
 else
  for pair in [[1,82],[2,84]] do
   i,j:=Explode(pair);a:=Index(xs,rs[i]);sourcekeys:=Sort([keys[a],keys[j]]);
   ka:=M4TwentyOrbitKeys(X,xs[a],K);kj:=M4TwentyOrbitKeys(X,xs[j],K);
   ims:=[Sort([ka[g],kj[g]]):g in [1..96]];
   stabilizer:=[g:g in [1..96]|ims[g] eq sourcekeys];
   Append(~records,[*2,d,i,j,sourcekeys,ims,stabilizer*]);
   print "M4N20_ROOT_PAIR_STABILIZER",i,j,#stabilizer,Cputime(t0);
  end for;
 end if;
end for;
Write("C4x20/data/M4N20_stabilizers_root_arithmetic.txt",Sprint([*M4TwentyAutParameters(),records*]):Overwrite:=true);
print "M4N20_ROOT_STABILIZERS_COMPLETE",Cputime(t0);quit;
