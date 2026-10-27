SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);t:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");
Attach("common/fast_diamond_hash_m2.m");
Attach("common/fast_diamond_hash_cyclic_bc.m");
X:=MDX1(40,GF(3));S:=NoncuspidalPlaces(X,4);assert #S eq 32;
reps:=PlacesUpToDiamondCyclicHashBC(X,S);
allorbits:=[];
for P in reps do
 orb:=Seqset([DiamondOperator(X,a,P):a in [1..20]|GCD(a,40) eq 1]);
 assert #orb eq 8 and orb subset Seqset(S);
 assert &and[orb meet other eq {}:other in allorbits];Append(~allorbits,orb);
end for;
assert &join allorbits eq Seqset(S) and #reps eq 4;
printf "C40_D8_DIAMOND_HASH_REFERENCE_DONE places=32 orbits=4 exact_original_diamond_partition=true time=%o\n",Cputime(t);
quit;
