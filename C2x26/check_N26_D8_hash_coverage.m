/* Fresh complete key coverage plus explicit moduli-inverse orbit checks. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
Attach("common/fast_diamond_hash_m2.m");
X:=MDX11(2,26,GF(3));assert X`_E[1] eq 0 and X`_E[3] eq 0;
for d in [5..8] do
 xs:=NoncuspidalPlaces(X,d);rs:=PlacesUpToDiamondM2HashDirect(X,xs);
 assert #xs eq [36,122,288,936][d-4] and #rs eq [6,26,48,159][d-4];
 K:=Parent(Evaluate(X`_coordinates[1],xs[1]));
 keys:=[MDM2FrobeniusKey([Evaluate(f,x):f in X`_coordinates],K,3,d):x in xs];
 allowed:=Seqset(keys);assert #allowed eq #xs;seen:=AssociativeArray();
 for j in [1..#rs] do
  for ky in MDM2DiamondKeysDirect(X,rs[j],K) do
   assert ky in allowed;
   if IsDefined(seen,ky) then assert seen[ky] eq j;else seen[ky]:=j;end if;
  end for;
 end for;
 assert #Keys(seen) eq #xs;
 x:=rs[1];actual:={};
 for a in [1,3,5,7,9,11] do
  y:=DiamondOperator(X,a,x);assert y in xs;
  Include(~actual,MDM2FrobeniusKey([Evaluate(f,y):f in X`_coordinates],K,3,d));
 end for;
 assert actual eq Seqset(MDM2DiamondKeysDirect(X,x,K));
 print "N26_D8_HASH_COVERAGE_VERIFIED",d,#xs,#rs,#actual,Cputime(t0);
end for;
print "N26_D8_HASH_COVERAGE_DONE complete=true",Cputime(t0);quit;
