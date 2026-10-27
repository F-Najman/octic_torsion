SetColumns(0);SetQuitOnError(true);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
Attach("common/fast_full_aut_hash_m3_15.m");
Attach("C3x15/fast_full_aut_hash_m3_15_native.m");
Attach("C3x15/fast_full_aut_hash_m3_15_keytable.m");
X:=MDX11(3,15,GF(19));K:=GF(19^2);xs:=NoncuspidalPlaces(X,2);assert #xs eq 144;
for i in [1..#xs] do
 x:=xs[i];ky:=M3FifteenNativeKey([Evaluate(f,x):f in X`_coordinates],19,2);
 pair:=M3FifteenPairFromNativeKey(ky,K,19,2);
 y:=MDPlaceCommonZeros(X`_coordinates,pair);assert y eq x;
end for;
print "M3N15_NATIVE_KEY_ALL_QUADRATIC_ROUNDTRIPS",#xs,Cputime(t0);quit;
