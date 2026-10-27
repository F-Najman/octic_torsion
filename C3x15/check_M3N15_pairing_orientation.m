SetColumns(0);SetQuitOnError(true);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
for pp in [7,13] do
 F:=GF(pp);z:=RootOfUnity(3,F);assert z^2+z+1 eq 0;
 X:=MDX11(3,15,F:zeta_M:=z);xs:=NoncuspidalPlaces(X,2);assert #xs gt 0;
 E:=EllipticCurve(X,xs[1]);L:=LevelStructure(X,xs[1]);
 e:=WeilPairing(L`P,5*L`Q,3);print "M3N15_PAIRING_ORIENTATION",pp,z,F!e,e eq z,e eq z^2,Cputime(t0);
end for;quit;
