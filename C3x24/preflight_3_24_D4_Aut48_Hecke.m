SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
Attach("common/fast_full_aut_hash_m3_15.m");
Attach("C3x24/fast_full_aut_hash_m3_24.m");
X:=MDX11(3,24,GF(7));assert Genus(X) eq 25;records:=[* *];coverage:=[* *];tests:=[* *];
print "M3N24_D4_PREFLIGHT_MODEL",Cputime(t0);
for d in [3,4] do
 xs:=NoncuspidalPlaces(X,d);assert #xs eq (d eq 3 select 48 else 456);
 rs,labels,keys:=PlacesUpToFullAutM3TwentyFour(X,xs);print "M3N24_D4_AUT48_PARTITION",d,#xs,#rs,Cputime(t0);
 x:=rs[1];E:=EllipticCurve(X,x);Ls:=LevelStructure(X,x);P:=Ls`P;Q:=Ls`Q;R3:=8*Q;R8:=3*Q;orb:={x};
 for mat in M3FifteenAutMatrices() do
  a,b,c,dd:=Explode(mat);
  for u in [1,3] do
   La:=rec<recformat<P:PtEll,Q:PtEll>|P:=a*P+b*R3,Q:=2*(c*P+dd*R3)+3*u*R8>;
   Include(~orb,ModuliPointFast(X,E,La));
  end for;
 end for;
 assert orb eq {xs[j]:j in [1..#xs]|labels[j] eq 1};print "M3N24_D4_AUT48_REFERENCE_EQUAL",d,#orb,Cputime(t0);
 K:=GF(7^d);orbs:=[M3TwentyFourOrbitKeys(X,x,K):x in rs];
 Append(~coverage,[*d,keys,labels,orbs,[Index(xs,x):x in rs]*]);
 for i in [1..#rs] do
  qt:=M3FifteenToCommon([Evaluate(f,rs[i]):f in X`_coordinates],K);assert MDPlaceCommonZeros(X`_coordinates,qt) eq rs[i];
  Append(~records,[*d,i,[Integers()!a:a in Coefficients(DefiningPolynomial(K,GF(7)))],[[Integers()!a:a in Eltseq(z,GF(7))]:z in qt]*]);
 end for;
 Write("C3x24/data/preflight_3_24_D4_Aut48_input.txt",Sprint([*records,coverage*]):Overwrite:=true);
 for i in [1..Minimum(2,#rs)] do
  A:=HeckeOperatorFast(X,13,rs[i])-Divisor(rs[i])-13*Divisor(DiamondOperator(X,13,rs[i]));
  pr:=IsPrincipal(A);Append(~tests,[*d,i,13,pr select 1 else 0*]);print "M3N24_D4_HECkE_PROBE",d,i,13,pr,Cputime(t0);
 end for;
end for;
Write("C3x24/data/preflight_3_24_D4_Aut48_Hecke_arithmetic.txt",Sprint(tests):Overwrite:=true);
print "M3N24_D4_AUT48_PROBE_DONE complete=true",Cputime(t0);quit;
