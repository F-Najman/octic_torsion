/* Fresh complete octic full24-aut input, with direct original-inverse checks. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
Attach("common/fast_diamond_hash_m2.m");
Attach("C2x28/fast_full_aut_hash_m2_28.m");
X:=MDX11(2,28,GF(3));assert Genus(X) eq 25 and X`_E[1] eq 0 and X`_E[3] eq 0;
records:=[* *];coverage:=[* *];
for dd in [5..8] do
 xs:=NoncuspidalPlaces(X,dd);assert #xs eq [48,64,312,606][dd-4];
 rs,labels,keys:=PlacesUpToFullAutM2TwentyEight(X,xs);
 x:=rs[1];E:=EllipticCurve(X,x);Ls:=LevelStructure(X,x);P:=Ls`P;Q:=Ls`Q;orb:={x};
 for par in M2TwentyEightAutParameters() do
  b,c,d:=Explode(par);La:=rec<recformat<P:PtEll,Q:PtEll>|P:=P+14*b*Q,Q:=c*P+d*Q>;
  y:=ModuliPointFast(X,E,La);assert Degree(y) eq dd;Include(~orb,y);
 end for;
 assert orb eq {xs[j]:j in [1..#xs]|labels[j] eq 1};
 print "N28_D8_FULLAUT_REFERENCE_EQUAL",dd,#xs,#rs,#orb,Cputime(t0);
 K:=Parent(Evaluate(X`_coordinates[1],xs[1]));
 orbs:=[M2TwentyEightOrbitKeys(X,x,K):x in rs];
 Append(~coverage,[*dd,keys,labels,orbs,[Index(xs,x):x in rs]*]);
 for i in [1..#rs] do
  bc:=MDM2ToCommonField([Evaluate(f,rs[i]):f in X`_coordinates],K);
  assert MDPlaceCommonZeros(X`_coordinates,bc) eq rs[i];
  def:=[Integers()!a:a in Coefficients(DefiningPolynomial(K,GF(3)))];
  vals:=[[Integers()!a:a in Eltseq(z,GF(3))]:z in bc];Append(~records,[*dd,i,def,vals*]);
 end for;
 print "N28_D8_FULLAUT_EXPORT_TYPE",dd,#rs,Cputime(t0);
end for;
Write("C2x28/data/N28_D8_fullaut_Hecke_input.txt",Sprint([*records,coverage*]):Overwrite:=true);
print "N28_D8_FULLAUT_EXPORT_DONE",#records,Cputime(t0);quit;
