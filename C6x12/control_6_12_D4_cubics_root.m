SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
SetLogFile("C6x12/logs/control_6_12_D4_cubics_root.log");
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
Attach("C6x12/model_6_12_root.m");
Attach("common/fast_diamond_hash_m2.m");
X:=MDX612(GF(7));C:=Curve(X);J:=JH(72,[13]);assert Genus(X) eq 9 and Dimension(J) eq 9;
frob:=[* *];
for p in [7,13] do
 f:=FrobeniusPolynomial(J,p);Append(~frob,[*p,Coefficients(f)*]);
end for;
xs:=NoncuspidalPlaces(X,3);rs:=Representatives612(X,xs);assert #xs eq 96 and #rs eq 48;
K:=GF(7^3);records:=[* *];keys:=[];
for x in xs do
 coords:=MDM2ToCommonField([Evaluate(f,x):f in X`_coordinates],K);
 Append(~keys,MDM2FrobeniusKey(coords,K,7,3));
end for;
for i in [1..#rs] do
 x:=rs[i];E:=EllipticCurve(X,x);LL:=LevelStructure612(X,x);HD:=HeckeOperator612(X,13,x);
 if i eq 1 then
  yy:=[];mm:=[];
  for phi in MDIsogenies(E,13) do
   La:=rec<recformat<P:PtEll,Q:PtEll>|P:=phi(LL`P),Q:=phi(LL`Q)>;
   y:=ModuliPoint612(X,Codomain(phi),La);j:=Index(yy,y);
   if j eq 0 then Append(~yy,y);Append(~mm,1);else mm[j]+:=1;end if;
  end for;
  ref:=DivisorGroup(C)!0;
  for j in [1..#yy] do mu:=mm[j]*3/Degree(yy[j]);assert Denominator(mu) eq 1;ref+:=(Integers()!mu)*Divisor(yy[j]);end for;
  assert ref eq HD;print "M6N12_CUBIC_GENERIC_REFERENCE",Degree(HD),Cputime(t0);
 end if;
 A:=HD-Divisor(x)-13*Divisor(Diamond612(X,13,x));assert Degree(A) eq 0 and not IsPrincipal(A);
 assert ModuliPoint612(X,E,LL) eq x;
 vals:=aInvariants(E) cat Eltseq(LL`P)[1..2] cat Eltseq(LL`Q)[1..2] cat [Evaluate(f,x):f in X`_coordinates];
 vals:=MDM2ToCommonField(vals,K);Append(~records,[*i,[Eltseq(v):v in vals]*]);
 print "M6N12_CUBIC_EXCLUDED",i,Cputime(t0);
end for;
Write("C6x12/data/M6N12_cubic_control_arithmetic_root.txt",Sprint([*frob,Coefficients(DefiningPolynomial(K)),keys,records*]):Overwrite:=true);
print "M6N12_CUBIC_CONTROL_COMPLETE",Cputime(t0);quit;
