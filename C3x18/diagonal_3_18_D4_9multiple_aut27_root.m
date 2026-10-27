/* Fresh exact nine-multiple cusp-diagonal inputs; no completion claim. */
SetColumns(0);SetQuitOnError(true);SetLogFile("C3x18/logs/diagonal_3_18_D4_9multiple_aut27_root.log");SetMemoryLimit(12*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
load "common/class_group_direct_RR_checked_v3.m";
Attach("common/fast_diamond_hash_m2.m");
Attach("C3x18/fast_full_aut_hash_m3_18_root.m");
for p in [13] do
 X:=MDX11(3,18,GF(p));C:=Curve(X);FF:=FunctionField(C);F0,afmap:=AlgorithmicFunctionField(FF);assert Genus(X) eq 10;
 E:=EllipticCurve([FF!f:f in X`_E]);P:=E![FF!f:f in X`_P];Q:=E![FF!f:f in X`_Q];
 ai:=aInvariants(E);fs:=[ai[3],ai[1]^3-27*ai[3]];
 for k in [1..9] do
  for j in [0..2] do
   A0:=k*Q+j*P;assert A0 ne E!0 and A0[1] ne 0 and A0[2] ne 0;fs cat:=[A0[1],A0[2]];
  end for;
  A0:=k*Q;B0:=A0+P;assert A0 ne B0 and A0 ne -B0;Append(~fs,A0[1]-B0[1]);
 end for;
 cs:=Cusps(X);assert {*Degree(c):c in cs*} eq {*1^^18,3^^6*};groups:=AssociativeArray();
 for c in cs do
  sig:=[Valuation(f,c):f in fs];if IsDefined(groups,sig) then Append(~groups[sig],c);else groups[sig]:=[c];end if;
 end for;
 sigs:=Sort(Setseq(Keys(groups)));assert #sigs eq 24;blocks:=[&+[Divisor(c):c in groups[s]]:s in sigs];ds:=[Degree(b):b in blocks];assert {*d:d in ds*} eq {*1^^18,3^^6*};
 c0:=Support(blocks[Index(ds,1)])[1];
 print "M3N18_9MULT_CUSP_PROFILES",p,#sigs,ds,Cputime(t0);
 G,unused,phi:=ClassGroup(F0:ReductionDivisor:=FunctionFieldDivisor(Divisor(c0)),Proof:=true);
 print "M3N18_9MULT_NATIVE_GROUP_PROVEN",p,Invariants(G),Cputime(t0);
 GG,gm,mg:=MDClassGroupDirectRRCheckedV3(F0,Sprintf("C3x18/data/M3N18_9multiple_aut27_p%o_direct_RR_identities.txt",p));assert GG eq G;
 m2:=map<DivisorGroup(C)->G|D:->mg(FunctionFieldDivisor(D))>;
 gs:=[m2(blocks[i]-ds[i]*Divisor(c0)):i in [1..#blocks]];H:=sub<G|gs>;Qs:=quo<G|H>;
 print "M3N18_9MULT_CUSP_GROUP",p,Invariants(H),Invariants(Qs),Cputime(t0);
 J:=JH(54,[19]);fp:=FrobeniusPolynomial(J,p);assert &*[n:n in Invariants(G)|n ne 0] eq Integers()!Evaluate(fp,1);
 rat:=[c:c in cs|Degree(c) eq 1];ratcl:=[m2(Divisor(c)-Divisor(c0)):c in rat];
 data:=[*p,sigs,ds,Invariants(G),[Eltseq(v):v in Rows(RelationMatrix(G))],[Eltseq(v):v in gs],Invariants(H),[Eltseq(v):v in ratcl],Coefficients(fp)*];
 Write(Sprintf("C3x18/data/M3N18_9multiple_aut27_p%o_group_data.m",p),Sprint(data):Overwrite:=true);
 for d in [2..4] do
  allxs:=NoncuspidalPlaces(X,d);xs:=allxs;KS:=GF(p^d);labels:=[i:i in [1..#xs]];images:=[* *];
  if d eq 4 then xs,labels:=PlacesUpToFullAutM3EighteenRootHash(X,allxs);images:=[*M3EighteenRootOrbitKeys(X,x,KS):x in xs*];end if;
  indices:=[Index(allxs,x):x in xs];keys:=[M3EighteenRootKey([Evaluate(f,x):f in X`_coordinates],KS,p,d):x in allxs];
  vals:=[];portable:=[];print "M3N18_9MULT_SUPPORT_COUNT",p,d,#xs,Cputime(t0);
  for i in [1..#xs] do
   cl:=m2(Divisor(xs[i])-d*Divisor(c0));Append(~vals,Eltseq(cl));
   vs:=MDM2ToCommonField([Evaluate(f,xs[i]):f in X`_coordinates],KS);Append(~portable,[Eltseq(v):v in vs]);
   if i mod 100 eq 0 then print "M3N18_9MULT_CLASSES_PROGRESS",p,d,i,#xs,Cputime(t0);end if;
  end for;
  Write(Sprintf("C3x18/data/M3N18_9multiple_aut27_p%o_degree%o_data.m",p,d),Sprint([*p,d,Coefficients(DefiningPolynomial(KS)),vals,portable,indices,keys,labels,images*]):Overwrite:=true);
  print "M3N18_9MULT_SUPPORT_COMPLETE",p,d,#xs,Cputime(t0);
 end for;
end for;
print "M3N18_9MULT_AUT27_INPUTS_COMPLETE",Cputime(t0);
quit;
