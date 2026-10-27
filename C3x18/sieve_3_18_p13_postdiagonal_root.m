/* Post-diagonal finite checks only. Not a complete exclusion certificate. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
SetLogFile("C3x18/logs/sieve_3_18_p13_postdiagonal_root.log");
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
Attach("C3x18/fast_full_aut_hash_m3_18_root.m");
X:=MDX11(3,18,GF(13));C:=Curve(X);assert Genus(X) eq 10;
sel:=eval Read("C3x18/data/M3N18_p13_postdiagonal_selected_root.m");
out:=[* *];remaining:=[];reference_done:=false;
for dd in [3,4] do
 data:=eval Read(Sprintf("C3x18/data/M3N18_9multiple_aut27_p13_degree%o_data.m",dd));
 xs:=NoncuspidalPlaces(X,dd);K:=GF(13^dd);assert Coefficients(DefiningPolynomial(K)) eq data[3];
 assert #xs eq #data[7];
 seen:={};
 for index in sel[dd-2] do
  j:=data[6][index];x:=xs[j];ky:=M3EighteenRootKey([Evaluate(f,x):f in X`_coordinates],K,13,dd);
  assert ky eq data[7][j];
  if ky in seen then continue;end if;
  ims:=M3EighteenRootOrbitKeys(X,x,K);orb:=Seqset(ims);seen join:=orb;
  assert #orb in [9,27];
  if dd eq 3 and #orb eq 9 then
   Append(~out,[dd,index,j,#orb,0,0]);
   print "M3N18_P13_POSTDIAG_STABILIZER",dd,index,j,#orb,Cputime(t0);
   continue;
  end if;
  killed:=false;
  for q in [7,19,31] do
   D:=Divisor(x);HD:=HeckeOperatorFast(X,q,D);assert Degree(HD) eq (q+1)*dd;
   if not reference_done and q eq 7 then
    E:=EllipticCurve(X,x);LL:=LevelStructure(X,x);yy:=[];mm:=[];
    for phi in MDIsogenies(E,q) do
     La:=rec<recformat<P:PtEll,Q:PtEll>|P:=phi(LL`P),Q:=phi(LL`Q)>;
     y:=ModuliPoint3(X,Codomain(phi),La);h:=Index(yy,y);
     if h eq 0 then Append(~yy,y);Append(~mm,1);else mm[h]+:=1;end if;
    end for;
    ref:=DivisorGroup(C)!0;
    for h in [1..#yy] do mu:=mm[h]*dd/Degree(yy[h]);assert Denominator(mu) eq 1;ref+:=(Integers()!mu)*Divisor(yy[h]);end for;
    assert ref eq HD;reference_done:=true;
    print "M3N18_P13_POSTDIAG_GENERIC_REFERENCE",q,dd,Degree(HD),Cputime(t0);
   end if;
   A:=HD-D-q*DiamondOperator(X,q,D);assert Degree(A) eq 0;yes:=IsPrincipal(A);
   Append(~out,[dd,index,j,#orb,q,yes select 1 else 0]);
   print "M3N18_P13_POSTDIAG_HECKE",dd,index,j,#orb,q,yes,Cputime(t0);
   if not yes then killed:=true;break;end if;
  end for;
  if not killed then Append(~remaining,[dd,index,j]);end if;
 end for;
end for;
Write("C3x18/data/M3N18_p13_postdiagonal_arithmetic_root.txt",Sprint([*sel,out,remaining*]):Overwrite:=true);
print "M3N18_P13_POSTDIAG_COMPLETE",remaining,Cputime(t0);quit;
