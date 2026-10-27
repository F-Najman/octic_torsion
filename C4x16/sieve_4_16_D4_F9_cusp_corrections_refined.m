/* Fresh relativeD4 correction sieve at inert3/residueF9. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
Attach("common/fast_diamond_hash_m2.m");
Attach("common/fast_full_aut_hash_m4_16.m");
Attach("common/fast_full_aut_hash_m4_16_relative.m");
function M4RationalCuspPairs(X)
 FF:=FunctionField(Curve(X));E:=EllipticCurve([FF!f:f in X`_E]);P:=E![FF!f:f in X`_P];Q:=E![FF!f:f in X`_Q];
 ai:=aInvariants(E);B:=ai[2];assert ai eq [FF!1,B,B,FF!0,FF!0] and P eq E![0,0];
 assert Discriminant(E) eq B^4*(1-16*B);fs:=[B,1-16*B];
 for k in [1..8] do for j in [0..3] do
  A:=k*Q+j*P;assert A ne E!0 and A[1] ne 0 and A[2] ne 0 and A[1]+B ne 0;
  fs cat:=[A[1],A[2],A[1]+B,A[2]+A[1]+B];
 end for;end for;
 for k in [1..8] do A:=k*Q;for j in [1..3] do R:=A+j*P;assert A ne R and A ne -R;Append(~fs,A[1]-R[1]);end for;end for;
 rawcs:=Cusps(X);groups:=AssociativeArray();
 for c in rawcs do sig:=[Valuation(f,c):f in fs];
  if IsDefined(groups,sig) then Append(~groups[sig],c);else groups[sig]:=[c];end if;
 end for;
 sigs:=Sort(Setseq(Keys(groups)));assert #sigs eq 24;
 pairs:=[groups[sig]:sig in sigs|#groups[sig] eq 2 and &and[Degree(c) eq 1:c in groups[sig]]];
 assert #pairs eq 4;
 fsd:=[B,1-16*B];Qd:=5*Q;
 for k in [1..8] do for j in [0..3] do A:=k*Qd+j*P;fsd cat:=[A[1],A[2],A[1]+B,A[2]+A[1]+B];end for;end for;
 for k in [1..8] do A:=k*Qd;for j in [1..3] do R:=A+j*P;Append(~fsd,A[1]-R[1]);end for;end for;
 ps:=[[Valuation(f,pair[1]):f in fs]:pair in pairs];perm:=[];
 for pair in pairs do
  ss:=[Valuation(f,pair[1]):f in fsd];assert ss eq [Valuation(f,pair[2]):f in fsd];
  jj:=Index(ps,ss);assert jj ne 0;Append(~perm,jj);
 end for;assert Seqset(perm) eq {1,2,3,4};
 return pairs,perm;
end function;
X:=MDX11(4,16,GF(9));C:=Curve(X);pairs,perm:=M4RationalCuspPairs(X);V:=[DivisorGroup(C)|];
for pair in pairs do v:=Divisor(pair[2])-Divisor(pair[1]);Append(~V,v);Append(~V,-v);end for;
q:=5;corr:=[DivisorGroup(C)!0];
for vi in [1..8] do for wi in [1..8] do if (wi-1) div 2+1 eq perm[(vi-1) div 2+1] then Append(~corr,V[vi]+q*V[wi]);end if;end for;end for;
assert #corr eq 17 and #Seqset(corr) eq 17;print "M4N16_D4_CUSP_PAIR_DIAMOND5_PERMUTATION",perm;
xs:=NoncuspidalPlaces(X,3);rs,labels,keys:=PlacesUpToFullAutM4SixteenRelative(X,xs);assert #rs eq 5 and #xs eq 320;
print "M4N16_D4_F9_CUBIC_CORRECTION_SETUP",#pairs,#corr,#rs,Cputime(t0);
tests:=[* *];remaining:=[];
for i in [1..5] do
 x:=rs[i];A:=HeckeOperatorFast(X,q,x)-Divisor(x)-q*Divisor(DiamondOperator(X,q,x));assert Degree(A) eq 0;
 flags:=[];for cc in corr do Append(~flags,IsPrincipal(A+cc) select 1 else 0);end for;
 Append(~tests,[*i,flags*]);if &+flags gt 0 then Append(~remaining,i);end if;
 print "M4N16_D4_F9_CUBIC_CORRECTION_TEST",i,&+flags,Cputime(t0);
end for;
Write("C4x16/data/sieve_4_16_D4_F9_cusp_corrections_refined_arithmetic.txt",Sprint([*keys,labels,[Index(xs,x):x in rs],tests,remaining,perm*]):Overwrite:=true);
print "M4N16_D4_F9_CUBIC_CORRECTIONS_DONE",remaining,Cputime(t0);quit;
