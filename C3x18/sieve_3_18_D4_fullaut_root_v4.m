/* Fresh complete27-automorphism partition and first split-Hecke timings. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
Attach("C3x18/fast_full_aut_hash_m3_18_root.m");
X:=MDX11(3,18,GF(7));C:=Curve(X);exports:=[* *];all_ps:=[* *];all_rs:=[* *];remaining:=[* *];tests:=[* *];
for dd in [2..4] do
 xs:=NoncuspidalPlaces(X,dd);rs,labels:=PlacesUpToFullAutM3EighteenRootHash(X,xs);
 print "M3N18_ROOT_AUT_PARTITION",dd,#xs,#rs,Cputime(tt);
 x:=rs[1];E:=EllipticCurve(X,x);Ls:=LevelStructure(X,x);orb:={x};P:=Ls`P;Q:=Ls`Q;
 for par in M3EighteenRootAutParameters() do
  a,b,c,d:=Explode(par);La:=rec<recformat<P:PtEll,Q:PtEll>|P:=a*P+6*b*Q,Q:=c*P+d*Q>;
  y:=ModuliPointFast(X,E,La);assert Degree(y) eq dd;Include(~orb,y);
 end for;
 assert orb eq {xs[i]:i in [1..#xs]|labels[i] eq 1};
 print "M3N18_ROOT_AUT_NATIVE_REFERENCE",dd,#orb,Cputime(tt);
 Append(~all_ps,xs);Append(~all_rs,rs);
 for i in [1..#rs] do
  D:=Divisor(rs[i]);killed:=false;
  for q in [13,19] do
   HD:=HeckeOperatorFast(X,q,D);assert Degree(HD) eq (q+1)*dd;
   if dd eq 2 and i eq 1 and q eq 13 then
    EE:=EllipticCurve(X,rs[i]);LL:=LevelStructure(X,rs[i]);yy:=[];mm:=[];
    for phi in MDIsogenies(EE,q) do
     La:=rec<recformat<P:PtEll,Q:PtEll>|P:=phi(LL`P),Q:=phi(LL`Q)>;
     y:=ModuliPoint3(X,Codomain(phi),La);j:=Index(yy,y);
     if j eq 0 then Append(~yy,y);Append(~mm,1);else mm[j]+:=1;end if;
    end for;
    ref:=DivisorGroup(C)!0;
    for j in [1..#yy] do mu:=mm[j]*dd/Degree(yy[j]);assert Denominator(mu) eq 1;ref+:=(Integers()!mu)*Divisor(yy[j]);end for;
    assert ref eq HD;
    print "M3N18_ROOT_HECKE_GENERIC_REFERENCE",dd,q,Degree(HD),Cputime(tt);
   end if;
   A:=HD-D-q*DiamondOperator(X,q,D);assert Degree(A) eq 0;yes:=IsPrincipal(A);
   Append(~tests,[dd,i,0,q,yes select 1 else 0]);print "M3N18_ROOT_SINGLE_HECKE",dd,i,q,yes,Cputime(tt);
   if not yes then killed:=true;break;end if;
  end for;
  if not killed then Append(~remaining,[dd,i,0]);end if;
 end for;
 KK:=GF(7^dd);keys:=[M3EighteenRootKey([Evaluate(f,x):f in X`_coordinates],KK,7,dd):x in xs];
 Append(~exports,[*dd,#xs,[Index(xs,r):r in rs],labels,keys*]);
end for;
// Complete coverage of2+2, including repetition. The27 automorphisms
// act transitively on27 quadratics, so fixing the first leaves27 choices.
xs:=all_ps[1];rs:=all_rs[1];assert #rs eq 1 and #xs eq 27;
for i in [1..#xs] do
 D:=Divisor(rs[1])+Divisor(xs[i]);assert Degree(D) eq 4;killed:=false;
 for q in [13,19] do
  A:=HeckeOperatorFast(X,q,D)-D-q*DiamondOperator(X,q,D);assert Degree(A) eq 0;yes:=IsPrincipal(A);
  Append(~tests,[22,1,i,q,yes select 1 else 0]);print "M3N18_ROOT_PAIR_HECKE",i,q,yes,Cputime(tt);
  if not yes then killed:=true;break;end if;
 end for;
 if not killed then Append(~remaining,[22,1,i]);end if;
end for;
Write("C3x18/data/sieve_3_18_D4_fullaut_root_v4_arithmetic.txt",Sprint([*exports,tests,remaining*]):Overwrite:=true);
print "M3N18_ROOT_AUT_SIEVE_DONE",remaining,"complete",true,Cputime(tt);quit;
