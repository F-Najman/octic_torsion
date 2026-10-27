SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
Attach("common/fast_full_aut_hash_m3_15.m");
Attach("C3x21/fast_full_aut_m3_21.m");
F<w>:=GF(4);X:=MDX11(3,21,F);C:=Curve(X);assert Genus(X) eq 25;
cs:=Cusps(X);assert {*Degree(c):c in cs*} eq {*1^^24,3^^8*};
xs:=NoncuspidalPlaces(X,4);assert #xs eq 144;
assert [#NoncuspidalPlaces(X,d):d in [1..3]] eq [0,0,0];
rs,labels,keys:=PlacesUpToFullAutM3TwentyOne(X,xs);
print "M3N21_D4_FULLAUT_SETUP",#xs,#rs,Cputime(t0);
K:=ext<F|4>;x:=rs[1];E:=EllipticCurve(X,x);Ls:=LevelStructure(X,x);P:=Ls`P;Q:=Ls`Q;
Q3:=7*Q;Q7:=3*Q;orb:={x};
for mat in M3FifteenAutMatrices() do
 a,b,c,d1:=Explode(mat);
 for u in [1,2,3] do
  La:=rec<recformat<P:PtEll,Q:PtEll>|P:=a*P+b*Q3,Q:=c*P+d1*Q3+5*u*Q7>;
  y:=ModuliPointFast(X,E,La);assert Degree(y) eq 4;Include(~orb,y);
 end for;
end for;
assert orb eq {xs[i]:i in [1..#xs]|labels[i] eq 1};
print "M3N21_D4_FULLAUT_REFERENCE_EQUAL",#orb,Cputime(t0);
q:=13;fast:=HeckeOperatorFast(X,q,x);ys:=[];ms:=[];
for phi in MDIsogenies(E,q) do
 La:=rec<recformat<P:PtEll,Q:PtEll>|P:=phi(Ls`P),Q:=phi(Ls`Q)>;
 y:=ModuliPoint3(X,Codomain(phi),La);ii:=Index(ys,y);
 if ii eq 0 then Append(~ys,y);Append(~ms,1);else ms[ii]+:=1;end if;
end for;
ref:=DivisorGroup(C)!0;
for i in [1..#ys] do
 mu:=ms[i]*4/Degree(ys[i]);assert Denominator(mu) eq 1;
 ref+:=Integers()!mu*Divisor(ys[i]);
end for;
assert ref eq fast;print "M3N21_D4_HECKE_REFERENCE_EQUAL",q,Degree(ref),Cputime(t0);
tests:=[];remaining:=[];
for i in [1..#rs] do
 x:=rs[i];passed:=true;
 for q in [13,19] do
  A:=HeckeOperatorFast(X,q,x)-Divisor(x)-q*Divisor(DiamondOperator(X,q,x));assert Degree(A) eq 0;
  pr:=IsPrincipal(A);Append(~tests,[i,q,pr select 1 else 0]);
  print "M3N21_D4_SIEVE_TEST",i,q,pr,Cputime(t0);
  if not pr then passed:=false;break;end if;
 end for;
 if passed then Append(~remaining,i);end if;
end for;
Write("C3x21/data/sieve_3_21_D4_arithmetic.txt",Sprint([*keys,labels,[Index(xs,x):x in rs],tests,remaining*]):Overwrite:=true);
printf "M3N21_D4_SIEVE_DONE tested=%o remaining=%o complete=true time=%o\n",#rs,#remaining,Cputime(t0);quit;
