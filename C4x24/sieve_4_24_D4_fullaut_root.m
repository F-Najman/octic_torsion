/* Bounded fresh p5 full64-automorphism+split-Hecke probe. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
Attach("common/fast_diamond_hash_m2.m");
Attach("C4x24/fast_full_aut_hash_m4_24_root.m");
X:=MDX11(4,24,GF(5));C:=Curve(X);assert Genus(X) eq 33;
cs:=Cusps(X);assert [#NoncuspidalPlaces(X,d):d in [1,2]] eq [0,0];
exports:=[* *];remaining:=[* *];
for dd in [4] do
 xs:=NoncuspidalPlaces(X,dd);assert #xs eq 96;
 rs,labels,keys:=PlacesUpToFullAutM4TwentyFourRoot(X,xs);
 print "M4N24_ROOT_D4_FULLAUT_SETUP",dd,#xs,#rs,Cputime(t0);
 x:=rs[1];E:=EllipticCurve(X,x);Ls:=LevelStructure(X,x);P:=Ls`P;Q:=Ls`Q;orb:={x};
 for par in M4TwentyFourRootAutParameters() do
  a,b,c,d:=Explode(par);La:=rec<recformat<P:PtEll,Q:PtEll>|P:=a*P+6*b*Q,Q:=c*P+d*Q>;
  y:=ModuliPointFast(X,E,La);assert Degree(y) eq dd;Include(~orb,y);
 end for;
 assert orb eq {xs[j]:j in [1..#xs]|labels[j] eq 1};
 print "M4N24_ROOT_D4_FULLAUT_REFERENCE_EQUAL",dd,#orb,Cputime(t0);
 tests:=[];
 for i in [1..#rs] do
  x:=rs[i];killed:=false;
  for q in [13,17,29] do
   A:=HeckeOperatorFast(X,q,x)-Divisor(x)-q*Divisor(DiamondOperator(X,q,x));assert Degree(A) eq 0;
   pr:=IsPrincipal(A);Append(~tests,[i,q,pr select 1 else 0]);
   print "M4N24_ROOT_D4_SPLIT_HECKE_CASE",dd,i,q,pr,Cputime(t0);
   if not pr then killed:=true;break;end if;
  end for;
  if not killed then Append(~remaining,[dd,i]);end if;
 end for;
 Append(~exports,[*dd,keys,labels,[Index(xs,x):x in rs],tests*]);
end for;
Write("C4x24/data/sieve_4_24_D4_fullaut_root_arithmetic.txt",Sprint([*exports,remaining*]):Overwrite:=true);
print "M4N24_ROOT_D4_FULLAUT_HECKE_DONE",remaining,Cputime(t0);quit;
