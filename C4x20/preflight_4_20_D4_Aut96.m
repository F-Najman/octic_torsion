SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
Attach("common/fast_diamond_hash_m2.m");
Attach("common/fast_full_aut_hash_m4_16.m");
Attach("common/fast_full_aut_hash_m4_16_relative.m");
Attach("C4x20/fast_full_aut_hash_m4_20.m");
X:=MDX11(4,20,GF(9));assert Genus(X) eq 25;coverage:=[* *];
for d in [2..4] do
 xs:=NoncuspidalPlaces(X,d);assert #xs eq [96,320,888][d-1];
 rs,labels,keys:=PlacesUpToFullAutM4Twenty(X,xs);print "M4N20_D4_AUT96_PARTITION",d,#xs,#rs,Cputime(t0);
 x:=rs[1];E:=EllipticCurve(X,x);Ls:=LevelStructure(X,x);P:=Ls`P;Q:=Ls`Q;Q4:=5*Q;Q5:=4*Q;orb:={x};
 for par in M4TwentyAutParameters() do
  a,b,c,dd,u:=Explode(par);La:=rec<recformat<P:PtEll,Q:PtEll>|P:=a*P+b*Q4,Q:=c*P+dd*Q4-u*Q5>;
  Include(~orb,ModuliPointFast(X,E,La));
 end for;
 assert orb eq {xs[j]:j in [1..#xs]|labels[j] eq 1};
 K:=ext<GF(9)|d>;orbs:=[M4TwentyOrbitKeys(X,x,K):x in rs];
 Append(~coverage,[*d,keys,labels,orbs,[Index(xs,x):x in rs]*]);
 print "M4N20_D4_AUT96_REFERENCE_EQUAL",d,#orb,Cputime(t0);
end for;
Write("C4x20/data/preflight_4_20_D4_Aut96_arithmetic.txt",Sprint(coverage):Overwrite:=true);
print "M4N20_D4_AUT96_DONE complete=true",Cputime(t0);quit;
