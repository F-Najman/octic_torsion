SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
X:=MDX11(4,24,GF(5));C:=Curve(X);assert Genus(X) eq 33;
cs:=Cusps(X);assert {*Degree(c):c in cs*} eq {*1^^16,2^^24*};
ns:=[#NoncuspidalPlaces(X,d):d in [1..4]];assert ns eq [0,0,0,96];
print "M4N24_ROOT_CONTROL_COUNTS",ns,Cputime(t0);
J:=JH(96,[25]);assert Dimension(J) eq 33;f:=FrobeniusPolynomial(J,5);g:=33;tr:=[];cl:=[];
for n in [1..4] do
 a:=-n*Coefficient(f,2*g-n);for j in [1..n-1] do a-:=Coefficient(f,2*g-j)*tr[n-j];end for;
 Append(~tr,Integers()!a);num:=5^n+1-tr[n];for e in Divisors(n) do if e lt n then num-:=e*cl[e];end if;end for;
 assert num mod n eq 0;Append(~cl,num div n);
end for;
assert cl eq [16,24,0,96];
print "M4N24_ROOT_CONTROL_NEWTON_COUNTS",cl,Cputime(t0);
Write("C4x24/data/control_4_24_D4_root_Frobenius.txt",Sprint([*Coefficients(f),cl*]):Overwrite:=true);
x:=NoncuspidalPlaces(X,4)[1];E:=EllipticCurve(X,x);Ls:=LevelStructure(X,x);
q:=13;fast:=HeckeOperatorFast(X,q,x);ys:=[];ms:=[];
for phi in MDIsogenies(E,q) do
 La:=rec<recformat<P:PtEll,Q:PtEll>|P:=phi(Ls`P),Q:=phi(Ls`Q)>;
 y:=ModuliPoint4(X,Codomain(phi),La);ii:=Index(ys,y);
 if ii eq 0 then Append(~ys,y);Append(~ms,1);else ms[ii]+:=1;end if;
end for;
ref:=DivisorGroup(C)!0;
for i in [1..#ys] do
 mu:=ms[i]*4/Degree(ys[i]);assert Denominator(mu) eq 1;
 ref+:=Integers()!mu*Divisor(ys[i]);
end for;
assert ref eq fast;print "M4N24_ROOT_CONTROL_HECKE_REFERENCE_EQUAL",q,Degree(ref),Cputime(t0);
print "M4N24_ROOT_CONTROLS_COMPLETE",Cputime(t0);quit;
