SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
X:=MDX11(3,24,GF(7));C:=Curve(X);a:=eval Read("C3x24/data/preflight_3_24_D4_Aut48_input.txt");
d,i,def,vals:=Explode(a[1][1]);assert d eq 3;R:=PolynomialRing(GF(7));K<z>:=ext<GF(7)|R!def>;
qt:=[Evaluate(R!v,z):v in vals];x:=MDPlaceCommonZeros(X`_coordinates,qt);assert Degree(x) eq 3;
T:=HeckeOperatorFast(X,13,x);E:=EllipticCurve(X,x);Ls:=LevelStructure(X,x);ys:=[];ms:=[];
for phi in MDIsogenies(E,13) do
 La:=rec<recformat<P:PtEll,Q:PtEll>|P:=phi(Ls`P),Q:=phi(Ls`Q)>;
 y:=ModuliPointFast(X,Codomain(phi),La);idx:=Index(ys,y);
 if idx eq 0 then Append(~ys,y);Append(~ms,1);else ms[idx]+:=1;end if;
end for;
ref:=DivisorGroup(C)!0;
for i in [1..#ys] do mu:=ms[i]*3/Degree(ys[i]);assert Denominator(mu) eq 1;ref+:=Integers()!mu*Divisor(ys[i]);end for;
assert T eq ref and Degree(ref) eq 42;
print "M3N24_D4_T13_GENERIC_REFERENCE_EQUAL",Degree(ref),Cputime(t0);print "M3N24_D4_T13_REFERENCE_DONE complete=true";quit;
