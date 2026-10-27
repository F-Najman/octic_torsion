/* Exact necessary class condition on the genus10 two-isogeny maps. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
Attach("common/fast_hecke_cyclic_commonzero.m");
Attach("common/class_rational_reduction.m");
Y:=MDX1(28,GF(3));C:=Curve(Y);assert Genus(Y) eq 10;
osraw:=CuspOrbitsQ(Y);os:=[Type(o) eq PlcCrvElt select Divisor(o) else o:o in osraw];
c0:=Support([o:o in os|Degree(o) eq 1][1])[1];
print "N28_D8_ISOGENY_QUOTIENT_CUSP_SETUP",[Degree(o):o in os],Cputime(t0);
G,m1,m2:=ClassGroupRationalReduction(C,c0:Proof:=true);
fp:=FrobeniusPolynomial(JOne(28),3);
assert Invariants(G) eq [4,4,24,936,0] and &*[n:n in Invariants(G)|n ne 0] eq Integers()!Evaluate(fp,1);
hg:=[m2(o-Degree(o)*Divisor(c0)):o in os];H:=sub<G|hg>;
assert Invariants(H) eq [2,4,12,936];
print "N28_D8_ISOGENY_QUOTIENT_GROUP",Invariants(G),Invariants(H),Cputime(t0);
X:=MDX11(2,28,GF(3));a:=eval Read("C2x28/data/N28_D8_fullaut_Hecke_input.txt");records:=a[1];
R:=PolynomialRing(GF(3));tests:=[* *];
for n in [35] do
 d,i,def,vals:=Explode(records[n]);assert d eq 8;K<z>:=ext<GF(3)|R!def>;
 bc:=[Evaluate(R!v,z):v in vals];x:=MDPlaceCommonZeros(X`_coordinates,bc);assert Degree(x) eq d;
 E:=EllipticCurve(X,x);ls:=LevelStructure(X,x);
 for j in [0,1] do
  T:=ls`P+j*14*ls`Q;assert T ne 14*ls`Q and 2*T eq E!0;
  S:=PolynomialRing(BaseRing(E));EG,phi:=IsogenyFromKernel(E,S![-T[1],1]);
  P:=phi(ls`Q);assert 28*P eq EG!0 and 14*P ne EG!0 and 4*P ne EG!0;
  level:=rec<recformat<P:PtEll>|P:=P>;y:=ModuliPointCyclicCommonZeros(Y,EG,level);
  ref:=ModuliPoint(Y,EG,level);assert y eq ref;
  e:=Degree(y);assert d mod e eq 0;mult:=d div e;
  cl:=mult*m2(Divisor(y)-e*Divisor(c0));inside:=cl in H;
  Append(~tests,[*n,j,d,e,mult,Eltseq(cl),inside*]);
  print "N28_D8_ISOGENY_QUOTIENT_TEST",n,j,d,e,mult,inside,Cputime(t0);
 end for;
end for;
Write("C2x28/data/N28_D8_isogeny_class_arithmetic.txt",Sprint([*Invariants(G),[Eltseq(r):r in Rows(RelationMatrix(G))],[Eltseq(h):h in hg],[Degree(o):o in os],tests,Coefficients(fp)*]):Overwrite:=true);
print "N28_D8_ISOGENY_QUOTIENT_CLASS_DONE",#tests,Cputime(t0);quit;
