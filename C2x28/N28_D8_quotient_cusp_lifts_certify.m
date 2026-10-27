/* Fresh actual-relation certificate linking the FINAL source orbit to
   the signed GLOBAL cusp divisors used in characteristic-zero RR. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
Attach("common/fast_hecke_cyclic_commonzero.m");
Attach("common/class_rational_reduction.m");
a:=eval Read("C2x28/data/N28_D8_quotient_cusp_RR_arithmetic.txt");sigs:=a[1];ds:=a[2];base:=a[3];lifts:=a[4];
Y:=MDX1(28,GF(3));C:=Curve(Y);raw:=CuspOrbitsQ(Y);os:=[Type(o) eq PlcCrvElt select Divisor(o) else o:o in raw];
ss:=[CuspSignature(Y,Support(o)[1]):o in os];order:=Sort([<ss[i],i>:i in [1..#os]]);os:=[os[x[2]]:x in order];
assert [x[1]:x in order] eq sigs and [Degree(o):o in os] eq ds;c0:=Support(os[base])[1];
G,m1,m2:=ClassGroupRationalReduction(C,c0:Proof:=true);hg:=[m2(os[i]-ds[i]*Divisor(c0)):i in [1..15]];
Free:=FreeAbelianGroup(15);h:=hom<Free->G|hg>;assert Invariants(Image(h)) eq [2,4,12,936];
X:=MDX11(2,28,GF(3));input:=eval Read("C2x28/data/N28_D8_fullaut_Hecke_input.txt");d,i,def,vals:=Explode(input[1][35]);assert d eq 8;
R:=PolynomialRing(GF(3));K<z>:=ext<GF(3)|R!def>;bc:=[Evaluate(R!v,z):v in vals];x:=MDPlaceCommonZeros(X`_coordinates,bc);E:=EllipticCurve(X,x);ls:=LevelStructure(X,x);tests:=[* *];
for kind in [0..3] do
 if kind lt 2 then EG:=E;P:=ls`Q+kind*ls`P;
 else T:=ls`P+(kind-2)*14*ls`Q;S:=PolynomialRing(BaseRing(E));EG,phi:=IsogenyFromKernel(E,S![-T[1],1]);P:=phi(ls`Q);end if;
 level:=rec<recformat<P:PtEll>|P:=P>;y:=ModuliPointCyclicCommonZeros(Y,EG,level);assert y eq ModuliPoint(Y,EG,level);
 assert 8 mod Degree(y) eq 0;mult:=8 div Degree(y);cl:=mult*m2(Divisor(y)-Degree(y)*Divisor(c0));
 assert h(Free!lifts[kind+1]) eq cl;
 Append(~tests,[*kind,Degree(y),mult,Eltseq(cl)*]);
 print "N28_D8_EXACT_GLOBAL_CUSP_LIFT_MATCH",kind,Degree(y),mult,Cputime(t0);
end for;
Write("C2x28/data/N28_D8_quotient_cusp_lifts_certify_arithmetic.txt",Sprint([*sigs,ds,base,lifts,Invariants(G),[Eltseq(r):r in Rows(RelationMatrix(G))],[Eltseq(v):v in hg],tests*]):Overwrite:=true);
print "N28_D8_EXACT_GLOBAL_CUSP_LIFTS_CERTIFIED",Cputime(t0);quit;
