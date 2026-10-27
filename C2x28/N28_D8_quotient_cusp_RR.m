/* Four rational maps to X1(28); short GLOBAL Q-cusp classes and cross-prime RR. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
Attach("common/fast_hecke_cyclic_commonzero.m");
Attach("common/class_rational_reduction.m");
function OrderedCuspOrbits(Y)
 raw:=CuspOrbitsQ(Y);os:=[Type(o) eq PlcCrvElt select Divisor(o) else o:o in raw];
 ss:=[CuspSignature(Y,Support(o)[1]):o in os];assert #os eq 15 and #Seqset(ss) eq 15;
 for i in [1..#os] do assert &and[CuspSignature(Y,c) eq ss[i]:c in Support(os[i])];end for;
 order:=Sort([<ss[i],i>:i in [1..#os]]);return [os[x[2]]:x in order],[x[1]:x in order];
end function;
Y:=MDX1(28,GF(3));C:=Curve(Y);os,sigs:=OrderedCuspOrbits(Y);ds:=[Degree(o):o in os];base:=Index(ds,1);c0:=Support(os[base])[1];
G,m1,m2:=ClassGroupRationalReduction(C,c0:Proof:=true);
hg:=[m2(os[i]-ds[i]*Divisor(c0)):i in [1..15]];Free:=FreeAbelianGroup(15);h:=hom<Free->G|hg>;
assert Invariants(Image(h)) eq [2,4,12,936];
KB:=Matrix(Integers(),[Eltseq(Free!v):v in Generators(Kernel(h))]);LL:=LatticeWithBasis(LLL(KB));
X:=MDX11(2,28,GF(3));input:=eval Read("C2x28/data/N28_D8_fullaut_Hecke_input.txt");d,i,def,vals:=Explode(input[1][35]);assert d eq 8;
R:=PolynomialRing(GF(3));K<z>:=ext<GF(3)|R!def>;bc:=[Evaluate(R!v,z):v in vals];x:=MDPlaceCommonZeros(X`_coordinates,bc);
E:=EllipticCurve(X,x);ls:=LevelStructure(X,x);lifts:=[];classes:=[];local3:=[* *];
for kind in [0..3] do
 if kind lt 2 then EG:=E;P:=ls`Q+kind*ls`P;
 else T:=ls`P+(kind-2)*14*ls`Q;S:=PolynomialRing(BaseRing(E));EG,phi:=IsogenyFromKernel(E,S![-T[1],1]);P:=phi(ls`Q);end if;
 y:=ModuliPointCyclicCommonZeros(Y,EG,rec<recformat<P:PtEll>|P:=P>);assert 8 mod Degree(y) eq 0;mult:=8 div Degree(y);
 cl:=mult*m2(Divisor(y)-Degree(y)*Divisor(c0));assert cl in Image(h);
 v:=Eltseq(cl@@h);w:=ClosestVector(LL,Vector(Rationals(),v));short:=[v[i]-Integers()!w[i]:i in [1..15]];assert h(Free!short) eq cl;
 D:=8*Divisor(c0)+&+[short[i]*(os[i]-ds[i]*Divisor(c0)):i in [1..15]];assert IsPrincipal(D-mult*Divisor(y));
 Append(~lifts,short);Append(~classes,Eltseq(cl));Append(~local3,[*kind,Degree(y),mult,Dimension(RiemannRochSpace(D))*]);
 print "N28_D8_QUOTIENT_SHORT_LIFT",kind,short,"local3",local3[#local3],Cputime(t0);
end for;
outcomes:=[* *];
for p in [5,11,13] do
 YY:=MDX1(28,GF(p));oo,ss:=OrderedCuspOrbits(YY);assert ss eq sigs and [Degree(o):o in oo] eq ds;
 c00:=Support(oo[base])[1];
 for k in [1..#lifts] do
  D:=8*Divisor(c00)+&+[lifts[k][i]*(oo[i]-ds[i]*Divisor(c00)):i in [1..15]];
  V,mp:=RiemannRochSpace(D);dim:=Dimension(V);degs:=[];mults:=[];allcusps:=false;
  if dim eq 1 then T:=D+Divisor(mp(V.1));assert IsEffective(T) and Degree(T) eq 8;sp,mults:=Support(T);degs:=[Degree(s):s in sp];allcusps:=&and[IsCusp(YY,s):s in sp];end if;
  Append(~outcomes,[*p,k-1,dim,degs,mults,allcusps*]);print "N28_D8_QUOTIENT_CROSS_PRIME_RR",outcomes[#outcomes],Cputime(t0);
 end for;
end for;
Write("C2x28/data/N28_D8_quotient_cusp_RR_arithmetic.txt",Sprint([*sigs,ds,base,lifts,classes,local3,outcomes*]):Overwrite:=true);
print "N28_D8_QUOTIENT_CUSP_RR_DONE",#outcomes,Cputime(t0);quit;
