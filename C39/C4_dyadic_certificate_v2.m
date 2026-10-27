/* Fresh genus9 quotient and 2-primary specialization certificate for D8. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt:=Cputime();
CQ:=SmallModularCurve(39);jraw:=jInvariant(CQ,39);CC:=ChangeRing(CQ,GF(2));
RC:=CoordinateRing(Ambient(CC));jn:=RC!Numerator(jraw);jd:=RC!Denominator(jraw);
AF,am:=AlgorithmicFunctionField(FunctionField(CC));jaf:=am(FunctionField(CC)!(jn/jd));
assert Genus(AF) eq 3;
cs:=Poles(jaf);widths:=[-Valuation(jaf,P):P in cs];
assert #cs eq 4 and &and[Degree(P) eq 1:P in cs] and Seqset(widths) eq {1,3,13,39};
printf "C39_C4_PREFLIGHT_BASE genus=3 degree=%o cusp_widths=%o jdegree=%o time=%o\n",Degree(AF),widths,&+widths,Cputime(tt);
fJ:=FrobeniusPolynomial(JZero(39),2);if Type(fJ) eq SeqEnum then fJ:=&*fJ;end if;
nj:=Integers()!Evaluate(fJ,1);
G,gm:=RayClassGroup(DivisorGroup(AF)!0);inv:=Invariants(G);
assert &*[a:a in inv|a ne 0] eq nj and #[a:a in inv|a eq 0] eq 1;
cols:=[Eltseq(Divisor(P)@@gm):P in cs];
expect:=[w eq 39 select 0 else (w eq 13 select 2 else (w eq 3 select 3 else 1)):w in widths];
chars:=[];
for r in CartesianPower({0..3},Ngens(G)) do
 rr:=[r[i]:i in [1..Ngens(G)]];
 if not &and[inv[i]*rr[i] mod 4 eq 0:i in [1..#inv]] then continue;end if;
 if &and[(&+[cols[i][j]*rr[j]:j in [1..#inv]]) mod 4 eq expect[i]:i in [1..4]] then Append(~chars,rr);end if;
end for;
printf "C39_C4_PREFLIGHT_CHARACTER Picard=%o Jorder=%o characters=%o cusp_classes=%o time=%o\n",inv,nj,chars,cols,Cputime(tt);
assert #chars eq 1;K:=BaseField(AF);Dzero:=DivisorGroup(AF)!0;
C4:=AbelianGroup([4]);hh:=hom<G->C4|[C4![a]:a in chars[1]]>;
assert #Image(hh) eq 4;
fields,auts,wits:=InternalRayClassField(Dzero,Kernel(hh):WithAut:=true);
assert #fields eq #auts and #fields eq #wits and #fields eq 1;
wt:=wits[1][1];wc:=Eltseq(wt);WT:=wits[1][2];CY:=WT;E0:=BaseField(CY);
assert #wc eq 2 and Degree(CY) eq 2 and Degree(E0) eq 2 and BaseField(E0) cmpeq AF;
assert BaseRing(Parent(wt)) cmpeq AF;
Z0<z0>:=PolynomialRing(AF);Z1<z1>:=PolynomialRing(E0);
assert DefiningPolynomial(E0) eq z0^2+z0+wc[1];
assert DefiningPolynomial(CY) eq z1^2+z1+(E0!wc[2])+E0.1*(E0!wc[1]);
sigma0:=hom<E0->E0|map<AF->E0|x:->E0!x>,E0.1+1>;
sigma:=hom<CY->CY|map<E0->CY|x:->CY!sigma0(x)>,CY.1+CY!E0.1>;
assert sigma(sigma(sigma(sigma(CY.1)))) eq CY.1;
assert sigma(sigma(CY.1)) eq CY.1+1;
assert &and[sigma(CY!v) eq CY!v:v in [AF.1,AF!K.1,AF!K!BaseField(K).1]];
printf "C39_C4_C4_CONSTRUCTED degree4=true cyclic4_verified=true time=%o\n",Cputime(tt);

// A fixed-generator simple model of the independently verified Witt tower.
e0:=wc[1];e1:=wc[2];assert e0 ne 0;
Q4<z>:=PolynomialRing(AF);
f4:=z^4+(1+e0)*z^2+e0*z+e1^2+e0*e1+e0^3;
assert IsIrreducible(f4);CY<zz>:=FunctionField(f4);
bb:=(zz^2+zz+CY!e1)/(CY!e0);
emb0:=hom<E0->CY|map<AF->CY|x:->CY!x>,bb>;
toSimple:=hom<WT->CY|emb0,zz>;
toWitt:=hom<CY->WT|map<AF->WT|x:->WT!x>,WT.1>;
assert toSimple(toWitt(zz)) eq zz;
assert toWitt(toSimple(WT.1)) eq WT.1;
assert toWitt(toSimple(WT!E0.1)) eq WT!E0.1;
sigma:=hom<CY->CY|map<AF->CY|x:->CY!x>,zz+bb>;
assert sigma(sigma(zz)) eq zz+1 and sigma(sigma(sigma(sigma(zz)))) eq zz;
printf "C39_C4_SIMPLE_WITT_MODEL exact_inverse_maps=true cyclic4=true time=%o\n",Cputime(tt);
// Independent length-two Witt trace: addition is (a,b)+(c,d)=(a+c,b+d+ac).
function WittTrace2(a,b)
 F:=Parent(a);assert Parent(b) cmpeq F;
 A:=F!0;B:=F!0;
 for i in [0..Degree(F,GF(2))-1] do
  c:=a^(2^i);d:=b^(2^i);B:=B+d+A*c;A:=A+c;
 end for;
 assert A^2 eq A and B^2 eq B;
 return ((Integers()!A)+2*(Integers()!B)) mod 4;
end function;

G4,q4:=quo<G|[4*g:g in Generators(G)]>;generated:=sub<G4|>;
units:={1,3};pe:=PlaceEnumInit(AF);tried:=0;tested:=0;
repeat
 P:=PlaceEnumNext(pe);tried+:=1;assert tried le 100;
 if &or[Valuation(x,P) lt 0:x in wc] then continue;end if;
 F:=ResidueClassField(P);a:=F!Evaluate(wc[1],P);b:=F!Evaluate(wc[2],P);
 actual:=WittTrace2(a,b);cl:=Divisor(P)@@gm;expected:=Eltseq(hh(cl))[1];
 units:={u:u in units|actual eq u*expected mod 4};assert #units gt 0;
 generated:=sub<G4|generated,q4(cl)>;tested+:=1;
 printf "C39_C4_WITT_FROBENIUS degree=%o class=%o actual=%o expected=%o generated_order=%o quotient_order=%o time=%o\n",Degree(P),Eltseq(cl),actual,expected,#generated,#G4,Cputime(tt);
until generated eq G4;
printf "C39_C4_WITT_CHARACTER_VERIFIED samples=%o orientation_units=%o all_Picard_mod4=true time=%o\n",tested,units,Cputime(tt);
function RelativeDecomposition(F,P)
 assert BaseField(F) cmpeq FunctionField(P);
 f,g:=TwoGenerators(P);unif:=UniformizingElement(P);
 assert Valuation(f,P) gt 0 and Valuation(g,P) gt 0 and Valuation(unif,P) eq 1;
 raw:=Decomposition(F,P);
 above:=[Q:Q in raw|Valuation(F!f,Q) gt 0 and Valuation(F!g,Q) gt 0];
 assert #above gt 0 and &+[Degree(Q):Q in above] eq Degree(F)*Degree(P);
 assert &and[Valuation(F!unif,Q) eq 1:Q in above];
 return above;
end function;
function AboveAF(F,P)
 out:=RelativeDecomposition(F,P);
 assert &+[Degree(Q):Q in out] eq 4*Degree(P);
 return out;
end function;

cys:=[];cywidths:=[];
for i in [1..4] do
 pp:=AboveAF(CY,cs[i]);
 want:=Order(hh(Divisor(cs[i])@@gm));
 assert #pp eq 4 div want and &and[Degree(P) eq want:P in pp];
 cys cat:=pp;cywidths cat:=[widths[i]:P in pp];
 printf "C39_C4_CUSP_FIBRE width=%o expected_frobenius_order=%o degrees=%o time=%o\n",widths[i],want,[Degree(P):P in pp],Cputime(tt);
end for;
assert {*Degree(P):P in cys*} eq {*1^^4,2^^2,4^^2*};
assert &+[cywidths[i]*Degree(cys[i]):i in [1..#cys]] eq 224;
assert Genus(CY) eq 9;
printf "C39_C4_COVER_GEOMETRY genus9=true split_rational_cusp=true unramified_by_RH=true time=%o\n",Cputime(tt);

GY,gym:=RayClassGroup(DivisorGroup(CY)!0);
assert &*[a:a in Invariants(GY)|a ne 0] eq 2016;
c0Y:=[Q:Q in cys|Degree(Q) eq 1][1];
gensY:=[(Divisor(Q)-Degree(Q)*Divisor(c0Y))@@gym:Q in cys];
HY:=sub<GY|gensY>;assert Valuation(#HY,2) eq 5;
printf "C39_C4_DYADIC_IMAGE Picard=%o cusp_image=%o two_primary_order=32 time=%o\n",Invariants(GY),Invariants(HY),Cputime(tt);
ar:=Open("C39/data/C4_cusp_arithmetic.py","w");
fprintf ar,"relations = %o\n",[Eltseq(r):r in Rows(RelationMatrix(GY))];
fprintf ar,"cusp_generators = %o\n",[Eltseq(g):g in gensY];
fprintf ar,"cusp_degrees = %o\n",[Degree(Q):Q in cys];
fprintf ar,"cusp_widths = %o\n",cywidths;
fprintf ar,"complete = True\n";delete ar;
print "C39_C4_DYADIC_CERTIFICATE_DONE complete=true";quit;
