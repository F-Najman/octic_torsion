SetColumns(0);SetQuitOnError(true);SetMemoryLimit(16*10^9);tt:=Cputime();
J:=JZero(55);fj:=FrobeniusPolynomial(J,2);if Type(fj) eq SeqEnum then fj:=&*fj;end if;
nj:=Integers()!Evaluate(fj,1);
printf "OCTIC55_X0_JACOBIAN order=%o factorization=%o\n",nj,Factorization(nj);
CQ:=SmallModularCurve(55);jraw:=jInvariant(CQ,55);
CC:=ChangeRing(CQ,GF(2));RC:=CoordinateRing(Ambient(CC));
jn:=RC!Numerator(jraw);jd:=RC!Denominator(jraw);
AF,am:=AlgorithmicFunctionField(FunctionField(CC));K:=BaseField(AF);
jaf:=am(FunctionField(CC)!(jn/jd));assert Degree(AF) eq 4 and Genus(AF) eq 5;
for ell in [5,11] do
 ww:=[RC!e:e in DefiningPolynomials(AtkinLehnerInvolution(CQ,55,ell))];
 je:=am(FunctionField(CC)!(Evaluate(jn,ww)/Evaluate(jd,ww)));
 assert je ne jaf and Evaluate(ClassicalModularPolynomial(ell),[jaf,je]) eq 0;
end for;
cs:=Poles(jaf);widths:=[-Valuation(jaf,c):c in cs];
assert #cs eq 4 and &and[Degree(c) eq 1:c in cs] and Seqset(widths) eq {1,5,11,55};
assert &+[Degree(c)*(-Valuation(jaf,c)):c in cs] eq 72;
printf "OCTIC55_X0_MODEL genus5=true rational_base_degree4=true modular_polynomials5and11_verified=true cusp_widths=%o time=%o\n",widths,Cputime(tt);
DG:=DivisorGroup(AF);Dzero:=DG!0;G,gm:=RayClassGroup(Dzero);
inv:=Invariants(G);assert &*[a:a in inv|a ne 0] eq nj and #[a:a in inv|a eq 0] eq 1;
printf "OCTIC55_PICARD invariants=%o time=%o\n",inv,Cputime(tt);
cols:=[];
for i in [1..4] do
 Append(~cols,Eltseq(Divisor(cs[i])@@gm));
 printf "OCTIC55_CUSP_CLASS index=%o class=%o time=%o\n",i,cols[i],Cputime(tt);
end for;
expect:=[w eq 55 select 0 else (w eq 11 select 1 else (w eq 5 select 2 else 3)):w in widths];
chars:=[];
for r in CartesianPower({0..3},Ngens(G)) do
 rr:=[r[i]:i in [1..Ngens(G)]];
 if not &and[inv[i]*rr[i] mod 4 eq 0:i in [1..#inv]] then continue;end if;
 if &and[(&+[cols[i][j]*rr[j]:j in [1..#inv]]) mod 4 eq expect[i]:i in [1..4]] then Append(~chars,rr);end if;
end for;
printf "OCTIC55_C4_CUSP_CHARACTERS count=%o images=%o time=%o\n",#chars,chars,Cputime(tt);
assert #chars eq 1;
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
printf "OCTIC55_C4_CONSTRUCTED degree4=true cyclic4_verified=true time=%o\n",Cputime(tt);

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
printf "OCTIC55_SIMPLE_WITT_MODEL exact_inverse_maps=true cyclic4=true time=%o\n",Cputime(tt);
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
 printf "OCTIC55_WITT_FROBENIUS degree=%o class=%o actual=%o expected=%o generated_order=%o quotient_order=%o time=%o\n",Degree(P),Eltseq(cl),actual,expected,#generated,#G4,Cputime(tt);
until generated eq G4;
printf "OCTIC55_WITT_CHARACTER_VERIFIED samples=%o orientation_units=%o all_Picard_mod4=true time=%o\n",tested,units,Cputime(tt);
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
 printf "OCTIC55_CUSP_FIBRE width=%o expected_frobenius_order=%o degrees=%o time=%o\n",widths[i],want,[Degree(P):P in pp],Cputime(tt);
end for;
assert {*Degree(P):P in cys*} eq {*1^^4,2^^2,4^^2*};
assert &+[cywidths[i]*Degree(cys[i]):i in [1..#cys]] eq 288;
assert Genus(CY) eq 17;
printf "OCTIC55_COVER_GEOMETRY genus17=true split_rational_cusp=true unramified_by_RH=true time=%o\n",Cputime(tt);

print "OCTIC55_VERIFIED_GENUS17_MODEL_DONE";
