/* Fresh degree8 use of the degree-independent exact A3 quotient model.
   Generation histories and earlier-degree exclusions are not reused. */
if not assigned A3p then A3p:=7;end if;A3X15:=MDX1(15,GF(A3p));A3C15:=Curve(A3X15);
A3F15:=FunctionField(A3C15);
A3A15,A3ma:=AlgorithmicFunctionField(A3F15);
assert Genus(A3A15) eq 1;
A3E15:=EllipticCurve([A3ma(a):a in A3X15`_E]);
A3P15:=A3E15![A3ma(a):a in A3X15`_P];assert Order(A3P15) eq 15;
A3disc15:=Discriminant(A3E15);A3j15:=jInvariant(A3E15);
A3R<A3v>:=PolynomialRing(A3A15);A3f:=A3v^2-A3disc15;
assert IsIrreducible(A3f) and GCD(A3f,Derivative(A3f)) eq 1;
A3K<A3root>:=FunctionField(A3f);assert Degree(A3K) eq 2 and Genus(A3K) eq 9;
function A3PlacesAboveBase(pp)
 aa,bb:=TwoGenerators(pp);
 raw:=Setseq(Seqset(Decomposition(A3K,pp)));
 // In a tower, Decomposition may return places over other intermediate
 // places above the same rational base place. Select the actual base
 // place using both ideal generators before using residue degrees.
 return [qq:qq in raw|Valuation(A3K!aa,qq) gt 0 and Valuation(A3K!bb,qq) gt 0];
end function;
A3basecs:=Poles(A3j15);A3cs:=[];
for pp in A3basecs do
 assert IsOdd(Valuation(A3disc15,pp));
 qq:=A3PlacesAboveBase(pp);assert #qq eq 1 and Degree(qq[1]) eq Degree(pp);
 assert Valuation(A3K!A3j15,qq[1]) lt 0;Append(~A3cs,qq[1]);
end for;
assert #Seqset(A3cs) eq #A3cs and &+[Degree(q):q in A3cs] eq 16;
assert &+[-Valuation(A3K!A3j15,q)*Degree(q):q in A3cs] eq 192;
// Every base cusp is totally ramified, with unchanged residue field.
// Consequently each rational base cusp orbit has its unique rational
// orbit upstairs. This is exactly one half its conorm divisor.
A3raworbits:=CuspOrbitsQ(A3X15);A3orbitdict:=AssociativeArray();
for oo in A3raworbits do
 ss:=Support(oo);sg:=CuspSignature(A3X15,ss[1]);assert &and[CuspSignature(A3X15,c) eq sg:c in ss];
 assert not IsDefined(A3orbitdict,sg);A3orbitdict[sg]:=oo;
end for;
A3signatures:=Sort(Setseq(Keys(A3orbitdict)));A3baseorbits:=[A3orbitdict[sg]:sg in A3signatures];A3orbits:=[];
for OO in A3baseorbits do
 DD:=FunctionFieldDivisor(OO);assert FunctionField(DD) eq A3A15;
 ss,mm:=Support(DD);UU:=DivisorGroup(A3K)!0;
 for i in [1..#ss] do
  qq:=A3PlacesAboveBase(ss[i]);assert #qq eq 1 and Degree(qq[1]) eq Degree(ss[i]);
  UU+:=mm[i]*Divisor(qq[1]);
 end for;
 assert Degree(UU) eq Degree(OO);Append(~A3orbits,UU);
end for;
assert Sort([Degree(OO):OO in A3orbits]) eq [1,1,1,1,2,2,4,4];
A3baseindex:=Index([Degree(OO):OO in A3orbits],1);
A3c0:=Support(A3orbits[A3baseindex])[1];assert Degree(A3c0) eq 1;
A3coordinates:=[A3K!A3ma(A3F15.1),A3K!A3ma(A3F15.2),A3root];
