/* Exact geometric snapshot for independent replay of the direct-RR witnesses.
   All prime ideals are in the native reversed factor-base order. */
Fgeom40:=RationalExtensionRepresentation(AlgorithmicFunctionField(FunctionField(C)));
assert Fgeom40`cld`proven and Fgeom40`cld`Cl cmpeq G;
Bgeom40:=BaseField(Fgeom40);kgeom40:=ConstantField(Fgeom40);
assert #kgeom40 eq 3;
function RatData40(a)
 return [[Integers()!c:c in Coefficients(Numerator(a))],
         [Integers()!c:c in Coefficients(Denominator(a))]];
end function;
function FFData40(a)
 data:=[RatData40(z):z in Eltseq(a)];
 P:=PolynomialRing(kgeom40);
 assert Fgeom40![Bgeom40!(P!v[1])/Bgeom40!(P!v[2]):v in data] eq a;
 return data;
end function;
fbgeom40:=Reverse([v`place:v in Fgeom40`cld`fb]);
OFgeom40:=MaximalOrderFinite(Fgeom40);
geomout40:=Open("C40/data/factorbase_geometry.py","w");
fprintf geomout40,"characteristic = 3\nfield_equation = %o\n",[RatData40(z):z in Coefficients(DefiningPolynomial(Fgeom40))];
fprintf geomout40,"finite_flags = %o\n",[Order(Ideal(P)) eq OFgeom40 select 1 else 0:P in fbgeom40];
fprintf geomout40,"factorbase_degrees = %o\n",[Degree(P):P in fbgeom40];
fprintf geomout40,"factorbase_ideal_bases = [\n";
for ii40 in [1..#fbgeom40] do
 ideal40:=Ideal(fbgeom40[ii40]);
 fprintf geomout40,"%o,\n",[FFData40(Fgeom40!z):z in Basis(ideal40)];
end for;
fprintf geomout40,"]\nfactorbase_classes = %o\n",[Eltseq(Fgeom40`cld`ZnToCl(Domain(Fgeom40`cld`ZnToCl).i)):i in [1..#fbgeom40]];
fprintf geomout40,"relations = %o\n",[Eltseq(r):r in Rows(RelationMatrix(G))];
delete geomout40;
print "C40_EXACT_FACTORBASE_GEOMETRY_EXPORTED",#fbgeom40;
