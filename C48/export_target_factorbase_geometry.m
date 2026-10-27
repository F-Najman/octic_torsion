/* Exact geometric snapshot for independent replay of the direct-RR witnesses.
   All prime ideals are in the native reversed factor-base order. */
Fgeom48:=RationalExtensionRepresentation(AlgorithmicFunctionField(FunctionField(C)));
assert Fgeom48`cld`proven and Fgeom48`cld`Cl cmpeq G;
Bgeom48:=BaseField(Fgeom48);kgeom48:=ConstantField(Fgeom48);
assert #kgeom48 eq 5;
function RatData48(a)
 return [[Integers()!c:c in Coefficients(Numerator(a))],
         [Integers()!c:c in Coefficients(Denominator(a))]];
end function;
function FFData48(a)
 data:=[RatData48(z):z in Eltseq(a)];
 P:=PolynomialRing(kgeom48);
 assert Fgeom48![Bgeom48!(P!v[1])/Bgeom48!(P!v[2]):v in data] eq a;
 return data;
end function;
fbgeom48:=Reverse([v`place:v in Fgeom48`cld`fb]);
OFgeom48:=MaximalOrderFinite(Fgeom48);
geomout48:=Open("C48/data/clean_run/factorbase_geometry.py","w");
fprintf geomout48,"characteristic = 5\nfield_equation = %o\n",[RatData48(z):z in Coefficients(DefiningPolynomial(Fgeom48))];
fprintf geomout48,"finite_flags = %o\n",[Order(Ideal(P)) eq OFgeom48 select 1 else 0:P in fbgeom48];
fprintf geomout48,"factorbase_degrees = %o\n",[Degree(P):P in fbgeom48];
fprintf geomout48,"factorbase_ideal_bases = [\n";
for ii48 in [1..#fbgeom48] do
 ideal48:=Ideal(fbgeom48[ii48]);
 fprintf geomout48,"%o,\n",[FFData48(Fgeom48!z):z in Basis(ideal48)];
end for;
fprintf geomout48,"]\nfactorbase_classes = %o\n",[Eltseq(Fgeom48`cld`ZnToCl(Domain(Fgeom48`cld`ZnToCl).i)):i in [1..#fbgeom48]];
fprintf geomout48,"relations = %o\n",[Eltseq(r):r in Rows(RelationMatrix(G))];
delete geomout48;
print "C48_EXACT_FACTORBASE_GEOMETRY_EXPORTED",#fbgeom48;
