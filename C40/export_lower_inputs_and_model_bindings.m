/* Small additional exports from the retained proven parent. No class search. */
lowerout40:=Open("C40/data/fixed_lower_inputs.py","w");
fprintf lowerout40,"inputs = [\n";
for dd40 in [4..7] do
 KK40:=Parent(Evaluate(X`_coordinates[1],reps40[dd40][1]));assert Degree(KK40) eq dd40;
 ff40:=[Integers()!a:a in Coefficients(DefiningPolynomial(KK40,GF(3)))];
 coords40:=[];
 for PP40 in reps40[dd40] do
  bc40:=MDM2ToCommonField([Evaluate(f,PP40):f in X`_coordinates],KK40);
  Append(~coords40,[[Integers()!a:a in Eltseq(z,GF(3))]:z in bc40]);
 end for;
 fprintf lowerout40,"[%o,%o,%o],\n",dd40,ff40,coords40;
end for;
fprintf lowerout40,"]\n";delete lowerout40;
modelout40:=Open("C40/data/function_field_bindings.py","w");
AFbind40,mpbind40:=AlgorithmicFunctionField(FunctionField(C));
assert RationalExtensionRepresentation(AFbind40) cmpeq Fgeom40;
fprintf modelout40,"sutherland_coordinates = %o\nuniversal_tate_coordinates = %o\nbase_cusp_signature = %o\n",[FFData40(Fgeom40!mpbind40(f)):f in X`_coordinates_xy],[FFData40(Fgeom40!mpbind40(f)):f in X`_coordinates],CuspSignature(X,c0);
fprintf modelout40,"retained_parent_pid = %o\n",Getpid();delete modelout40;
print "C40_LOWER_INPUTS_AND_EXACT_MODEL_BINDINGS_EXPORTED";
