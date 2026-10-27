/* Bind the abstract AF generator to the actual target modular functions. */
function C48SerializeBinding(a)
 return [[[Integers()!z:z in Coefficients(Numerator(v))],
          [Integers()!z:z in Coefficients(Denominator(v))]]:v in Eltseq(a)];
end function;
ff48,mp48:=AlgorithmicFunctionField(FunctionField(Curve(Y)));
assert RationalExtensionRepresentation(ff48) cmpeq FR;
binding48:=Open("C48/data/clean_run_v2/function_field_bindings.py","w");
fprintf binding48,"target_modular_coordinates = %o\n",[C48SerializeBinding(FR!mp48(z)):z in Y`_coordinates];
fprintf binding48,"base_cusp_signature = %o\n",[Valuation(z,c0):z in fs];
fprintf binding48,"base_cusp_degree = %o\n",Degree(c0);
delete binding48;print "C48_TARGET_FUNCTION_FIELD_BINDINGS_EXPORTED";
