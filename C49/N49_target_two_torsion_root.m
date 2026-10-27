SetColumns(0);SetQuitOnError(true);
E:=SmallModularCurve(49);
assert Type(E) eq CrvEll and aInvariants(E) eq [1,-1,0,-2,-1];
assert Discriminant(E) eq -343 and jInvariant(E) eq -3375;
j49:=jInvariant(E,49);FF:=FunctionField(E);AF,am:=AlgorithmicFunctionField(FF);assert Degree(am(FF!j49)) eq 56;
T:=E![2,-1];assert T ne E!0 and 2*T eq E!0;
E2:=ChangeRing(E,GF(2));T2:=E2![0,1];assert T2 ne E2!0 and 2*T2 eq E2!0;
print "N49_X0_TARGET_TWO_TORSION_VERIFIED native_modular_model=true jdegree56=true discriminant_minus343=true rational_T_order2=true reduced_T_nonzero=true";quit;
