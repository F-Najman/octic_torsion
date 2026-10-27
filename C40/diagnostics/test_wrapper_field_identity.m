SetColumns(0);SetQuitOnError(true);SetMemoryLimit(2*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/class_rational_reduction.m");
load "common/class_group_unflattened.m";
X:=MDX1(21,GF(5));C:=Curve(X);c0:=[c:c in Cusps(X)|Degree(c) eq 1][1];
Fbefore:=RationalExtensionRepresentation(AlgorithmicFunctionField(FunctionField(C)));
assert not assigned Fbefore`cld;
G,gm,mp:=ClassGroupRationalReduction(C,c0:Proof:=true);
Fafter:=RationalExtensionRepresentation(AlgorithmicFunctionField(FunctionField(C)));
assert Fbefore cmpeq Fafter;
assert assigned Fafter`cld and Fafter`cld`proven and assigned Fafter`cld`Cl;
assert Fafter`cld`Cl cmpeq G;
G2,gm2,mp2:=MDClassGroupUnflattened(AlgorithmicFunctionField(FunctionField(C)));
assert G2 cmpeq G and RelationMatrix(G2) eq RelationMatrix(G);
for d in [1..4] do
 ps:=Places(C,d);
 for i in [1..Minimum(10,#ps)] do
  D:=Divisor(ps[i])-d*Divisor(c0);assert mp(D) eq mp2(FunctionFieldDivisor(D));
 end for;
end for;
printf "WRAPPER_FIELD_CACHE_IDENTITY_DONE same_parent=true proven=true same_group=true same_relations=true all_place_maps=true time=%o\n",Cputime(tt);
quit;
