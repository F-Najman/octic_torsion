SetColumns(0);SetQuitOnError(true);SetMemoryLimit(2*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
load "common/class_group_unflattened.m";
X:=MDX1(21,GF(5));C:=Curve(X);c0:=[c:c in Cusps(X)|Degree(c) eq 1][1];
F:=AlgorithmicFunctionField(FunctionField(C));A0:=FunctionFieldDivisor(Divisor(c0));
assert RationalExtensionRepresentation(F) cmpeq RationalExtensionRepresentation(F);
G,gm,mp:=ClassGroup(F:ReductionDivisor:=A0,Proof:=true);
assert assigned RationalExtensionRepresentation(F)`cld`proven;
G2,gm2,mp2:=MDClassGroupUnflattened(F);assert G cmpeq G2;
assert RelationMatrix(G) eq RelationMatrix(G2);
tested:=0;
for d in [1..4] do
 ps:=Places(F,d);
 for i in [1..Minimum(10,#ps)] do
  dd:=Divisor(ps[i])-d*A0;assert mp(dd) eq mp2(dd);tested+:=1;
 end for;
end for;
for g in Generators(G) do
 assert mp(gm2(g)) eq g and mp2(gm(g)) eq g;
end for;
printf "UNFLATTENED_MAP_NATIVE_COMPARISON_DONE curve=X1_21 p=5 genus=%o invariants=%o place_tests=%o inverse_generators=%o time=%o\n",Genus(F),Invariants(G),tested,Ngens(G),Cputime(tt);
quit;
