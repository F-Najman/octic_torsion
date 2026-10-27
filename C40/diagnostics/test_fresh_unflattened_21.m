SetColumns(0);SetQuitOnError(true);SetMemoryLimit(2*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
load "common/class_group_unflattened.m";
load "common/native_group_generators.m";
X:=MDX1(21,GF(5));C:=Curve(X);c0:=[c:c in Cusps(X)|Degree(c) eq 1][1];
F:=AlgorithmicFunctionField(FunctionField(C));FR:=RationalExtensionRepresentation(F);
A0:=FunctionFieldDivisor(Divisor(c0));
iv:=ClassGroupAbelianInvariants(F:ReductionDivisor:=A0,Proof:=true);
assert FR`cld`proven and not assigned FR`cld`DivToClFlat;
MDProvenClassGroupGenerators(FR);
G,gm,mp:=MDClassGroupUnflattened(F);
assert Invariants(G) eq iv;
printf "FRESH_UNFLATTENED_GROUP_READY invariants=%o no_native_flat_maps=true time=%o\n",iv,Cputime(tt);
Gnative,gmnative,mpnative:=ClassGroup(F:Proof:=true);
assert RelationMatrix(G) eq RelationMatrix(Gnative);
tested:=0;
for d in [1..4] do
 ps:=Places(F,d);
 for i in [1..Minimum(10,#ps)] do
  D:=Divisor(ps[i])-d*A0;
  assert Gnative!Eltseq(mp(D)) eq mpnative(D);tested+:=1;
 end for;
end for;
for g in Generators(G) do assert mpnative(gm(g)) eq Gnative!Eltseq(g);end for;
printf "FRESH_UNFLATTENED_NATIVE_COMPARISON_DONE proven_before_maps=true same_relations=true place_tests=%o exact_inverse_generators=true time=%o\n",tested,Cputime(tt);
quit;
