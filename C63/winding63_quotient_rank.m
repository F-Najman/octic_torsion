SetColumns(0);SetQuitOnError(true);SetMemoryLimit(2*10^9);
N:=63;D:=FullDirichletGroup(N);total:=0;orbits:=0;
for chi0 in GaloisConjugacyRepresentatives(D) do
 if IsOdd(chi0) or Evaluate(chi0,8) ne 1 then continue;end if;
 chi:=MinimalBaseRingCharacter(chi0);S:=CuspidalSubspace(ModularSymbols(chi,2,1));
 dd:=Dimension(S);ww:=Dimension(WindingSubmodule(S));assert dd eq ww;
 total+:=EulerPhi(Order(chi))*dd;orbits+:=1;
 printf "Y63_RANK0_CHARACTER order=%o conductor=%o invariant8=%o dimension=%o winding_dimension=%o\n",Order(chi),Conductor(chi),Evaluate(chi,8),dd,ww;
end for;
assert total eq 49;
printf "Y63_RANK0_DONE invariant_character_orbits=%o dimension_sum=%o all_winding_full=true\n",orbits,total;
quit;
