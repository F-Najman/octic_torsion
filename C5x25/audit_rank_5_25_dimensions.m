SetColumns(0);SetQuitOnError(true);SetMemoryLimit(2*10^9);tt:=Cputime();weighted25:=0;weighted125:=0;
for N in [25,125] do
 D:=FullDirichletGroup(N);
 for chi in GaloisConjugacyRepresentatives(D) do
  if IsOdd(chi) then continue;end if;
  if N eq 125 and (Evaluate(chi,26) ne 1 or Order(chi) notin [5,10]) then continue;end if;
  S:=CuspidalSubspace(ModularSymbols(MinimalBaseRingCharacter(chi),2,1));d:=Dimension(S);w:=Dimension(WindingSubmodule(S));
  printf "AUDIT_J25_CYCLOTOMIC_WINDING level=%o character_order=%o conductor=%o dimension=%o winding=%o coefficient_degree=%o\n",N,Order(chi),Conductor(chi),d,w,Degree(BaseRing(S));
  if N eq 25 then
   if d ne 0 then assert Order(chi) in [5,10] and Conductor(chi) eq 25;assert d eq (Order(chi) eq 5 select 1 else 2);end if;
   weighted25+:=d*EulerPhi(Order(chi));
  else
   assert Conductor(chi) eq 25 and d eq w;assert d eq (Order(chi) eq 5 select 7 else 8);
   weighted125+:=d*EulerPhi(Order(chi));
  end if;
 end for;
end for;
assert weighted25 eq 12 and weighted125 eq 60;
print "AUDIT_J25_CYCLOTOMIC_RANK_INPUTS_PASS",weighted25,weighted125,Cputime(tt);quit;
