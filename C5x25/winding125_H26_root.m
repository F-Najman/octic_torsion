/* Characteristic-zero diagnostic only; no degree8 exclusion is asserted. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);t0:=Cputime();
N:=125;D:=FullDirichletGroup(N);total:=0;bad:=[];
for chi in GaloisConjugacyRepresentatives(D) do
 if IsOdd(chi) or Evaluate(chi,26) ne 1 then continue;end if;
 S:=CuspidalSubspace(ModularSymbols(MinimalBaseRingCharacter(chi),2,1));W:=WindingSubmodule(S);
 dim:=Dimension(S);wind:=Dimension(W);total+:=dim*EulerPhi(Order(chi));
 print "M5N25_H26_WINDING",Order(chi),Conductor(chi),dim,wind,Degree(BaseRing(S)),Cputime(t0);
 if dim ne wind then
  Append(~bad,Order(chi));
  print "M5N25_H26_BAD_COMPONENT",Order(chi),CharacteristicPolynomial(HeckeOperator(S,2));
 end if;
end for;
J:=JH(125,[26]);assert Dimension(J) eq total;
print "M5N25_H26_WINDING_DIAGNOSTIC_COMPLETE",total,bad,Cputime(t0);quit;
