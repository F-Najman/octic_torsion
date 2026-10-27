SetColumns(0); SetQuitOnError(true); SetMemoryLimit(10^9);
D:=FullDirichletGroup(175); H:=[1,36,71,106,141];
assert Seqset(H) eq Seqset([Integers()!((Integers(175)!36)^k):k in [0..4]]);
assert (Integers(175)!24)^3 eq -1;
assert not 24 in H cat [175-h:h in H];
assert not Integers()!((Integers(175)!24)^2) in H cat [175-h:h in H];
cnt:=0; moved:=0; orders:=[]; movedorders:=[];
for chi in GaloisConjugacyRepresentatives(D) do
 if IsOdd(chi) then continue; end if;
 if Evaluate(chi,36) ne 1 then continue; end if;
 assert 35 mod Conductor(chi) eq 0;
 cnt+:=EulerPhi(Order(chi)); Append(~orders,Order(chi));
 nontrivial:=Evaluate(chi,24) ne 1;
 if nontrivial then
  moved+:=EulerPhi(Order(chi)); Append(~movedorders,Order(chi));
  assert Order(chi) in [3,6,12];
 end if;
 printf "H_FIXED_CHARACTER order=%o conductor=%o delta24_nontrivial=%o\n",
  Order(chi),Conductor(chi),nontrivial;
end for;
assert Seqset(orders) eq {1,2,3,4,6,12};
assert Seqset(movedorders) eq {3,6,12};
assert cnt eq 12 and moved eq 8;
for psi in GaloisConjugacyRepresentatives(FullDirichletGroup(5)) do
 assert Evaluate(psi,24)^2 eq 1;
 printf "TWIST_CHARACTER order=%o conductor=%o psi24_squared_one=true\n",
  Order(psi),Conductor(psi);
end for;
print "FIVE_THIRTY_FIVE_CHARACTER_CHECK_DONE";
quit;
