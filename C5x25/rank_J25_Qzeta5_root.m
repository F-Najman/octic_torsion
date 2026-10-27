/* Exact modular-symbol support for the restriction-of-scalars rank lemma. */
SetColumns(0);SetQuitOnError(true);tt:=Cputime();
for N in [25,125] do
 D:=FullDirichletGroup(N);orders:=[];dims:=[];total:=0;
 for chi in GaloisConjugacyRepresentatives(D) do
  if IsOdd(chi) or Evaluate(chi,26) ne 1 then continue;end if;
  S:=CuspidalSubspace(ModularSymbols(MinimalBaseRingCharacter(chi),2,1));
  W:=WindingSubmodule(S);o:=Order(chi);d:=Dimension(S);w:=Dimension(W);
  if N eq 25 then
   assert o in [1,2,5,10];expected:=o eq 5 select 1 else (o eq 10 select 2 else 0);assert d eq expected;
  elif o in [5,10] then assert w eq d and d eq (o eq 5 select 7 else 8);
  end if;
  if o in [5,10] then assert Evaluate(chi,6) ne 1 and Conductor(chi) eq 25;end if;
  Append(~orders,o);Append(~dims,d);total+:=EulerPhi(o)*d;
  printf "J25_CYCLO5_RANK_SUPPORT level=%o order=%o conductor=%o dimension=%o winding=%o coefficient_degree=%o eigen6_nontrivial=%o\n",N,o,Conductor(chi),d,w,Degree(BaseRing(S)),Evaluate(chi,6) ne 1;
 end for;
 assert Sort(orders) eq [1,2,5,10] and total eq (N eq 25 select 12 else 76);
end for;
J:=JH(25,[1]);f:=FrobeniusPolynomial(J,11);assert Evaluate(f,1) eq 5014032203125 and IsOdd(Integers()!Evaluate(f,1));
A:=CompanionMatrix(FrobeniusPolynomial(J,2));N16:=17-Trace(A^4);
assert N16 eq 18;
printf "J25_CYCLO5_SPLIT11_ODD_ORDER %o\n",Evaluate(f,1);
printf "J25_CYCLO5_FINITE16_TOTAL_POINTS %o\n",N16;
print "J25_CYCLO5_RANK_SUPPORT_COMPLETE",Cputime(tt);quit;
