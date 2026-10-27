SetColumns(0);SetQuitOnError(true);SetMemoryLimit(3*10^9);tt75:=Cputime();
N:=75;DD:=FullDirichletGroup(N);total75:=0;bad75:=[];
for chi in GaloisConjugacyRepresentatives(DD) do
 if IsOdd(chi) then continue;end if;
 S:=CuspidalSubspace(ModularSymbols(MinimalBaseRingCharacter(chi),2,1));
 W:=WindingSubmodule(S);d:=Dimension(S);w:=Dimension(W);
 total75+:=d*EulerPhi(Order(chi));
 printf "OCTIC75_WINDING order=%o conductor=%o dimension=%o winding=%o coefficient_degree=%o time=%o\n",Order(chi),Conductor(chi),d,w,Degree(BaseRing(S)),Cputime(tt75);
 if d ne w then
  Append(~bad75,chi);
  if Order(chi) eq 1 then
   BS:=Matrix([Eltseq(v):v in Basis(S)]);BW:=Matrix([Eltseq(v):v in Basis(W)]);
   VV:=VectorSpace(Rationals(),d);bW:=[VV!Solution(BS,r):r in Rows(BW)];BB:=bW;
   for v in Basis(VV) do if v notin sub<VV|BB> then Append(~BB,v);end if;end for;
   B:=Matrix(BB);
   for q in [2,3,5,7,11] do
    T:=B*HeckeOperator(S,q)*B^-1;
    assert Submatrix(T,1,w+1,w,d-w) eq 0;
    printf "OCTIC75_UPPER_TRIVIAL_HECKE q=%o quotient=%o\n",q,Submatrix(T,w+1,w+1,d-w,d-w);
   end for;
  end if;
 end if;
end for;
assert total75 eq 145;
assert #bad75 eq 0;
print "C75_RANK_ZERO_WINDING_COMPLETE",total75,Cputime(tt75);quit;
