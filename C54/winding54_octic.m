SetColumns(0);SetQuitOnError(true);SetMemoryLimit(3*10^9);tt54:=Cputime();
N:=54;DD:=FullDirichletGroup(N);total54:=0;bad54:=[];
for chi in GaloisConjugacyRepresentatives(DD) do
 if IsOdd(chi) then continue;end if;
 S:=CuspidalSubspace(ModularSymbols(MinimalBaseRingCharacter(chi),2,1));
 W:=WindingSubmodule(S);d:=Dimension(S);w:=Dimension(W);
 total54+:=d*EulerPhi(Order(chi));
 printf "OCTIC54_WINDING order=%o conductor=%o dimension=%o winding=%o coefficient_degree=%o time=%o\n",Order(chi),Conductor(chi),d,w,Degree(BaseRing(S)),Cputime(tt54);
 if d ne w then
  Append(~bad54,chi);
  if Order(chi) eq 1 then
   BS:=Matrix([Eltseq(v):v in Basis(S)]);BW:=Matrix([Eltseq(v):v in Basis(W)]);
   VV:=VectorSpace(Rationals(),d);bW:=[VV!Solution(BS,r):r in Rows(BW)];BB:=bW;
   for v in Basis(VV) do if v notin sub<VV|BB> then Append(~BB,v);end if;end for;
   B:=Matrix(BB);
   for q in [2,3,5,7,11] do
    T:=B*HeckeOperator(S,q)*B^-1;
    assert Submatrix(T,1,w+1,w,d-w) eq 0;
    printf "OCTIC54_UPPER_TRIVIAL_HECKE q=%o quotient=%o\n",q,Submatrix(T,w+1,w+1,d-w,d-w);
   end for;
  end if;
 end if;
end for;
assert total54 eq 52;
admissible54:=[a:a in [2..52]|GCD(a,54) eq 1 and &and[Evaluate(chi,a) eq 1:chi in bad54]];
printf "OCTIC54_WINDING_DONE weighted_dimension=%o bad_orders=%o finite_image_diamonds=%o time=%o\n",total54,[Order(chi):chi in bad54],admissible54,Cputime(tt54);
quit;
