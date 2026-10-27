SetColumns(0);SetQuitOnError(true);SetMemoryLimit(2*10^9);tt:=Cputime();
N:=57;D:=FullDirichletGroup(N);total:=0;bad:=[];
for chi in GaloisConjugacyRepresentatives(D) do
 if IsOdd(chi) then continue;end if;
 S:=CuspidalSubspace(ModularSymbols(MinimalBaseRingCharacter(chi),2,1));W:=WindingSubmodule(S);
 dim:=Dimension(S);wind:=Dimension(W);total+:=dim*EulerPhi(Order(chi));
 printf "OCTIC57_WINDING order=%o conductor=%o relative_dim=%o winding_dim=%o coefficient_degree=%o time=%o\n",Order(chi),Conductor(chi),dim,wind,Degree(BaseRing(S)),Cputime(tt);
 if dim ne wind then
  Append(~bad,Order(chi));assert Order(chi) eq 1 and dim-wind eq 1;
  BS:=Matrix([Eltseq(b):b in Basis(S)]);BW:=Matrix([Eltseq(b):b in Basis(W)]);T:=HeckeOperator(S,3);
  VS:=VectorSpace(Rationals(),Nrows(BS));ww:=[VS!Solution(BS,r):r in Rows(BW)];UW:=sub<VS|ww>;
  assert &and[r*T in UW:r in ww];bb:=[VS!v:v in Basis(UW)];
  for v in Basis(VS) do if v notin sub<VS|bb> then Append(~bb,v);end if;end for;
  B:=Matrix(bb);M:=B*T*B^-1;Q:=Submatrix(M,wind+1,wind+1,1,1);assert Q eq -1;
  print "OCTIC57_POSITIVE_UPPER_U3_EXACTLY_MINUS1",Eltseq(Q);
 end if;
end for;
assert total eq 85 and bad eq [1];
print "OCTIC57_FINITE_IMAGE_U3_PLUS1_CERTIFIED",total,bad,Cputime(tt);quit;
