try
 rt65:=Cputime();VO65:=VectorSpace(GF(2),61);BM65:=BasisMatrix(U65);
 JetFile65:="C65/data/C65_root_allbad_rr_kernels.m";
 Write(JetFile65,"N:=65;D:=8;gx:=61;udim:=48;JetFactor:=1;sourcekernels:=AssociativeArray();kernelspaces:=AssociativeArray();basisjets:=AssociativeArray();localjets:=AssociativeArray();directtargetkernels:=AssociativeArray();");
 Write(JetFile65,Sprintf("CartierMatrix:=Matrix(GF(2),61,61,%o);ImageBasis:=Matrix(GF(2),48,61,%o);",Eltseq(CM65),Eltseq(BM65)));
 Write(JetFile65,Sprintf("Diamond5Matrix:=Matrix(GF(2),61,61,%o);Diamond13Matrix:=Matrix(GF(2),61,61,%o);",Eltseq(M5),Eltseq(M13)));
 Write(JetFile65,Sprintf("ds:=%o;CuspWidths:=%o;CuspBaseIndices:=%o;CuspBValuations:=%o;BasePermutation13:=%o;CuspPermutation5:=%o;CuspPermutation13:=%o;AllowedCusps:=%o;",qdegree65,qwidth65,qbase65,qbval65,baseperm65,pp5,pp13,allowed65));
 sourceA65:=AssociativeArray();
 reps65:=[ratrep65] cat [o[1]:o in orbits65|qdegree65[o[1]] ne 1];
 assert #reps65 eq 3;
 for rep in reps65 do
  deg:=qdegree65[rep];step:=steps65[rep];cols:=[];annprev:=sub<VO65|>;
  for m in [step..(8 div deg) by step] do
   key:=<rep,m>;
   printf "C65_RR_SOURCE_START cusp=%o degree=%o multiplicity=%o time=%o\n",rep,deg,m,Cputime(rt65);
   if IsDefined(sourceW65,key) then W:=sourceW65[key];
   else
    VV,vm:=DifferentialSpace(m*Divisor(cusps65[rep]));
    W:=sub<VO65|[VO65!Eltseq(vm(v)@@om65):v in Basis(VV)]>;sourceW65[key]:=W;
   end if;
   ann:=Nullspace(Transpose(BasisMatrix(W)));assert annprev subset ann;
   extra:=[];
   for v in Basis(ann) do
    if v notin annprev then Append(~extra,VO65!v);annprev:=sub<VO65|[VO65!a:a in Basis(annprev)] cat [VO65!v]>;end if;
   end for;
   assert #extra le deg*step;cols cat:=extra;while #cols lt m*deg do Append(~cols,VO65!0);end while;
   A:=Transpose(Matrix(GF(2),[Eltseq(v):v in cols]));assert Ncols(A) eq m*deg and Nullspace(A) eq W;
   sourceA65[key]:=A;
   Write(JetFile65,Sprintf("sourcekernels[<%o,%o>]:=Matrix(GF(2),%o,61,%o);",rep,m,Dimension(W),Eltseq(BasisMatrix(W))));
   printf "C65_RR_SOURCE_DONE cusp=%o multiplicity=%o kernel_dim=%o U_rank=%o time=%o\n",rep,m,Dimension(W),Rank(BM65*A),Cputime(rt65);
  end for;
 end for;
 local65:=AssociativeArray();transport65:=[];repOf65:=[];
 for i in allowed65 do
  rep:=[r:r in reps65|qdegree65[r] eq qdegree65[i]][1];
  matches:=[];
  for a in [0,1] do for b in [0..5] do
   k:=rep;for r in [1..a] do k:=pp5[k];end for;for r in [1..b] do k:=pp13[k];end for;
   if k eq i then Append(~matches,<a,b>);end if;
  end for;end for;
  assert #matches gt 0;power:=matches[1];M:=M5^power[1]*M13^power[2];
  Append(~transport65,[i,rep,power[1],power[2]]);Append(~repOf65,rep);
  for m in [steps65[i]..(8 div qdegree65[i]) by steps65[i]] do
   W:=RowSpace(BasisMatrix(sourceW65[<rep,m>])*M^-1);A:=M*sourceA65[<rep,m>];
   assert Nullspace(A) eq W;
   for power2 in matches do
    M2:=M5^power2[1]*M13^power2[2];assert RowSpace(BasisMatrix(sourceW65[<rep,m>])*M2^-1) eq W;
   end for;
   J:=BM65*A;local65[<i,m>]:=J;
   Write(JetFile65,Sprintf("kernelspaces[<%o,%o>]:=Matrix(GF(2),%o,61,%o);basisjets[<%o,%o>]:=Matrix(GF(2),61,%o,%o);localjets[<%o,%o>]:=Matrix(GF(2),48,%o,%o);",i,m,Dimension(W),Eltseq(BasisMatrix(W)),i,m,Ncols(A),Eltseq(A),i,m,Ncols(J),Eltseq(J)));
  end for;
 end for;
 Write(JetFile65,Sprintf("SourceRepresentatives:=%o;DeckTransports:=%o;",reps65,transport65));
 for orb in orbits65 do
  target:=[i:i in orb|i ne orb[1]][1];m:=steps65[target];
  VV,vm:=DifferentialSpace(m*Divisor(cusps65[target]));
  Wdirect:=sub<VO65|[VO65!Eltseq(vm(v)@@om65):v in Basis(VV)]>;
  transport:=[v:v in transport65|v[1] eq target][1];rep:=transport[2];M:=M5^transport[3]*M13^transport[4];
  assert Wdirect eq RowSpace(BasisMatrix(sourceW65[<rep,m>])*M^-1);
  Write(JetFile65,Sprintf("directtargetkernels[<%o,%o>]:=Matrix(GF(2),%o,61,%o);",target,m,Dimension(Wdirect),Eltseq(BasisMatrix(Wdirect))));
  printf "C65_RR_DIRECT_TARGET_CHECK cusp=%o multiplicity=%o time=%o\n",target,m,Cputime(rt65);
 end for;
 multiplicities65:=[i in allowed65 select [0..(8 div qdegree65[i]) by steps65[i]] else [0]:i in [1..22]];
 Write(JetFile65,Sprintf("AdmissibleMultiplicities:=%o;",multiplicities65));
 Write(JetFile65,"print \"C65_RR_KERNEL_DATA_READY\";");
 function Divisors65(d,i,ds,mults)
  if i gt #ds then return d eq 0 select [[]] else [];end if;
  out:=[];
  for m in mults[i] do
   if m*ds[i] le d then for v in Divisors65(d-m*ds[i],i+1,ds,mults) do Append(~out,[m] cat v);end for;end if;
  end for;
  return out;
 end function;
 vectors65:=Divisors65(8,1,qdegree65,multiplicities65);assert #vectors65 eq 83949;
 bad65:=[];tested65:=0;
 for vector in vectors65 do
  A:=ZeroMatrix(GF(2),48,0);
  for i in allowed65 do if vector[i] gt 0 then A:=HorizontalJoin(A,local65[<i,vector[i]>]);end if;end for;
  assert Ncols(A) eq 8;rank:=Rank(A);tested65+:=1;
  if rank ne 8 then Append(~bad65,vector);printf "C65_FI_EXCEPTION index=%o rank=%o vector=%o\n",tested65,rank,vector;break;end if;
 end for;
 printf "C65_GENUS61_FI_DONE tested=%o total=83949 exceptions=%o complete=%o time=%o\n",tested65,#bad65,tested65 eq 83949,Cputime(rt65);
catch err
 print "C65_GENUS61_RR_TRANSPORT_FAILED";print err`Object;
end try;
