AttachSpec("mdmagma/v2/mdmagma.spec");
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(24*10^9);tt65:=Cputime();
X13:=MDX1(13,GF(2));C13:=Curve(X13);F13:=FunctionField(C13);
assert Genus(C13) eq 2;
E13:=EllipticCurve(X13`_E);P13:=E13!X13`_P;
j13scheme:=jInvariant(E13);
K13,m13:=AlgorithmicFunctionField(F13);j13:=m13(j13scheme);
assert &+[-Valuation(j13,p)*Degree(p):p in Poles(j13)] eq 84;
R65<B65>:=PolynomialRing(K13);
f65:=(B65+1)^12+j13*B65^5*(B65^2+B65+1);
assert Degree(f65) eq 12;print "C65_GENUS61_IRREDUCIBILITY_START",Cputime(tt65);assert IsIrreducible(f65);print "C65_GENUS61_IRREDUCIBILITY_DONE",Cputime(tt65);
Y65<b65>:=FunctionField(f65);j65:=Y65!j13;
assert j65 eq (b65+1)^12/(b65^5*(b65^2+b65+1));
print "C65_GENUS61_NORMALIZATION_START",Cputime(tt65);assert Genus(Y65) eq 61;
printf "C65_GENUS61_FIBRE_MODEL_VERIFIED irreducible12=true genus61=true index1008=true time=%o\n",Cputime(tt65);
print "C65_GENUS61_CUSP_START",Cputime(tt65);
basecusps65:=Poles(j13);assert Sort([Degree(p):p in basecusps65]) eq [1:i in [1..6]] cat [6];
cusps65:=[];
for baseIndex65 in [1..#basecusps65] do
 p:=basecusps65[baseIndex65];f,g:=TwoGenerators(p);u:=UniformizingElement(p);
 printf "C65_GENUS61_BASE_CUSP_START index=%o degree=%o time=%o\n",baseIndex65,Degree(p),Cputime(tt65);
 raw:=Decomposition(Y65,p);above:=[q:q in raw|Valuation(Y65!f,q) gt 0 and Valuation(Y65!g,q) gt 0];
 assert &+[Valuation(Y65!u,q)*Degree(q):q in above] eq 12*Degree(p);
 assert &and[Valuation(Y65!u,q) in [1,5]:q in above];
 cusps65 cat:=above;
 printf "C65_GENUS61_BASE_CUSP_DONE index=%o output=%o time=%o\n",baseIndex65,[<Valuation(Y65!u,q),Degree(q)>:q in above],Cputime(tt65);
end for;
assert #Seqset(cusps65) eq #cusps65;
sig65:=Sort([Degree(p):p in cusps65]);
assert sig65 eq [1:i in [1..12]] cat [2:i in [1..6]] cat [6:i in [1..4]];
assert &+[-Valuation(j65,p)*Degree(p):p in cusps65] eq 1008;
printf "C65_GENUS61_CUSPS_VERIFIED count=%o degrees=%o widths=%o time=%o\n",#cusps65,sig65,Sort([<-Valuation(j65,p),Degree(p)>:p in cusps65]),Cputime(tt65);
uv2:=TateNormalForm_xy(E13,2*P13);uvi2:=TateNormalForm_xy(E13,7*P13);
phi13:=iso<C13->C13|[uv2[1],uv2[2],F13!1],[uvi2[1],uvi2[2],F13!1]:Check:=false>;
s13:=map<K13->K13|a:->m13(Pullback(phi13,a@@m13))>;
dg13:=hom<Y65->Y65|map<K13->Y65|a:->Y65!s13(a)>,b65>;
dg5:=hom<Y65->Y65|map<K13->Y65|a:->Y65!a>,1/b65>;
fieldgens65:=[Y65!K13!BaseField(K13).1,Y65!K13.1,b65];
assert &and[(dg5*dg13)(a) eq (dg13*dg5)(a):a in fieldgens65];
assert &and[(dg5^2)(a) eq a and (dg13^6)(a) eq a:a in fieldgens65];
assert &and[(dg5*dg13)(j65) eq j65:a in [1]];
assert #{[ (dg5^i*dg13^k)(a):a in fieldgens65]:i in [0..1],k in [0..5]} eq 12;
print "C65_GENUS61_DECK_GROUP_VERIFIED C2xC6=true order12=true",Cputime(tt65);
Omega65,om65:=SpaceOfHolomorphicDifferentials(Y65);basis65:=[om65(v):v in Basis(Omega65)];
CM65:=Matrix(GF(2),[Eltseq(Cartier(w)@@om65):w in basis65]);I65:=IdentityMatrix(GF(2),61);
eta65:=Differential(fieldgens65[1]);
eta5:=Differential(dg5(fieldgens65[1]));eta13:=Differential(dg13(fieldgens65[1]));
M5:=Matrix(GF(2),[Eltseq((dg5(w/eta65)*eta5)@@om65):w in basis65]);
M13:=Matrix(GF(2),[Eltseq((dg13(w/eta65)*eta13)@@om65):w in basis65]);
assert M5^2 eq I65 and M13^6 eq I65 and M5*M13 eq M13*M5;
assert M5*CM65 eq CM65*M5 and M13*CM65 eq CM65*M13;
U5:=RowSpace((CM65-I65)*(M5-I65));U13:=RowSpace((CM65-I65)*(M13-I65));
U6:=RowSpace((CM65-I65)*(CM65^6-I65));U65:=U5+U13+U6;
printf "C65_GENUS61_COTANGENT_VERIFIED Cartier_rank=%o diamond5_image=%o diamond13_image=%o Hecke6_image=%o joint_image=%o time=%o\n",Rank(CM65),Dimension(U5),Dimension(U13),Dimension(U6),Dimension(U65),Cputime(tt65);
if not assigned MatrixFile65 then MatrixFile65:="C65/data/C65_clean_matrices.m";end if;
Write(MatrixFile65,Sprintf("CartierMatrix:=Matrix(GF(2),61,61,%o);",Eltseq(CM65)) : Overwrite:=true);
Write(MatrixFile65,Sprintf("Diamond5Matrix:=Matrix(GF(2),61,61,%o);",Eltseq(M5)));
Write(MatrixFile65,Sprintf("Diamond13Matrix:=Matrix(GF(2),61,61,%o);",Eltseq(M13)));
Write(MatrixFile65,Sprintf("ImageBasis:=Matrix(GF(2),48,61,%o);",Eltseq(BasisMatrix(U65))));

 tr65:=Cputime();assert #cusps65 eq 22 and Dimension(U65) eq 48;
 s13inverse:=map<K13->K13|a:->m13(Pullback(phi13^-1,a@@m13))>;
 dg13inverse:=hom<Y65->Y65|map<K13->Y65|a:->Y65!s13inverse(a)>,b65>;
 assert &and[dg13(dg13inverse(a)) eq a:a in fieldgens65];
 baseperm65:=[];
 for P in basecusps65 do
  f,g:=TwoGenerators(P);
  hit:=[i:i in [1..#basecusps65]|Degree(basecusps65[i]) eq Degree(P) and Valuation(s13inverse(f),basecusps65[i]) gt 0 and Valuation(s13inverse(g),basecusps65[i]) gt 0];
  assert #hit eq 1;Append(~baseperm65,hit[1]);
 end for;
 assert Sort(baseperm65) eq [1..7];
 qbase65:=[];off65:=0;
 for i in [1..#basecusps65] do
  n:=Degree(basecusps65[i]) eq 1 select 3 else 4;
  for k in [1..n] do
   Append(~qbase65,i);assert Valuation(Y65!UniformizingElement(basecusps65[i]),cusps65[off65+k]) gt 0;
  end for;
  off65+:=n;
 end for;
 assert off65 eq 22;
 qwidth65:=[-Valuation(j65,P):P in cusps65];qdegree65:=[Degree(P):P in cusps65];
 qbval65:=[Valuation(b65,P):P in cusps65];
 pp5:=[];pp13:=[];
 for i in [1..22] do
  same:=[k:k in [1..22]|qbase65[k] eq qbase65[i] and qwidth65[k] eq qwidth65[i]];
  if qbval65[i] ne 0 then
   hs:=[k:k in same|qbval65[k] eq -qbval65[i]];
  elif #same eq 1 then hs:=same;
  else hs:=[k:k in same|k ne i];end if;
  assert #hs eq 1;Append(~pp5,hs[1]);
  otherbase:=baseperm65[qbase65[i]];
  hs:=[k:k in [1..22]|qbase65[k] eq otherbase and qwidth65[k] eq qwidth65[i] and qbval65[k] eq qbval65[i]];
  if #hs ne 1 then
   f,g:=TwoGenerators(cusps65[i]);
   hs:=[k:k in hs|Valuation(dg13inverse(f),cusps65[k]) gt 0 and Valuation(dg13inverse(g),cusps65[k]) gt 0];
  end if;
  assert #hs eq 1;Append(~pp13,hs[1]);
 end for;
 assert Sort(pp5) eq [1..22] and Sort(pp13) eq [1..22];
 assert &and[pp5[pp5[i]] eq i and pp5[pp13[i]] eq pp13[pp5[i]]:i in [1..22]];
 for i in [1..22] do k:=i;for r in [1..6] do k:=pp13[k];end for;assert k eq i;end for;
 allowed65:=[i:i in [1..22]|qwidth65[i] ne 5];
 steps65:=[qdegree65[i] eq 2 select 2 else 1:i in [1..22]];
 remaining65:=Seqset(allowed65);orbits65:=[];
 while #remaining65 gt 0 do
  r:=Minimum(remaining65);orb:={r};
  repeat old:=orb;orb join:={pp5[i]:i in orb} join {pp13[i]:i in orb};until orb eq old;
  assert orb subset remaining65;Append(~orbits65,Sort(Setseq(orb)));remaining65 diff:=orb;
 end while;
 assert Sort([#o:o in orbits65]) eq [2,6,12];
 BM65:=BasisMatrix(U65);assert RowSpace(BM65*M5) eq U65 and RowSpace(BM65*M13) eq U65;
 printf "C65_GENUS61_TRANSPORT_SETUP_DONE orbits=%o source_representatives=%o source_degrees=%o steps=%o time=%o\n",[#o:o in orbits65],[o[1]:o in orbits65],[qdegree65[o[1]]:o in orbits65],[steps65[o[1]]:o in orbits65],Cputime(tr65);
 print "C65_GENUS61_CUSP_PERMUTATIONS",pp5,pp13;
if not assigned MaxFile65 then MaxFile65:="C65/data/C65_clean_max8.m";end if;
 rrmax65:=Cputime();VO65:=VectorSpace(GF(2),61);
 ratrep65:=[o[1]:o in orbits65|qdegree65[o[1]] eq 1][1];
 printf "C65_GENUS61_RR_MAX8_START cusp=%o time=%o\n",ratrep65,Cputime(rrmax65);
 VV65,vmap65:=DifferentialSpace(8*Divisor(cusps65[ratrep65]));
 Wmax65:=sub<VO65|[VO65!Eltseq(vmap65(v)@@om65):v in Basis(VV65)]>;
 Amax65:=Transpose(BasisMatrix(Nullspace(Transpose(BasisMatrix(Wmax65)))));
 assert Nullspace(Amax65) eq Wmax65;
 sourceW65:=AssociativeArray();sourceW65[<ratrep65,8>]:=Wmax65;
 Write(MaxFile65,Sprintf("SourceCusp:=%o;SourceKernel:=Matrix(GF(2),%o,61,%o);RestrictionMatrix:=Matrix(GF(2),61,%o,%o);",ratrep65,Dimension(Wmax65),Eltseq(BasisMatrix(Wmax65)),Ncols(Amax65),Eltseq(Amax65)) : Overwrite:=true);
 printf "C65_GENUS61_RR_MAX8_DONE kernel_dim=%o full_rank=%o rank_on_U=%o time=%o\n",Dimension(Wmax65),Rank(Amax65),Rank(BM65*Amax65),Cputime(rrmax65);
 rt65:=Cputime();VO65:=VectorSpace(GF(2),61);BM65:=BasisMatrix(U65);
 if not assigned JetFile65 then JetFile65:="C65/data/C65_clean_rr_kernels.m";end if;
 Write(JetFile65,"N:=65;D:=8;gx:=61;udim:=48;JetFactor:=1;sourcekernels:=AssociativeArray();kernelspaces:=AssociativeArray();basisjets:=AssociativeArray();localjets:=AssociativeArray();directtargetkernels:=AssociativeArray();" : Overwrite:=true);
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
quit;
