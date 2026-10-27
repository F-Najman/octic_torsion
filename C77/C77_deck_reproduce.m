SetColumns(0);SetQuitOnError(true);SetMemoryLimit(24*10^9);tt:=Cputime();
CQ:=SmallModularCurve(11);jraw:=jInvariant(CQ,11);
CC:=ChangeRing(CQ,GF(2));assert IsNonsingular(CC) and Genus(CC) eq 1;
RC:=CoordinateRing(Ambient(CC));jn:=RC!Numerator(jraw);jd:=RC!Denominator(jraw);
K,kmap:=AlgorithmicFunctionField(FunctionField(CC));jK:=kmap(FunctionField(CC)!(jn/jd));
// The classical X0(7) parameter has j=(v^2+v+1)^4/v^7 in characteristic2.
Fv<v>:=FunctionField(GF(2));j7:=(v^2+v+1)^4/v^7;j7dual:=(v^2+v+1)^4/v;
Phi7:=ClassicalModularPolynomial(7);
assert Evaluate(Phi7,[j7,j7dual]) eq 0;
R<x>:=PolynomialRing(K);f:=x^8+jK*x^7+x^4+1;
assert IsIrreducible(f);AF<u>:=FunctionField(f);assert Genus(AF) eq 7;
jaf:=AF!jK;
cs:=Poles(jaf);widths:=[-Valuation(jaf,c):c in cs];
assert #cs eq 4 and &and[Degree(c) eq 1:c in cs] and Seqset(widths) eq {1,7,11,77};
printf "OCTIC77_X0_MODEL genus7=true modular_polynomial_verified=true cusp_widths=%o time=%o\n",widths,Cputime(tt);
load "C77/C77_AS_component.m";
RF<Y>:=PolynomialRing(AF);f2:=Y^2+Y+a2;assert IsIrreducible(f2);
// Direct cyclic cubic: the fine X1(7) family over its genus-zero parameter.
FB<t>:=FunctionField(GF(2));b7:=t^3-t^2;c7:=t^2-t;
E7u:=EllipticCurve([1-c7,-b7,-b7,FB!0,FB!0]);P7u:=E7u![0,0];assert Order(P7u) eq 7;
ju:=jInvariant(E7u);vu:=(t^2+t)/(t^3+t+1);
assert vu^8+ju*vu^7+vu^4+1 eq 0;
assert Evaluate(Numerator(ju),1/(t+1))/Evaluate(Denominator(ju),1/(t+1)) eq ju;
// j on X011 is j0(x)+j1(x)y, with y quadratic over F2(x).
jparts:=Eltseq(jK);assert #jparts eq 2 and jparts[2] ne 0;
kpoly:=DefiningPolynomial(K);assert Degree(kpoly) eq 2;
FX<X>:=FunctionField(FB);rx:=hom<BaseField(K)->FX|X>;
j0:=rx(jparts[1]);j1:=rx(jparts[2]);aa:=rx(Coefficient(kpoly,1));bb:=rx(Coefficient(kpoly,0));
normrel:=((FX!ju)+j0)^2+aa*j1*((FX!ju)+j0)+bb*j1^2;
f3:=Numerator(normrel);assert Degree(f3) eq 12 and IsIrreducible(f3);
E3<s>:=FunctionField(f3);
kx:=hom<BaseField(K)->E3|s>;
yval:=((E3!ju)+kx(jparts[1]))/kx(jparts[2]);
kemb:=hom<K->E3|kx,yval>;
afemb:=hom<AF->E3|kemb,E3!vu>;
assert afemb(jaf) eq E3!ju;
tau:=hom<E3->E3|hom<FB->E3|1/((E3!t)+1)>,s>;
assert tau(tau(tau(E3!t))) eq E3!t and tau(tau(tau(s))) eq s;
assert &and[tau(afemb(v)) eq afemb(v):v in [AF.1,AF!K.1,AF!K!BaseField(K).1]];
assert Genus(E3) eq 19;
printf "OCTIC77_DIRECT_CUBIC_VERIFIED degree12=true genus19=true exact_base_embedding=true deck_order3=true time=%o\n",Cputime(tt);
QY<z>:=PolynomialRing(E3);fy:=QY![afemb(c):c in Coefficients(f2)];assert IsIrreducible(fy);
CY<zz>:=FunctionField(fy);jY:=CY!E3!ju;
assert Genus(CY) eq 37;
printf "OCTIC77_DIRECT_GENUS37_VERIFIED time=%o\n",Cputime(tt);
cys:=Poles(jY);cywidths:=[-Valuation(jY,P):P in cys];
assert {*<cywidths[i],Degree(cys[i])>:i in [1..#cys]*} eq {*<77,1>^^6,<7,2>^^3,<11,3>^^2,<1,6>^^1*};
assert &+[cywidths[i]*Degree(cys[i]):i in [1..#cys]] eq 576;
printf "OCTIC77_GLOBAL_CUSPS_VERIFIED signature=%o width_degree_pairs=%o jdegree576=true time=%o\n",[Degree(P):P in cys],[<cywidths[i],Degree(cys[i])>:i in [1..#cys]],Cputime(tt);
print "OCTIC77_DIRECT_GLOBAL_GENUS37_MODEL_DONE";

try
 SetColumns(0);dt:=Cputime();assert Genus(CY) eq 37;
 assert &+[cywidths[i]*Degree(cys[i]):i in [1..#cys]] eq 576;
 cs77:=[cys[i]:i in [1..#cys]|cywidths[i] eq 77];assert #cs77 eq 6;
 assert BaseField(CY) cmpeq E3;
 deck2:=hom<CY->CY|map<E3->CY|a:->CY!a>,CY.1+1>;
 deck3:=hom<CY->CY|map<E3->CY|a:->CY!tau(a)>,CY.1>;
 generatorsY:=[CY.1,CY!E3.1,CY!FB.1];
 for a in generatorsY do
  assert deck2(deck2(a)) eq a and deck3(deck3(deck3(a))) eq a;
  assert deck2(deck3(a)) eq deck3(deck2(a));
 end for;
 labels77:=[];
 for P in cs77 do
  tv:=Valuation(CY!FB.1,P) lt 0 select 2 else Integers()!Evaluate(CY!FB.1,P);
  assert tv in [0,1,2] and Valuation(CY.1,P) ge 0;
  zv:=Integers()!Evaluate(CY.1,P);assert zv in [0,1];
  Append(~labels77,<tv,zv>);
 end for;
 assert Seqset(labels77) eq {<i,j>:i in [0,1,2],j in [0,1]};
 pp2:=[Index(labels77,<q[1],1-q[2]>):q in labels77];
 pp3:=[Index(labels77,<q[1] eq 0 select 1 else (q[1] eq 1 select 2 else 0),q[2]>):q in labels77];
 printf "OCTIC77_DECK_LABELS labels=%o perm2=%o perm3=%o time=%o\n",labels77,pp2,pp3,Cputime(dt);
 Omega77,om77:=SpaceOfHolomorphicDifferentials(CY);B77:=[om77(v):v in Basis(Omega77)];
 CM77:=Matrix(GF(2),[Eltseq(Cartier(w)@@om77):w in B77]);
 U77:=RowSpace(CM77*(CM77-IdentityMatrix(GF(2),37)));assert Dimension(U77) eq 32;
 eta77:=Differential(CY!FB.1);hs77:=[w/eta77:w in B77];MM77:=[];
 for k in [1,2] do
  sigma:=k eq 1 select deck2 else deck3;de:=Differential(sigma(CY!FB.1));rows:=[];
  for j in [1..37] do
   w:=sigma(hs77[j])*de;Append(~rows,Eltseq(w@@om77));
   if j mod 5 eq 0 then printf "OCTIC77_DECK_DIFFERENTIAL_PROGRESS deck=%o basis=%o/37 time=%o\n",k,j,Cputime(dt);end if;
  end for;
  Append(~MM77,Matrix(GF(2),rows));
 end for;
 M2:=MM77[1];M3:=MM77[2];I37:=IdentityMatrix(GF(2),37);
 assert M2^2 eq I37 and M3^3 eq I37 and M2*M3 eq M3*M2;
 assert Rank(M2) eq 37 and Rank(M3) eq 37;
 assert M2*CM77 eq CM77*M2 and M3*CM77 eq CM77*M3;
 assert #Seqset([M2^a*M3^b:a in [0,1],b in [0,1,2]]) eq 6;
 assert &and[v*M in U77:v in Basis(U77),M in [M2,M3]];
 printf "OCTIC77_DECK_DIFFERENTIALS_VERIFIED group_order6=true commutes_Cartier=true image_dimension32=true time=%o\n",Cputime(dt);
catch err
 print "OCTIC77_DECK_PREP_FAILED";print err`Object;
end try;

try
 SetColumns(0);rt:=Cputime();gx:=37;D:=8;udim:=32;
 if not assigned DeckJetFile then DeckJetFile:="C77/data/C77_deck_reproduction_scratch_kernels.m";end if;
 VO:=VectorSpace(GF(2),gx);BM:=BasisMatrix(U77);
 Write(DeckJetFile,"// Exact source RR kernels and verified C6 deck transport.\nN:=77;As:=[];D:=8;gx:=37;o:=1;udim:=32;JetFactor:=1;ds:=[1,1,1,1,1,1];perms:=[];\nlocaljets:=AssociativeArray();basisjets:=AssociativeArray();kernelspaces:=AssociativeArray();sourcekernels:=AssociativeArray();directtargetkernels:=AssociativeArray();\nCartierImagePolynomial:=\"C*(C-1)\";" : Overwrite:=true);
 Write(DeckJetFile,Sprintf("CartierMatrix:=Matrix(GF(2),37,37,%o);ImageBasis:=Matrix(GF(2),32,37,%o);",Eltseq(CM77),Eltseq(BM)));
 Write(DeckJetFile,Sprintf("MatrixC2:=Matrix(GF(2),37,37,%o);MatrixC3:=Matrix(GF(2),37,37,%o);",Eltseq(M2),Eltseq(M3)));
 Write(DeckJetFile,Sprintf("CuspLabels:=%o;CuspPermutation2:=%o;CuspPermutation3:=%o;SourceCusp:=1;",labels77,pp2,pp3));
 sourceW:=AssociativeArray();sourceA:=AssociativeArray();cols:=[];annprev:=sub<VO|>;
 for m in [1..8] do
  printf "C77_DECK_RR_SOURCE_START multiplicity=%o time=%o\n",m,Cputime(rt);
  VV,vvmap:=DifferentialSpace(m*Divisor(cs77[1]));
  W:=sub<VO|[VO!Eltseq(vvmap(b)@@om77):b in Basis(VV)]>;
  sourceW[m]:=W;
  ann:=Nullspace(Transpose(BasisMatrix(W)));assert annprev subset ann;
  extra:=[];
  for v in Basis(ann) do
   if v notin annprev then Append(~extra,VO!v);annprev:=sub<VO|[VO!b:b in Basis(annprev)] cat [VO!v]>;end if;
  end for;
  assert #extra le 1;cols cat:=extra;while #cols lt m do Append(~cols,VO!0);end while;
  assert annprev eq ann;A:=Transpose(Matrix(GF(2),[Eltseq(v):v in cols]));assert Nullspace(A) eq W;
  sourceA[m]:=A;
  Write(DeckJetFile,Sprintf("sourcekernels[%o]:=Matrix(GF(2),%o,37,%o);",m,Dimension(W),Eltseq(BasisMatrix(W))));
  printf "C77_DECK_RR_SOURCE_DONE multiplicity=%o kernel_dim=%o rank_on_U=%o time=%o\n",m,Dimension(W),Rank(BM*A),Cputime(rt);
 end for;
 // Physical forward action phi is paired with differential pullback M.
 // Therefore W_{phi(P)}=W_P*M^(-1), and A_{phi(P)}=M*A_P.
 deckpowers:=[];deckmatrices:=[];localjets:=AssociativeArray();
 for i in [1..6] do
  matches:=[];
  for e2 in [0,1] do for e3 in [0,1,2] do
   k:=1;for r in [1..e2] do k:=pp2[k];end for;for r in [1..e3] do k:=pp3[k];end for;
   if k eq i then Append(~matches,<e2,e3>);end if;
  end for;end for;
  assert #matches eq 1;q:=matches[1];M:=M2^q[1]*M3^q[2];
  Append(~deckpowers,q);Append(~deckmatrices,M);
  for m in [1..8] do
   W:=RowSpace(BasisMatrix(sourceW[m])*M^-1);A:=M*sourceA[m];assert Nullspace(A) eq W;
   J:=BM*A;localjets[<i,m>]:=J;
   Write(DeckJetFile,Sprintf("kernelspaces[<%o,%o>]:=Matrix(GF(2),%o,37,%o);",i,m,Dimension(W),Eltseq(BasisMatrix(W))));
   Write(DeckJetFile,Sprintf("localjets[<%o,%o>]:=Matrix(GF(2),32,%o,%o);",i,m,m,Eltseq(J)));
  end for;
  Write(DeckJetFile,Sprintf("basisjets[%o]:=Matrix(GF(2),37,8,%o);",i,Eltseq(M*sourceA[8])));
  if i gt 1 then
   printf "C77_DECK_RR_TARGET_START cusp=%o time=%o\n",i,Cputime(rt);
   VV,vvmap:=DifferentialSpace(Divisor(cs77[i]));
   Wcheck:=sub<VO|[VO!Eltseq(vvmap(b)@@om77):b in Basis(VV)]>;
   assert Wcheck eq RowSpace(BasisMatrix(sourceW[1])*M^-1);
   Write(DeckJetFile,Sprintf("directtargetkernels[%o]:=Matrix(GF(2),%o,37,%o);",i,Dimension(Wcheck),Eltseq(BasisMatrix(Wcheck))));
   printf "C77_DECK_RR_TARGET_VERIFIED cusp=%o time=%o\n",i,Cputime(rt);
  end if;
 end for;
 Write(DeckJetFile,Sprintf("DeckPowers:=%o;",deckpowers));
 Write(DeckJetFile,"print \"FI_JET_DATA_READY\";");
 function Compositions77(d,k)
  if k eq 1 then return [[d]];end if;
  result:=[];for a in [0..d] do for v in Compositions77(d-a,k-1) do Append(~result,[a] cat v);end for;end for;return result;
 end function;
 vs:=Compositions77(8,6);assert #vs eq 1287;bad:=[];
 for idx in [1..#vs] do
  v:=vs[idx];A:=ZeroMatrix(GF(2),32,0);
  for i in [1..6] do if v[i] gt 0 then A:=HorizontalJoin(A,localjets[<i,v[i]>]);end if;end for;
  assert Ncols(A) eq 8;r:=Rank(A);
  if r ne 8 then Append(~bad,v);printf "C77_DECK_FI_EXCEPTION index=%o rank=%o vector=%o\n",idx,r,v;end if;
 end for;
 printf "C77_DECK_FI_DONE tested=%o total=1287 exceptions=%o complete=true time=%o\n",#vs,#bad,Cputime(rt);
catch err
 print "C77_DECK_RR_FAILED";print err`Object;
end try;

quit;
