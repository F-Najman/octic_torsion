SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt:=Cputime();
try
 CQ:=SmallModularCurve(19);jraw:=jInvariant(CQ,19);
 CC:=ChangeRing(CQ,GF(2));assert IsNonsingular(CC) and Genus(CC) eq 1 and #CC eq 3;
 RC:=CoordinateRing(Ambient(CC));jn:=RC!Numerator(jraw);jd:=RC!Denominator(jraw);
 K,kmap:=AlgorithmicFunctionField(FunctionField(CC));jK:=kmap(FunctionField(CC)!(jn/jd));
 xK:=K!BaseField(K).1;yK:=K.1;EK:=ChangeRing(CC,K);Pgen:=EK![xK,yK,1];
 c19:=[P:P in Poles(jK)|-Valuation(jK,P) eq 19];assert #c19 eq 1;c19:=c19[1];
 if Valuation(xK,c19) lt 0 then origin19:=CC!0;
 else origin19:=CC![Evaluate(xK,c19),Evaluate(yK,c19),1];end if;
 Lgen:=EK![xK^2,yK^2,1]-Pgen+EK!origin19;
 langbase:=hom<BaseField(K)->K|Lgen[1]/Lgen[3]>;
 langK:=hom<K->K|langbase,Lgen[2]/Lgen[3]>;
 assert langK(jK) ne jK;
 point3:=[P:P in Points(CC)|P ne CC!0][1];assert Order(point3) eq 3;
 Rgen:=Pgen+EK!point3;
 deckbase:=hom<BaseField(K)->K|Rgen[1]/Rgen[3]>;
 deckK:=hom<K->K|deckbase,Rgen[2]/Rgen[3]>;
 assert &and[deckK(deckK(deckK(a))) eq a:a in [xK,yK]];
 assert &and[deckK(langK(a)) eq langK(a):a in [xK,yK]];
 langj:=langK(jK);
 assert &+[-Valuation(langj,P)*Degree(P):P in Poles(langj)] eq 60;
 printf "C57_LANG_ISOGENY_VERIFIED degree3=true cyclic3=true split_width19_origin=%o time=%o\n",origin19,Cputime(tt);
 R<v>:=PolynomialRing(K);f:=v^4+v^3+1/jK;assert IsIrreducible(f);
 AF<t3>:=FunctionField(f);assert Genus(AF) eq 5;jAF:=AF!jK;
 cs:=Poles(jAF);widths:=[-Valuation(jAF,P):P in cs];assert #cs eq 4 and &and[Degree(P) eq 1:P in cs] and Seqset(widths) eq {1,3,19,57};
 g:=R![langK(a):a in Coefficients(f)];assert IsIrreducible(g);
 E3<s3>:=FunctionField(g);assert Genus(E3) eq 13;
 afemb:=hom<AF->E3|map<K->E3|a:->E3!langK(a)>,s3>;
 tau:=hom<E3->E3|map<K->E3|a:->E3!deckK(a)>,s3>;
 assert &and[tau(afemb(a)) eq afemb(a):a in [AF.1,AF!xK,AF!yK]];
 printf "C57_CUBIC_PULLBACK_VERIFIED genus13=true exact_base_embedding=true cyclic3=true time=%o\n",Cputime(tt);
 DG:=DivisorGroup(AF);DZ:=DG!0;G,gm:=RayClassGroup(DZ);inv:=Invariants(G);assert inv eq [15,30,0];
 cls:=[Eltseq(Divisor(P)@@gm):P in cs];expected:=[w in [3,19] select 1 else 0:w in widths];rr:=[];
 for r in CartesianPower({0,1},Ngens(G)) do
  a:=[r[i]:i in [1..Ngens(G)]];
  if &and[inv[i]*a[i] mod 2 eq 0:i in [1..#inv]] and &and[(&+[a[k]*cls[i][k]:k in [1..#inv]]) mod 2 eq expected[i]:i in [1..4]] then Append(~rr,a);end if;
 end for;
 assert #rr eq 1;C2:=AbelianGroup([2]);h2:=hom<G->C2|[C2![a]:a in rr[1]]>;
 printf "C57_QUADRATIC_CHARACTER_UNIQUE invariants=%o images=%o widths=%o expected=%o time=%o\n",inv,rr[1],widths,expected,Cputime(tt);
// Exact AS coefficient over AF.
x0:=AF!K!BaseField(K).1;y0:=AF!K.1;u0:=AF.1;
aAS:=((x0^11 + x0^10 + x0^9 + x0^2 + x0)/(x0^11 + x0^10 + x0^9 + x0^8 + x0^3 + x0^2 + x0 + 1)*y0 + (x0^12 + x0^11 + x0^8 + x0^6 + x0^5 + x0^3 + 1)/(x0^8 + 1))*u0^3 + ((x0^13 + x0^12 + x0^10 + x0^9 + x0^8 + x0^7 + x0^6 + x0^5 + x0^4 + x0 + 1)/(x0^11 + x0^10 + x0^9 + x0^8 + x0^3 + x0^2 + x0 + 1)*y0 + (x0^11 + x0^8 + x0^7 + x0^5 + x0^4 + x0^2 + 1)/(x0^8 + 1))*u0^2 + ((x0^8 + x0^6 + x0^5 + x0 + 1)/(x0^6 + x0^4 + x0^2 + 1)*y0 + (x0^5 + x0^4 + x0^2 + x0 + 1)/(x0^3 + x0^2 + x0 + 1))*u0 + (x0^5 + x0^3 + x0)/(x0^3 + x0^2 + x0 + 1)*y0 + (x0^6 + x0^4 + x0)/(x0^3 + x0^2 + x0 + 1);assert Parent(aAS) cmpeq AF;

 RY<z>:=PolynomialRing(E3);fy:=z^2+z+afemb(aAS);assert IsIrreducible(fy);
 CY<zz>:=FunctionField(fy);jY:=CY!afemb(jAF);assert Genus(CY) eq 25;
 printf "C57_GENUS25_MODEL_VERIFIED good_modular_identification_pending=true time=%o\n",Cputime(tt);
 print "C57_FIXED_AS_SMALL_MODEL_READY";
catch err
 print "C57_FIXED_AS_MODEL_FAILED";print err`Object;
end try;

try
 G2,q2:=quo<G|[2*g:g in Generators(G)]>;span2:=sub<G2|>;pe:=PlaceEnumInit(AF);num57:=0;tried57:=0;
 repeat
  PP:=PlaceEnumNext(pe);tried57+:=1;assert tried57 le 200;
  if Valuation(aAS,PP) lt 0 then continue;end if;
  FF:=ResidueClassField(PP);tr:=Integers()!Trace(FF!Evaluate(aAS,PP),GF(2));cl:=Divisor(PP)@@gm;
  assert tr eq Eltseq(h2(cl))[1];span2:=sub<G2|span2,q2(cl)>;num57+:=1;
  printf "C57_AS_FROBENIUS degree=%o class=%o value=%o generated=%o total=%o\n",Degree(PP),Eltseq(cl),tr,#span2,#G2;
 until span2 eq G2;
 printf "C57_AS_CHARACTER_VERIFIED samples=%o all_Picard_mod2=true\n",num57;
 print "C57_MODULAR_GENUS25_VERIFIED";
catch err
 print "C57_AS_VALIDATION_FAILED";print err`Object;
end try;

Omega,om:=SpaceOfHolomorphicDifferentials(CY);B:=[om(v):v in Basis(Omega)];CM:=Matrix(GF(2),[Eltseq(Cartier(w)@@om):w in B]);U:=RowSpace(CM*(CM-IdentityMatrix(GF(2),25)));

try
 if not assigned QExpFile then QExpFile:="C57/data/C57_joint_reproduction_qexp_scratch.m";end if;
 qt:=Cputime();checkn:=20;need:=3*checkn+2;workp:=256;
 F8<a8>:=GF(8);EF8:=ChangeRing(CC,F8);
 c1:=[P:P in Poles(jK)|-Valuation(jK,P) eq 1];assert #c1 eq 1;c1:=c1[1];
 if Valuation(xK,c1) lt 0 then target1:=CC!0;else target1:=CC![Evaluate(xK,c1),Evaluate(yK,c1),1];end if;
 src:=[];
 for R in Points(EF8) do
  FR:=EF8![x^2:x in Eltseq(R)];
  if FR-R+EF8!origin19 eq EF8!target1 then Append(~src,R);end if;
 end for;
 assert #src eq 3 and &and[R[3] ne 0:R in src];Rsrc:=src[1];ax:=Rsrc[1]/Rsrc[3];yResid:=Rsrc[2]/Rsrc[3];
 assert Degree(MinimalPolynomial(ax,GF(2))) eq 3;
 LS<zloc>:=LaurentSeriesRing(F8,workp);xs:=ax+zloc+O(zloc^workp);
 function EvalRat57(r,xser)
  return Evaluate(PolynomialRing(F8)!Numerator(r),xser)/Evaluate(PolynomialRing(F8)!Denominator(r),xser);
 end function;
 kpoly:=DefiningPolynomial(K);assert Coefficient(kpoly,1) eq 1;
 rhs:=EvalRat57(Coefficient(kpoly,0),xs);rpos:=rhs-(yResid^2+yResid);assert Valuation(rpos) ge 1;
 ys:=LS!yResid;term:=rpos;
 while Valuation(term) lt workp do ys+:=term;term:=term^2;end while;
 assert Valuation(ys^2+ys+rhs) ge need;
 function EvalK57(a)
  co:=Eltseq(K!a);return &+[EvalRat57(co[i],xs)*ys^(i-1):i in [1..#co]];
 end function;
 jinv:=1/EvalK57(langj);assert Valuation(jinv) eq 1 and AbsolutePrecision(jinv) ge need;
 vs:=LS!1+O(zloc^workp);
 for step in [1..12] do vs+:=(vs^4+vs^3+jinv)/vs^2;end for;
 assert Valuation(vs-1) eq 1 and Valuation(vs^4+vs^3+jinv) ge need;
 function EvalE357(a)
  co:=Eltseq(E3!a);return &+[EvalK57(co[i])*vs^(i-1):i in [1..#co]];
 end function;
 aser:=EvalE357(afemb(aAS));hsAS:=LS!0;
 while Valuation(aser) lt 0 do
  v:=Valuation(aser);assert IsEven(v);hh:=Coefficient(aser,v)^4*zloc^(v div 2);
  hsAS+:=hh;aser+:=hh^2+hh;
 end while;
 R8<X8>:=PolynomialRing(F8);rootsAS:=Roots(X8^2+X8+Coefficient(aser,0));assert #rootsAS eq 2;
 z0:=rootsAS[1][1];zs:=LS!z0;term:=aser-(z0^2+z0);
 assert Valuation(term) ge 1;
 while Valuation(term) lt workp do zs+:=term;term:=term^2;end while;
 zs+:=hsAS;
 assert Valuation(zs^2+zs+EvalE357(afemb(aAS))) ge need;
 function EvalCY57(a)
  co:=Eltseq(CY!a);return &+[EvalE357(co[i])*zs^(i-1):i in [1..#co]];
 end function;
 printf "C57_EXPLICIT_HENSEL_EMBEDDING_VERIFIED K=true cubic=true AS=true source_point=%o precision=%o time=%o\n",Rsrc,workp,Cputime(qt);
 ss:=vs+1;local_inverse:=Reversion(ss);assert AbsolutePrecision(local_inverse) ge need;
 qQ:=zloc;Delta:=qQ+O(qQ^workp);
 for n in [1..workp div 8] do Delta*:=(1+qQ^(8*n))^3;end for;
 targetS:=Evaluate(Reversion(qQ*(1+qQ)^3+O(qQ^workp)),Delta);
 assert Valuation(targetS*(1+targetS)^3-Delta) ge need;dS:=Derivative(targetS);
 etaX:=Differential(CY!E3!xK);formsQ:=[];
 for i in [1..25] do
  printf "C57_EXPLICIT_U3_FORM_START basis=%o time=%o\n",i,Cputime(qt);
  hi:=EvalCY57(B[i]/etaX);assert Valuation(hi) ge 0 and AbsolutePrecision(hi) ge need;
  // omega=hi(t) dt and t=t(s); then s=s(q).
  hs:=Evaluate(hi,local_inverse)*Derivative(local_inverse);
  fq:=qQ*dS*Evaluate(hs,targetS);assert AbsolutePrecision(fq) ge need and Valuation(fq) ge 1;
  Append(~formsQ,fq);
  printf "C57_EXPLICIT_U3_FORM_DONE basis=%o precision=%o time=%o\n",i,AbsolutePrecision(fq),Cputime(qt);
 end for;
 Bq:=Matrix(GF(2),[&cat[Eltseq(Coefficient(fq,n),GF(2)):n in [1..checkn]]:fq in formsQ]);
 Tq:=Matrix(GF(2),[&cat[Eltseq(Coefficient(fq,3*n),GF(2)):n in [1..checkn]]:fq in formsQ]);
 assert Rank(Bq) eq 25;assert &and[r in RowSpace(Bq):r in Rows(Tq)];
 H3:=Matrix(GF(2),[Eltseq(Solution(Bq,r)):r in Rows(Tq)]);assert H3*Bq eq Tq and H3*CM eq CM*H3;
 Cq:=Matrix(GF(2),[&cat[Eltseq(Coefficient(fq,2*n)^4,GF(2)):n in [1..checkn]]:fq in formsQ]);assert CM*Bq eq Cq;
 U3extra:=RowSpace((CM-IdentityMatrix(GF(2),25))*(H3+IdentityMatrix(GF(2),25)));Uplus3:=U+U3extra;
 printf "C57_EXPLICIT_U3_MATRIX_VERIFIED coefficients=%o F2matrix=true Cartier_normalization=true commutesCartier=true oldU=%o extra=%o joint=%o time=%o\n",checkn,Dimension(U),Dimension(U3extra),Dimension(Uplus3),Cputime(qt);
 Write(QExpFile,Sprintf("N:=57;gx:=25;CheckedCoefficients:=%o;OriginalCoefficients:=%o;",checkn,3*checkn):Overwrite:=true);
 Write(QExpFile,Sprintf("CartierMatrix:=Matrix(GF(2),25,25,%o);U3Matrix:=Matrix(GF(2),25,25,%o);",Eltseq(CM),Eltseq(H3)));
 Write(QExpFile,Sprintf("SourceQExp:=Matrix(GF(2),25,%o,%o);TargetQExp:=Matrix(GF(2),25,%o,%o);",Ncols(Bq),Eltseq(Bq),Ncols(Tq),Eltseq(Tq)));
 fullQ:=Matrix(GF(2),[&cat[Eltseq(Coefficient(fq,n),GF(2)):n in [1..3*checkn]]:fq in formsQ]);
 Write(QExpFile,Sprintf("FullSourceQExp:=Matrix(GF(2),25,%o,%o);F8ModulusCoefficients:=%o;",Ncols(fullQ),Eltseq(fullQ),[Integers()!a:a in Eltseq(DefiningPolynomial(F8))]));
 Write(QExpFile,Sprintf("ImageBasis:=Matrix(GF(2),%o,25,%o);",Dimension(Uplus3),Eltseq(BasisMatrix(Uplus3))));
 print "C57_EXPLICIT_U3_DONE";
catch err
 print "C57_EXPLICIT_U3_FAILED";print err`Object;
end try;

cys:=Poles(jY);cywidths:=[-Valuation(jY,P):P in cys];assert {*<cywidths[i],Degree(cys[i])>:i in [1..#cys]*} eq {*<57,1>^^6,<19,2>^^3,<3,6>^^1,<1,3>^^2*};print "C57_FIXED_AS_GLOBAL_CUSPS_VERIFIED";
try
SetColumns(0);SetQuitOnError(true);N:=57;D:=8;gx:=25;gy:=25;t0:=Cputime();Limit:=0;DirectCheck:=1;StopAfterExceptions:=1;if not assigned JetFile then JetFile:="C57/data/C57_joint_reproduction_scratch_kernels.m";end if;
AFY:=CY;
U:=Uplus3;udim:=Dimension(U);assert udim eq 23;
printf "C57_JOINT_CARTIER_FORMS N=%o image_dimension=%o time=%o\n",N,udim,Cputime(t0);
if udim lt D then print "C57_JOINT_CARTIER_INSUFFICIENT_COTANGENT_DIMENSION"; quit; end if;
cs:=[cys[i]:i in [1..#cys]|cywidths[i] in [57,19]];ds:=[Degree(p):p in cs];
assert {*d:d in ds*} eq {*1^^6,2^^3*};
function WeightedVectors(ds,d)
    if #ds eq 0 then return d eq 0 select [[]] else []; end if;
    result:=[];
    for k in [0..d div ds[1]] do
        for v in WeightedVectors(ds[2..#ds],d-k*ds[1]) do Append(~result,[k] cat v); end for;
    end for;
    return result;
end function;
vs:=WeightedVectors(ds,D);
printf "C57_JOINT_CARTIER_COUNT N=%o total=%o time=%o\n",N,#vs,Cputime(t0);
assert #vs eq 3654;
/* Exact local restriction kernels from DifferentialSpace, no series. */

Write(JetFile,Sprintf("// Direct differential-space restriction certificate.\nN:=%o; As:=[]; D:=%o; gx:=%o; o:=1; udim:=%o; JetFactor:=1; ds:=%o; perms:=[];\nlocaljets:=AssociativeArray(); basisjets:=AssociativeArray(); kernelspaces:=AssociativeArray();",N,D,gx,udim,ds):Overwrite:=true);
Write(JetFile,"JointHeckeImage:=true;CartierImagePolynomial:=\"C*(C-1)+(U3+1)*(C-1)\";");
Write(JetFile,Sprintf("U3Matrix:=Matrix(GF(2),25,25,%o);",Eltseq(H3)));
Write(JetFile,Sprintf("CartierMatrix:=Matrix(GF(2),%o,%o,%o); ImageBasis:=Matrix(GF(2),%o,%o,%o);",gx,gx,Eltseq(CM),udim,gx,Eltseq(BasisMatrix(U))));
VO:=VectorSpace(GF(2),gx);localjets:=AssociativeArray();
for i in [1..#cs] do
    maxm:=D div ds[i];cols:=[];annprev:=sub<VO|>;
    for m in [1..maxm] do
        printf "C57_JOINT_RR_KERNEL_START cusp=%o/%o multiplicity=%o time=%o\n",i,#cs,m,Cputime(t0);
        VV,vvmap:=DifferentialSpace(m*Divisor(cs[i]));
        W:=sub<VO|[VO!Eltseq(vvmap(b) @@ om):b in Basis(VV)]>;
        Write(JetFile,Sprintf("kernelspaces[<%o,%o>]:=Matrix(GF(2),%o,%o,%o);",i,m,Dimension(W),gx,Eltseq(BasisMatrix(W))));
        ann:=Nullspace(Transpose(BasisMatrix(W)));assert annprev subset ann;
        extra:=[];
        for v in Basis(ann) do
            if v notin annprev then Append(~extra,v);annprev:=sub<VO|Basis(annprev) cat [v]>;end if;
        end for;
        assert #extra le ds[i];cols cat:=extra;
        while #cols lt m*ds[i] do Append(~cols,VO!0);end while;
        assert annprev eq ann;
        A:=Transpose(Matrix(GF(2),[Eltseq(v):v in cols]));
        assert Nullspace(A) eq W;
        J:=BasisMatrix(U)*A;localjets[<i,m>]:=J;
        Write(JetFile,Sprintf("localjets[<%o,%o>]:=Matrix(GF(2),%o,%o,%o);",i,m,Nrows(J),Ncols(J),Eltseq(J)));
        printf "C57_JOINT_RR_KERNEL_DONE cusp=%o/%o multiplicity=%o kernel_dim=%o rank_on_U=%o time=%o\n",i,#cs,m,Dimension(W),Rank(J),Cputime(t0);
    end for;
    Write(JetFile,Sprintf("basisjets[%o]:=Matrix(GF(2),%o,%o,%o);",i,Nrows(A),Ncols(A),Eltseq(A)));
end for;
Write(JetFile,"print \"FI_JET_DATA_READY\";");

VO:=VectorSpace(GF(2),gx); bad:=[]; tested:=0;
for v in vs do
    if Limit gt 0 and tested ge Limit then break; end if;
    A:=ZeroMatrix(GF(2),udim,0);
    for i in [1..#cs] do if v[i] gt 0 then A:=HorizontalJoin(A,localjets[<i,v[i]>]); end if; end for;
    assert Ncols(A) eq D; r:=Rank(A);
    if tested lt DirectCheck or r ne D then
        Ediv:=&+[v[i]*Divisor(cs[i]):i in [1..#cs]]; assert Degree(Ediv) eq D;
        V,vm:=DifferentialSpace(Ediv);
        V_original:=sub<VO|[VO!Eltseq(vm(b) @@ om):b in Basis(V)]>;
        direct_rank:=udim-Dimension(U meet V_original); assert direct_rank eq r;
        printf "C57_JOINT_CARTIER_DIRECT N=%o index=%o jet_rank=%o direct_rank=%o time=%o\n",N,tested+1,r,direct_rank,Cputime(t0);
    end if;
    tested+:=1;
    if r ne D then Append(~bad,v); printf "C57_JOINT_CARTIER_EXCEPTION N=%o index=%o rank=%o vector=%o time=%o\n",N,tested,r,v,Cputime(t0); end if;
    if tested mod 1000 eq 0 then printf "C57_JOINT_CARTIER_PROGRESS N=%o tested=%o total=%o exceptions=%o time=%o\n",N,tested,#vs,#bad,Cputime(t0); end if;
    if StopAfterExceptions gt 0 and #bad ge StopAfterExceptions then print "C57_JOINT_CARTIER_EXCEPTION_LIMIT_REACHED"; break; end if;
end for;
printf "C57_JOINT_CARTIER_DONE N=%o D=%o tested=%o total=%o exceptions=%o complete=%o time=%o\n",N,D,tested,#vs,#bad,tested eq #vs,Cputime(t0);

catch err
 print "C57_JOINT_CARTIER_RR_FAILED";print err`Object;
end try;

quit;
