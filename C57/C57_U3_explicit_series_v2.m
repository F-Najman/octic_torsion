try
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
 Write("C57/data/C57_U3_explicit_data.m",Sprintf("N:=57;gx:=25;CheckedCoefficients:=%o;OriginalCoefficients:=%o;",checkn,3*checkn):Overwrite:=true);
 Write("C57/data/C57_U3_explicit_data.m",Sprintf("CartierMatrix:=Matrix(GF(2),25,25,%o);U3Matrix:=Matrix(GF(2),25,25,%o);",Eltseq(CM),Eltseq(H3)));
 Write("C57/data/C57_U3_explicit_data.m",Sprintf("SourceQExp:=Matrix(GF(2),25,%o,%o);TargetQExp:=Matrix(GF(2),25,%o,%o);",Ncols(Bq),Eltseq(Bq),Ncols(Tq),Eltseq(Tq)));
 fullQ:=Matrix(GF(2),[&cat[Eltseq(Coefficient(fq,n),GF(2)):n in [1..3*checkn]]:fq in formsQ]);
 Write("C57/data/C57_U3_explicit_data.m",Sprintf("FullSourceQExp:=Matrix(GF(2),25,%o,%o);F8ModulusCoefficients:=%o;",Ncols(fullQ),Eltseq(fullQ),[Integers()!a:a in Eltseq(DefiningPolynomial(F8))]));
 Write("C57/data/C57_U3_explicit_data.m",Sprintf("ImageBasis:=Matrix(GF(2),%o,25,%o);",Dimension(Uplus3),Eltseq(BasisMatrix(Uplus3))));
 print "C57_EXPLICIT_U3_DONE";
catch err
 print "C57_EXPLICIT_U3_FAILED";print err`Object;
end try;
