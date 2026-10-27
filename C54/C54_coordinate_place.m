/* Exact coordinate ideal whenever at least one coordinate generates its
   residue field. In particular this holds for a full pair of degree8,
   because all proper subfields of F_(5^8) lie in F_(5^4). */
intrinsic MDPlaceCoordinateFullDegree(coordinates::[FldFunFracSchElt],values::[FldFinElt]) -> PlcCrvElt
{ Recover a normalized place from a coordinate which generates the residue field }
 require #coordinates eq 2 and #values eq 2 : "Require two coordinates and values";
 x,y:=Explode(values);K:=Parent(x);require Parent(y) cmpeq K : "Require common explicit field";
 F:=PrimeField(K);d:=Degree(K,F);fx:=MinimalPolynomial(x,F);fy:=MinimalPolynomial(y,F);
 require Degree(fx) eq d or Degree(fy) eq d : "One coordinate must generate the residue field";
 R<u,v>:=PolynomialRing(F,2);
 if Degree(fx) eq d then
  M:=Matrix(F,[Eltseq(x^i,F):i in [0..d-1]]);assert Rank(M) eq d;
  coeff:=Solution(M,Vector(F,Eltseq(y,F)));assert &+[K!coeff[i+1]*x^i:i in [0..d-1]] eq y;
  f:=Evaluate(fx,u);g:=v-&+[coeff[i+1]*u^i:i in [0..d-1]];
 else
  M:=Matrix(F,[Eltseq(y^i,F):i in [0..d-1]]);assert Rank(M) eq d;
  coeff:=Solution(M,Vector(F,Eltseq(x,F)));assert &+[K!coeff[i+1]*y^i:i in [0..d-1]] eq x;
  f:=Evaluate(fy,v);g:=u-&+[coeff[i+1]*v^i:i in [0..d-1]];
 end if;
 assert Evaluate(f,[x,y]) eq 0 and Evaluate(g,[x,y]) eq 0;
 f1:=Evaluate(f,coordinates);g1:=Evaluate(g,coordinates);FF:=Parent(coordinates[1]);C:=Curve(FF);AF,ma:=AlgorithmicFunctionField(FF);
 O:=MaximalOrderFinite(AF);I:=ideal<O|ma(f1),ma(g1)>;fs:=Factorization(I);ys:=[CurvePlace(C,Place(z[1])):z in fs|z[2] gt 0];
 if #ys eq 0 then
  Oi:=MaximalOrderInfinite(AF);Ii:=ideal<Oi|ma(f1),ma(g1)>;fsi:=Factorization(Ii);ys:=[CurvePlace(C,Place(z[1])):z in fsi|z[2] gt 0];
 end if;
 assert #ys eq 1;P:=ys[1];assert Degree(P) eq d and Valuation(f1,P) gt 0 and Valuation(g1,P) gt 0;
 return P;
end intrinsic;
