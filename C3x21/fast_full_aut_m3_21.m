/* Canonical relative-base full automorphisms of X1(3,21). */
intrinsic M3TwentyOneNormalizedPair(E::CrvEll,P::PtEll,Q::PtEll) -> SeqEnum
{ Canonical q,t coordinates for an independent pair of orders3,21 }
 assert P ne E!0 and 3*P eq E!0 and 21*Q eq E!0;
 a1,a2,a3,a4,a6:=Explode(aInvariants(E));r:=P[1];tt:=P[2];
 sl:=(3*r^2+2*a2*r+a4-a1*tt)/(2*tt+a1*r+a3);
 a1p:=a1+2*sl;a2p:=a2-sl*a1+3*r-sl^2;a3p:=a3+r*a1+2*tt;
 a4p:=a4-sl*a3+2*r*a2-(tt+r*sl)*a1+3*r^2-2*sl*tt;
 a6p:=a6+r*a4+r^2*a2+r^3-tt*a3-r*tt*a1-tt^2;
 assert a2p eq 0 and a4p eq 0 and a6p eq 0;
 xx:=Q[1]-r;yy:=Q[2]-tt-sl*xx;scale:=yy/xx^2;
 t:=-scale^2*xx;aa1:=scale*a1p;aa3:=scale^3*a3p;
 q:=(aa1-t-2)/(t-1);assert aa3 eq q*t^2+(1-q)*t;
 return [q,t];
end intrinsic;

intrinsic M3TwentyOneToCommon(v::SeqEnum,K::FldFin,F::FldFin) -> SeqEnum
{ Explicit residue-field map fixing the chosen quadratic base F }
 L:=Parent(v[1]);if L cmpeq K then return v;end if;
 d:=Degree(K,F);assert #L eq #K;
 f:=MinimalPolynomial(L.1,F);assert Degree(f) eq d;
 rr:=Roots(PolynomialRing(K)![K!a:a in Coefficients(f)]);assert #rr eq d;
 mp:=hom<L->K|rr[1][1]>;assert mp(L!F.1) eq K!F.1;
 return [mp(z):z in v];
end intrinsic;
intrinsic M3TwentyOneKey(v::SeqEnum,K::FldFin,F::FldFin,d::RngIntElt) -> SeqEnum
{ Canonical Frobenius-over-F key with an explicitly fixed base embedding }
 v:=M3TwentyOneToCommon(v,K,F);keys:=[];
 for i in [1..d] do
  ky:=&cat[[Integers()!a:a in Eltseq(z,PrimeField(K))]:z in v];
  Append(~keys,ky);v:=[z^(#F):z in v];
 end for;return Sort(keys)[1];
end intrinsic;
intrinsic M3TwentyOneOrbitKeys(X::MDCrvMod11,x::PlcCrvElt,K::FldFin) -> SeqEnum
{ All72 F-defined full3-basis and7-generator automorphisms }
 F:=BaseRing(X);d:=Degree(x);q,t:=Explode(M3TwentyOneToCommon([Evaluate(f,x):f in X`_coordinates],K,F));
 E:=EllipticCurve([q*t-q+t+2,K!0,q*t^2-q*t+t,K!0,K!0]);P:=E![0,0];Q:=E![-t,t^2];
 assert 3*P eq E!0 and 21*Q eq E!0;Q3:=7*Q;Q7:=3*Q;assert Q eq Q3+5*Q7;keys:=[];
 for mat in M3FifteenAutMatrices() do
  a,b,c,d1:=Explode(mat);PP:=a*P+b*Q3;QQ3:=c*P+d1*Q3;
  for u in [1,2,3] do
   QQ:=QQ3+5*u*Q7;
   Append(~keys,M3TwentyOneKey(M3TwentyOneNormalizedPair(E,PP,QQ),K,F,d));
  end for;
 end for;assert #keys eq 72;return keys;
end intrinsic;
intrinsic PlacesUpToFullAutM3TwentyOne(X::MDCrvMod11,S::SeqEnum) -> SeqEnum,SeqEnum,SeqEnum
{ Complete partition checked against all actual normalized-curve keys }
 require X`M eq 3 and Level(X) eq 21:"Require X1(3,21)";
 d:=Degree(S[1]);F:=BaseRing(X);K:=ext<F|d>;
 keys:=[M3TwentyOneKey([Evaluate(f,x):f in X`_coordinates],K,F,d):x in S];
 allowed:=Seqset(keys);assert #allowed eq #S;seen:=AssociativeArray();rs:=[];labels:=[];
 for i in [1..#S] do
  ky:=keys[i];
  if not IsDefined(seen,ky) then
   Append(~rs,S[i]);j:=#rs;
   for imagekey in M3TwentyOneOrbitKeys(X,S[i],K) do
    assert imagekey in allowed;
    if IsDefined(seen,imagekey) then assert seen[imagekey] eq j;else seen[imagekey]:=j;end if;
   end for;
  end if;
  Append(~labels,seen[ky]);
 end for;
 assert #Keys(seen) eq #S;return rs,labels,keys;
end intrinsic;
