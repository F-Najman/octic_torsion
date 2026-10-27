intrinsic M4SixteenRelativeToCommon(v::SeqEnum,K::FldFin,F::FldFin) -> SeqEnum
{ Explicit residue-field map fixing the chosen quadratic base F }
 L:=Parent(v[1]);if L cmpeq K then return v;end if;
 d:=Degree(K,F);assert #L eq #K;
 f:=MinimalPolynomial(L.1,F);assert Degree(f) eq d;
 rr:=Roots(PolynomialRing(K)![K!a:a in Coefficients(f)]);assert #rr eq d;
 // A bare finite-field hom may conjugate F when L.1 has smaller
 // absolute degree. Evaluate the relative polynomial coefficients
 // explicitly, so the prescribed embedding of F is fixed.
 alpha:=rr[1][1];answer:=[];
 for z in v do
  coeff:=Eltseq(z,F);
  assert Evaluate(PolynomialRing(L)![L!a:a in coeff],L.1) eq z;
  Append(~answer,Evaluate(PolynomialRing(K)![K!a:a in coeff],alpha));
 end for;
 assert Evaluate(PolynomialRing(K)![K!a:a in Eltseq(L!F.1,F)],alpha) eq K!F.1;
 return answer;
end intrinsic;
intrinsic M4SixteenRelativeKey(v::SeqEnum,K::FldFin,F::FldFin,d::RngIntElt) -> SeqEnum
{ Canonical Frobenius-over-F key with an explicitly fixed base embedding }
 v:=M4SixteenRelativeToCommon(v,K,F);keys:=[];
 for i in [1..d] do
  ky:=&cat[[Integers()!a:a in Eltseq(z,PrimeField(K))]:z in v];
  Append(~keys,ky);v:=[z^(#F):z in v];
 end for;return Sort(keys)[1];
end intrinsic;
intrinsic M4SixteenRelativeOrbitKeys(X::MDCrvMod11,x::PlcCrvElt,K::FldFin) -> SeqEnum
{ All64 keys fixing the prescribed quadratic residue-field embedding }
 F:=BaseRing(X);dd:=Degree(x);
 q,t:=Explode(M4SixteenRelativeToCommon([Evaluate(f,x):f in X`_coordinates],K,F));
 B:=(q^2-1)*(t^2-1)/16;E:=EllipticCurve([K!1,B,B,K!0,K!0]);
 P:=E![0,0];Q:=E![(q+1)*(t^2-1)/8,(q+1)^2*(t-1)^2*(t+1)/32];
 assert 4*P eq E!0 and 2*P ne E!0 and 16*Q eq E!0 and 8*Q ne E!0;
 keys:=[];
 for par in M4SixteenAutParameters() do
  a,b,c,d:=Explode(par);PP:=a*P+4*b*Q;QQ:=c*P+d*Q;
  Append(~keys,M4SixteenRelativeKey(M4SixteenNormalizedPair(E,PP,QQ),K,F,dd));
 end for;assert #keys eq 64;return keys;
end intrinsic;
intrinsic PlacesUpToFullAutM4SixteenRelative(X::MDCrvMod11,S::SeqEnum) -> SeqEnum,SeqEnum,SeqEnum
{ Complete relative-base64-map partition, preserving the chosen base embedding }
 dd:=Degree(S[1]);F:=BaseRing(X);K:=ext<F|dd>;
 keys:=[M4SixteenRelativeKey([Evaluate(f,x):f in X`_coordinates],K,F,dd):x in S];
 allowed:=Seqset(keys);assert #allowed eq #S;seen:=AssociativeArray();rs:=[];labels:=[];
 for i in [1..#S] do
  ky:=keys[i];
  if not IsDefined(seen,ky) then
   Append(~rs,S[i]);j:=#rs;
   for im in M4SixteenRelativeOrbitKeys(X,S[i],K) do
    assert im in allowed;
    if IsDefined(seen,im) then assert seen[im] eq j;else seen[im]:=j;end if;
   end for;
  end if;
  Append(~labels,seen[ky]);
 end for;assert #Keys(seen) eq #S;return rs,labels,keys;
end intrinsic;
