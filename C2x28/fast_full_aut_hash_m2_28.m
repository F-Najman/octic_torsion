/* The24 Q-defined level automorphisms on X1(2,28):
   P2 -> P2+14bQ28, Q28 -> cP2+dQ28,
   b,c mod2, d in U28/+-1. Hecke at every odd prime commutes.
   For a1=a3=0 translate newP2=(r,0), then canonically rescale at newQ. */
intrinsic M2TwentyEightAutParameters() -> SeqEnum
{ All24 representatives modulo simultaneous negation }
 return [[b,c,d]:b,c in [0,1],d in [1,3,5,9,11,13]];
end intrinsic;
intrinsic M2TwentyEightNormalizedPair(E::CrvEll,P::PtEll,Q::PtEll) -> SeqEnum
{ Canonical b,c coordinates after a2-preserving short-model translation }
 a1,a2,a3,a4,a6:=Explode(aInvariants(E));assert a1 eq 0 and a3 eq 0 and a6 eq 0;
 assert 2*P eq E!0 and P ne E!0 and P[2] eq 0;
 r:=P[1];xx:=Q[1]-r;yy:=Q[2];assert xx ne 0 and yy ne 0;
 bn:=xx^3/yy^2;cn:=(a2+3*r)*xx^2/yy^2;
 assert (a4+2*a2*r+3*r^2)*xx^4/yy^4 eq bn*(1-bn-cn);
 return [bn,cn];
end intrinsic;
intrinsic M2TwentyEightOrbitKeys(X::MDCrvMod11,x::PlcCrvElt,K::FldFin) -> SeqEnum
{ All24 keys for the exact Q-defined level automorphism group }
 require X`M eq 2 and Level(X) eq 28:"Require X1(2,28)";
 p:=#BaseRing(X);dd:=Degree(x);assert IsPrime(p) and #K eq p^dd;
 b,cc:=Explode(MDM2ToCommonField([Evaluate(f,x):f in X`_coordinates],K));
 E:=EllipticCurve([K!0,cc,K!0,b*(1-b-cc),K!0]);P:=E![0,0];Q:=E![b,b];Q2:=14*Q;
 assert 2*P eq E!0 and Q2 ne P and Q2 ne E!0 and 2*Q2 eq E!0;
 keys:=[];
 for par in M2TwentyEightAutParameters() do
  b1,c1,a:=Explode(par);PP:=P+b1*Q2;QQ:=c1*P+a*Q;
  Append(~keys,MDM2FrobeniusKey(M2TwentyEightNormalizedPair(E,PP,QQ),K,p,dd));
 end for;assert #keys eq 24;return keys;
end intrinsic;
intrinsic PlacesUpToFullAutM2TwentyEight(X::MDCrvMod11,S::SeqEnum) -> SeqEnum,SeqEnum,SeqEnum
{ Checked complete partition under all24 actual level transformations }
 dd:=Degree(S[1]);p:=#BaseRing(X);K:=Parent(Evaluate(X`_coordinates[1],S[1]));
 keys:=[MDM2FrobeniusKey([Evaluate(f,x):f in X`_coordinates],K,p,dd):x in S];
 allowed:=Seqset(keys);assert #allowed eq #S;seen:=AssociativeArray();rs:=[];labels:=[];
 for i in [1..#S] do
  ky:=keys[i];
  if not IsDefined(seen,ky) then
   Append(~rs,S[i]);j:=#rs;
   for im in M2TwentyEightOrbitKeys(X,S[i],K) do
    assert im in allowed;
    if IsDefined(seen,im) then assert seen[im] eq j;else seen[im]:=j;end if;
   end for;
  end if;
  Append(~labels,seen[ky]);
 end for;assert #Keys(seen) eq #S;return rs,labels,keys;
end intrinsic;
