/* SL2(F3) x U8 / simultaneous -1, preserving the chosen pairing. */
intrinsic M3TwentyFourNormalizedPair(E::CrvEll,P::PtEll,Q::PtEll) -> SeqEnum
{ Canonical q,t for an independent order3 point and order24 point }
 assert P ne E!0 and Q ne E!0 and 3*P eq E!0;
 a1,a2,a3,a4,a6:=Explode(aInvariants(E));r:=P[1];tt:=P[2];
 sl:=(3*r^2+2*a2*r+a4-a1*tt)/(2*tt+a1*r+a3);
 a1p:=a1+2*sl;a2p:=a2-sl*a1+3*r-sl^2;a3p:=a3+r*a1+2*tt;
 a4p:=a4-sl*a3+2*r*a2-(tt+r*sl)*a1+3*r^2-2*sl*tt;
 a6p:=a6+r*a4+r^2*a2+r^3-tt*a3-r*tt*a1-tt^2;
 assert a2p eq 0 and a4p eq 0 and a6p eq 0;
 xx:=Q[1]-r;yy:=Q[2]-tt-sl*xx;scale:=yy/xx^2;
 t:=-scale^2*xx;aa1:=scale*a1p;aa3:=scale^3*a3p;
 q:=(aa1-t-2)/(t-1);assert aa3 eq q*t^2+(1-q)*t;return [q,t];
end intrinsic;
intrinsic M3TwentyFourOrbitKeys(X::MDCrvMod11,x::PlcCrvElt,K::FldFin) -> SeqEnum
{ All48 exact image keys, over the fixed prime residue base }
 require X`M eq 3 and Level(X) eq 24:"Require X1(3,24)";
 p:=#BaseRing(X);d:=Degree(x);assert IsPrime(p) and #K eq p^d;
 q,t:=Explode(M3FifteenToCommon([Evaluate(f,x):f in X`_coordinates],K));
 E:=EllipticCurve([q*t-q+t+2,K!0,q*t^2-q*t+t,K!0,K!0]);P:=E![0,0];Q:=E![-t,t^2];
 assert 3*P eq E!0 and 24*Q eq E!0 and 12*Q ne E!0 and 8*Q ne E!0;
 R3:=8*Q;R8:=3*Q;assert Q eq 2*R3+3*R8;keys:=[];
 for mat in M3FifteenAutMatrices() do
  a,b,c,d1:=Explode(mat);PP:=a*P+b*R3;RR3:=c*P+d1*R3;
  for u in [1,3] do
   QQ:=2*RR3+3*u*R8;Append(~keys,M3FifteenKey(M3TwentyFourNormalizedPair(E,PP,QQ),K,p,d));
  end for;
 end for;assert #keys eq 48;return keys;
end intrinsic;
intrinsic PlacesUpToFullAutM3TwentyFour(X::MDCrvMod11,S::SeqEnum) -> SeqEnum,SeqEnum,SeqEnum
{ Full checked48-map partition of every actual support place }
 d:=Degree(S[1]);p:=#BaseRing(X);assert IsPrime(p);K:=GF(p^d);
 keys:=[M3FifteenKey([Evaluate(f,x):f in X`_coordinates],K,p,d):x in S];
 allowed:=Seqset(keys);assert #allowed eq #S;seen:=AssociativeArray();rs:=[];labels:=[];
 for i in [1..#S] do
  ky:=keys[i];
  if not IsDefined(seen,ky) then
   Append(~rs,S[i]);j:=#rs;
   for im in M3TwentyFourOrbitKeys(X,S[i],K) do
    assert im in allowed;
    if IsDefined(seen,im) then assert seen[im] eq j;else seen[im]:=j;end if;
   end for;
  end if;
  Append(~labels,seen[ky]);
 end for;assert #Keys(seen) eq #S;return rs,labels,keys;
end intrinsic;
