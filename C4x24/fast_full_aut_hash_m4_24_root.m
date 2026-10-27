/* Fresh full64 automorphisms of X1(4,24).
   P->aP+6bQ, Q->cP+dQ, with ad-6bc=1 mod4.
   b,c mod4 and d=1,5,7,11 represent U24/+-1.
   Rational moduli normalization is the standard validated M4 inverse. */
intrinsic M4TwentyFourRootNormalizedPair(E::CrvEll,P::PtEll,Q::PtEll) -> SeqEnum
{ Canonical q,t coordinates for an independent level(4,24) pair }
 a1,a2,a3,a4,a6:=Explode(aInvariants(E));r:=P[1];t0:=P[2];
 s:=(3*r^2+2*a2*r+a4-a1*t0)/(2*t0+a1*r+a3);
 a1p:=a1+2*s;a2p:=a2-s*a1+3*r-s^2;a3p:=a3+r*a1+2*t0;
 a4p:=a4-s*a3+2*r*a2-(t0+r*s)*a1+3*r^2-2*s*t0;
 a6p:=a6+r*a4+r^2*a2+r^3-t0*a3-r*t0*a1-t0^2;
 assert a4p eq 0 and a6p eq 0 and a3p eq a1p*a2p and a1p ne 0;
 scale:=1/a1p;b:=scale^2*a2p;
 xx:=scale^2*(Q[1]-r);yy:=scale^3*(Q[2]-t0-s*(Q[1]-r));
 q:=1+2*b/xx;t:=1+4*yy/((q+1)*xx);
 assert b eq (q^2-1)*(t^2-1)/16 and xx eq (q+1)*(t^2-1)/8;
 assert yy eq (q+1)^2*(t-1)^2*(t+1)/32;
 return [q,t];
end intrinsic;
intrinsic M4TwentyFourRootAutParameters() -> SeqEnum
{ Representatives modulo simultaneous negation for all64 automorphisms }
 return [[(InverseMod(d,4)*(1+6*b*c)) mod 4,b,c,d]:b,c in [0..3],d in [1,5,7,11]];
end intrinsic;
intrinsic M4TwentyFourRootOrbitKeys(X::MDCrvMod11,x::PlcCrvElt,K::FldFin) -> SeqEnum
{ Canonical keys of all64 automorphism images, over a prime base field }
 require X`M eq 4 and Level(X) eq 24:"Require X1(4,24)";
 p:=#BaseRing(X);dd:=Degree(x);assert IsPrime(p) and #K eq p^dd;
 q,t:=Explode(MDM2ToCommonField([Evaluate(f,x):f in X`_coordinates],K));
 B:=(q^2-1)*(t^2-1)/16;E:=EllipticCurve([K!1,B,B,K!0,K!0]);
 P:=E![0,0];Q:=E![(q+1)*(t^2-1)/8,(q+1)^2*(t-1)^2*(t+1)/32];
 assert 4*P eq E!0 and 2*P ne E!0 and 24*Q eq E!0 and 12*Q ne E!0 and 8*Q ne E!0;
 keys:=[];
 for par in M4TwentyFourRootAutParameters() do
  a,b,c,d:=Explode(par);PP:=a*P+6*b*Q;QQ:=c*P+d*Q;
  assert WeilPairing(PP,6*QQ,4) eq WeilPairing(P,6*Q,4);
  Append(~keys,MDM2FrobeniusKey(M4TwentyFourRootNormalizedPair(E,PP,QQ),K,p,dd));
 end for;assert #keys eq 64;return keys;
end intrinsic;
intrinsic PlacesUpToFullAutM4TwentyFourRoot(X::MDCrvMod11,S::SeqEnum) -> SeqEnum,SeqEnum,SeqEnum
{ Complete checked partition of actual closed places under the64 maps }
 dd:=Degree(S[1]);p:=#BaseRing(X);K:=Parent(Evaluate(X`_coordinates[1],S[1]));
 keys:=[MDM2FrobeniusKey([Evaluate(f,x):f in X`_coordinates],K,p,dd):x in S];
 allowed:=Seqset(keys);assert #allowed eq #S;seen:=AssociativeArray();rs:=[];labels:=[];
 for i in [1..#S] do
  ky:=keys[i];
  if not IsDefined(seen,ky) then
   Append(~rs,S[i]);j:=#rs;
   for im in M4TwentyFourRootOrbitKeys(X,S[i],K) do
    assert im in allowed;
    if IsDefined(seen,im) then assert seen[im] eq j;else seen[im]:=j;end if;
   end for;
  end if;
  Append(~labels,seen[ky]);
 end for;assert #Keys(seen) eq #S;return rs,labels,keys;
end intrinsic;
