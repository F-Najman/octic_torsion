/* Full F=Q(zeta3)-defined automorphisms of X1(3,15).
   Q3=5Q, Q5=3Q, Q=2Q3+2Q5. For M in SL2(F3) andu in U5,
   (P,Q3,Q5)->(aP+bQ3,cP+dQ3,uQ5). Modulo simultaneous -1,
   chooseu=1,2 and all24 matrices:48 automorphisms.
   Canonical(q,t) normal form is unique for an ordered level structure.
   This helper requires prime base field, avoiding any ambiguity about
   the embedding of the fixed Weil pairing into the residue field. */
intrinsic M3FifteenAutMatrices() -> SeqEnum
{ All24 matrices in SL2(F3), encoded by four integers }
 return [[a,b,c,d]:a,b,c,d in [0..2]|(a*d-b*c) mod 3 eq 1];
end intrinsic;

intrinsic M3FifteenNormalizedPair(E::CrvEll,P::PtEll,Q::PtEll) -> SeqEnum
{ Canonical q,t coordinates for an independent pair of orders3,15 }
 assert P ne E!0 and 3*P eq E!0 and 15*Q eq E!0;
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

intrinsic M3FifteenToCommon(v::SeqEnum,K::FldFin) -> SeqEnum
{ Explicit root-defined prime-field map between equal-size residue fields }
 L:=Parent(v[1]);if L cmpeq K then return v;end if;
 assert #L eq #K and Characteristic(L) eq Characteristic(K);
 if Degree(L) eq 1 then return [K!(Integers()!z):z in v];end if;
 fk:=PolynomialRing(K)![K!(Integers()!a):a in Coefficients(DefiningPolynomial(L,PrimeField(L)))];
 rr:=Roots(fk);assert #rr eq Degree(L);mp:=hom<L->K|rr[1][1]>;
 return [mp(z):z in v];
end intrinsic;

intrinsic M3FifteenKey(v::SeqEnum,K::FldFin,p::RngIntElt,d::RngIntElt) -> SeqEnum
{ Frobenius-canonical q,t key over the prime field }
 v:=M3FifteenToCommon(v,K);keys:=[];
 for i in [1..d] do
  ky:=&cat[[Integers()!a:a in Eltseq(z,PrimeField(K))]:z in v];
  assert #ky eq 2*d;Append(~keys,ky);v:=[z^p:z in v];
 end for;
 return Sort(keys)[1];
end intrinsic;

intrinsic M3FifteenOrbitKeys(X::MDCrvMod11,x::PlcCrvElt,K::FldFin) -> SeqEnum
{ All48 modular-automorphism keys of the closed point x }
 require X`M eq 3 and Level(X) eq 15 : "Require X1(3,15)";
 p:=#BaseRing(X);assert IsPrime(p);d:=Degree(x);assert #K eq p^d;
 q,t:=Explode(M3FifteenToCommon([Evaluate(f,x):f in X`_coordinates],K));
 E:=EllipticCurve([q*t-q+t+2,K!0,q*t^2-q*t+t,K!0,K!0]);
 P:=E![0,0];Q:=E![-t,t^2];assert 3*P eq E!0 and 15*Q eq E!0;
 Q3:=5*Q;Q5:=3*Q;assert Q eq 2*Q3+2*Q5;keys:=[];
 for mat in M3FifteenAutMatrices() do
  a,b,c,d1:=Explode(mat);PP:=a*P+b*Q3;QQ3:=c*P+d1*Q3;
  for u in [1,2] do
   QQ:=2*QQ3+2*u*Q5;
   Append(~keys,M3FifteenKey(M3FifteenNormalizedPair(E,PP,QQ),K,p,d));
  end for;
 end for;
 assert #keys eq 48;return keys;
end intrinsic;

intrinsic PlacesUpToFullAutM3FifteenHash(X::MDCrvMod11,S::SeqEnum[PlcCrvElt]) -> SeqEnum,SeqEnum
{ Exact complete orbit partition, checked against all source q,t keys }
 require X`M eq 3 and Level(X) eq 15 : "Require X1(3,15)";
 if #S eq 0 then return S,[];end if;
 d:=Degree(S[1]);assert &and[Degree(x) eq d:x in S];p:=#BaseRing(X);assert IsPrime(p);
 K:=GF(p^d);inputkeys:=[M3FifteenKey([Evaluate(f,x):f in X`_coordinates],K,p,d):x in S];
 allowed:=Seqset(inputkeys);assert #allowed eq #S;seen:=AssociativeArray();reps:=[];labels:=[];
 for i in [1..#S] do
  ky:=inputkeys[i];
  if not IsDefined(seen,ky) then
   Append(~reps,S[i]);j:=#reps;
   for imagekey in M3FifteenOrbitKeys(X,S[i],K) do
    assert imagekey in allowed;
    if IsDefined(seen,imagekey) then assert seen[imagekey] eq j;
    else seen[imagekey]:=j;end if;
   end for;
  end if;
  Append(~labels,seen[ky]);
 end for;
 assert #Keys(seen) eq #S;return reps,labels;
end intrinsic;
