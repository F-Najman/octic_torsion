/* Full pairing-preserving group SL2(Z/4) x U5 / simultaneous sign. */
intrinsic M4TwentyAutParameters() -> SeqEnum
{ All96 representatives, retaining the fixed fourth-root Weil pairing }
 return [[a,b,c,d,u]:a,b,c,d in [0..3],u in [1,2]|(a*d-b*c) mod 4 eq 1];
end intrinsic;
intrinsic M4TwentyOrbitKeys(X::MDCrvMod11,x::PlcCrvElt,K::FldFin) -> SeqEnum
{ All96 image keys with the prescribed base-field embedding fixed }
 require X`M eq 4 and Level(X) eq 20:"Require X1(4,20)";
 F:=BaseRing(X);dd:=Degree(x);q,t:=Explode(M4SixteenRelativeToCommon([Evaluate(f,x):f in X`_coordinates],K,F));
 B:=(q^2-1)*(t^2-1)/16;E:=EllipticCurve([K!1,B,B,K!0,K!0]);
 P:=E![0,0];Q:=E![(q+1)*(t^2-1)/8,(q+1)^2*(t-1)^2*(t+1)/32];
 assert 4*P eq E!0 and 2*P ne E!0 and 20*Q eq E!0 and 10*Q ne E!0 and 4*Q ne E!0;
 Q4:=5*Q;Q5:=4*Q;assert Q eq Q4-Q5;keys:=[];
 for par in M4TwentyAutParameters() do
  a,b,c,d,u:=Explode(par);PP:=a*P+b*Q4;QQ:=c*P+d*Q4-u*Q5;
  Append(~keys,M4SixteenRelativeKey(M4SixteenNormalizedPair(E,PP,QQ),K,F,dd));
 end for;assert #keys eq 96;return keys;
end intrinsic;
intrinsic PlacesUpToFullAutM4Twenty(X::MDCrvMod11,S::SeqEnum) -> SeqEnum,SeqEnum,SeqEnum
{ Complete checked96-map partition over the chosen finite base }
 dd:=Degree(S[1]);F:=BaseRing(X);K:=ext<F|dd>;
 keys:=[M4SixteenRelativeKey([Evaluate(f,x):f in X`_coordinates],K,F,dd):x in S];
 allowed:=Seqset(keys);assert #allowed eq #S;seen:=AssociativeArray();rs:=[];labels:=[];
 for i in [1..#S] do
  ky:=keys[i];
  if not IsDefined(seen,ky) then
   Append(~rs,S[i]);j:=#rs;
   for im in M4TwentyOrbitKeys(X,S[i],K) do
    assert im in allowed;
    if IsDefined(seen,im) then assert seen[im] eq j;else seen[im]:=j;end if;
   end for;
  end if;
  Append(~labels,seen[ky]);
 end for;assert #Keys(seen) eq #S;return rs,labels,keys;
end intrinsic;
