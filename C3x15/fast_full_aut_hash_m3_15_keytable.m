/* Integer key tables contain every actual native source place, exported
   in bounded fresh processes. This stage needs one common field only. */
intrinsic M3FifteenPairFromNativeKey(key::SeqEnum,K::FldFin,p::RngIntElt,d::RngIntElt) -> SeqEnum
{ Recover the q,t pair up to Frobenius from its invariant native key }
 assert #key eq 2*d+2 and #K eq p^d and key[1] in [0,1];
 RK:=PolynomialRing(K);f:=RK![K!a:a in key[2..d+2]];
 assert Degree(f) eq d and IsIrreducible(PolynomialRing(GF(p))!key[2..d+2]);
 roots:=Roots(f);assert #roots eq d;r:=roots[1][1];
 z:=Evaluate(RK![K!a:a in key[d+3..2*d+2]],r);
 pair:=key[1] eq 0 select [r,z] else [z,r];
 assert M3FifteenNativeKey(pair,p,d) eq key;return pair;
end intrinsic;

intrinsic PartitionFullAutM3FifteenKeyTable(inputkeys::SeqEnum,p::RngIntElt,d::RngIntElt) -> SeqEnum,SeqEnum
{ All48 automorphism orbits of a complete table, with exact coverage }
 allowed:=Seqset(inputkeys);assert #allowed eq #inputkeys;
 K:=GF(p^d);seen:=AssociativeArray();repkeys:=[];labels:=[];
 for i in [1..#inputkeys] do
  ky:=inputkeys[i];
  if not IsDefined(seen,ky) then
   Append(~repkeys,ky);j:=#repkeys;
   pair:=M3FifteenPairFromNativeKey(ky,K,p,d);
   for imagekey in M3FifteenOrbitKeysNative(pair,p,d) do
    assert imagekey in allowed;
    if IsDefined(seen,imagekey) then assert seen[imagekey] eq j;
    else seen[imagekey]:=j;end if;
   end for;
  end if;
  Append(~labels,seen[ky]);
 end for;
 assert #Keys(seen) eq #inputkeys;return repkeys,labels;
end intrinsic;
