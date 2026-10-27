/* Embedding-free keys. For d=2,3,4 all proper finite subfields are
   nested. A generating pair therefore has a primitive first or second
   coordinate. Its minimum polynomial plus the other coordinate in
   its power basis determines exactly its Frobenius orbit. */
intrinsic M3FifteenNativeKey(v::SeqEnum,p::RngIntElt,d::RngIntElt) -> SeqEnum
{ Frobenius-orbit key without constructing residue-field embeddings }
 require #v eq 2 and d in [2,3,4]:"Require pair in degree2,3 or4";
 L:=Parent(v[1]);Fp:=PrimeField(L);assert #Fp eq p and #L eq p^d;
 q,t:=Explode(v);fq:=MinimalPolynomial(q,Fp);
 if Degree(fq) eq d then r:=q;z:=t;flag:=0;fr:=fq;
 else r:=t;z:=q;flag:=1;fr:=MinimalPolynomial(t,Fp);end if;
 assert Degree(fr) eq d;
 M:=Matrix(Fp,d,d,&cat[Eltseq(r^i,Fp):i in [0..d-1]]);
 coeffs:=Eltseq(Solution(M,Vector(Fp,Eltseq(z,Fp))));
 assert &+[coeffs[i+1]*r^i:i in [0..d-1]] eq z;
 return [flag] cat [Integers()!a:a in Coefficients(fr)] cat [Integers()!a:a in coeffs];
end intrinsic;

intrinsic M3FifteenOrbitKeysNative(v::SeqEnum,p::RngIntElt,d::RngIntElt) -> SeqEnum
{ Native keys of all48 modular automorphisms }
 q,t:=Explode(v);L:=Parent(q);
 E:=EllipticCurve([q*t-q+t+2,L!0,q*t^2-q*t+t,L!0,L!0]);
 P:=E![0,0];Q:=E![-t,t^2];assert 3*P eq E!0 and 15*Q eq E!0;
 Q3:=5*Q;Q5:=3*Q;assert Q eq 2*Q3+2*Q5;keys:=[];
 for mat in M3FifteenAutMatrices() do
  a,b,c,d1:=Explode(mat);PP:=a*P+b*Q3;QQ3:=c*P+d1*Q3;
  for u in [1,2] do
   QQ:=2*QQ3+2*u*Q5;
   Append(~keys,M3FifteenNativeKey(M3FifteenNormalizedPair(E,PP,QQ),p,d));
  end for;
 end for;
 assert #keys eq 48;return keys;
end intrinsic;

intrinsic PlacesUpToFullAutM3FifteenNative(coordinates::SeqEnum,S::SeqEnum,p::RngIntElt) -> SeqEnum,SeqEnum
{ Exact native-place partition with embedding-free keys and full coverage }
 if #S eq 0 then return S,[];end if;
 d:=Degree(S[1]);assert &and[Degree(x) eq d:x in S];
 inputkeys:=[M3FifteenNativeKey([Evaluate(f,x):f in coordinates],p,d):x in S];
 allowed:=Seqset(inputkeys);assert #allowed eq #S;
 seen:=AssociativeArray();reps:=[];labels:=[];
 for i in [1..#S] do
  ky:=inputkeys[i];
  if not IsDefined(seen,ky) then
   Append(~reps,S[i]);j:=#reps;
   for imagekey in M3FifteenOrbitKeysNative([Evaluate(f,S[i]):f in coordinates],p,d) do
    assert imagekey in allowed;
    if IsDefined(seen,imagekey) then assert seen[imagekey] eq j;
    else seen[imagekey]:=j;end if;
   end for;
  end if;
  Append(~labels,seen[ky]);
 end for;
 assert #Keys(seen) eq #S;return reps,labels;
end intrinsic;
