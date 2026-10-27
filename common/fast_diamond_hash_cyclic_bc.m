/* Canonical Tate(b,c) hash for cyclic level structures, N>=10.
   Attach deps/fast_diamond_hash_m2.m first for explicit root-defined
   finite-field conversions; no generic Embed is called. */
intrinsic MDCyclicBCHashKeys(X::MDCrvMod1,x::PlcCrvElt,K::FldFin) -> SeqEnum
{ Frobenius keys of all cyclic diamond translates, using canonical Tate b,c }
 p:=#BaseRing(X);d:=Degree(x);assert #K eq p^d;
 bc:=MDM2ToCommonField([Evaluate(f,x):f in X`_coordinates],K);
 b,c:=Explode(bc);E:=EllipticCurve([1-c,-b,-b,K!0,K!0]);P:=E![0,0];
 keys:=[];n:=Level(X);
 for a in [1..n div 2] do
  if GCD(a,n) eq 1 then
   Append(~keys,MDM2FrobeniusKey(TateNormalForm_bc(E,a*P),K,p,d));
  end if;
 end for;
 assert #keys eq EulerPhi(n) div 2;return keys;
end intrinsic;

intrinsic PlacesUpToDiamondCyclicHashBC(X::MDCrvMod1,S::SeqEnum[PlcCrvElt]) -> SeqEnum[PlcCrvElt]
{ Complete cyclic diamond partition by unique Tate b,c and Frobenius keys }
 if #S eq 0 then return S;end if;
 p:=#BaseRing(X);assert IsPrime(p);d:=Degree(S[1]);
 assert &and[Degree(x) eq d:x in S];
 K:=Parent(Evaluate(X`_coordinates[1],S[1]));assert #K eq p^d;
 seen:=AssociativeArray();reps:=[];
 for x in S do
  bc:=[Evaluate(f,x):f in X`_coordinates];key:=MDM2FrobeniusKey(bc,K,p,d);
  if IsDefined(seen,key) then continue;end if;
  Append(~reps,x);
  for ky in MDCyclicBCHashKeys(X,x,K) do
   if IsDefined(seen,ky) then assert seen[ky] eq #reps;
   else seen[ky]:=#reps;end if;
  end for;
  assert IsDefined(seen,key) and seen[key] eq #reps;
 end for;
 assert #Keys(seen) eq #S;return reps;
end intrinsic;
