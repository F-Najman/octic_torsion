/* Explicit finite-field isomorphism on coefficients; creates no hom object.
   Different roots change only the Frobenius embedding and hence not keys. */
function C48CommonFieldByCoefficients(values,K)
 L:=Parent(values[1]);assert &and[Parent(v) cmpeq L:v in values];
 if L cmpeq K then return values;end if;
 assert #L eq #K and Characteristic(L) eq Characteristic(K);
 p:=Characteristic(K);d:=Degree(L);assert d eq Degree(K);
 f:=DefiningPolynomial(L,PrimeField(L));assert Degree(f) eq d;
 fk:=PolynomialRing(K)![K!(Integers()!a):a in Coefficients(f)];
 rr:=Roots(fk);assert #rr eq d;r:=rr[1][1];assert Evaluate(fk,r) eq 0;
 result:=[&+[K!(Integers()!a[i])*r^(i-1):i in [1..#a]] where a is Eltseq(v,PrimeField(L)):v in values];
 assert &and[Degree(MinimalPolynomial(result[i],PrimeField(K))) eq Degree(MinimalPolynomial(values[i],PrimeField(L))):i in [1..#values]];
 return result;
end function;
