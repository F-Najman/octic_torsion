/* Exact fixed-tower serialization, checked on every ideal generator. */
function C242SerializeFunction(a,K,A)
 a:=K!a;
 return [[[[Integers()!x:x in Coefficients(Numerator(z))],[Integers()!x:x in Coefficients(Denominator(z))]]:z in Eltseq(A!c)]:c in Eltseq(a)];
end function;
function C242DeserializeFunction(s,K,A)
 B:=BaseField(A);R:=PolynomialRing(ConstantField(A));
 return K![A![B!(R![Integers()!x:x in z[1]])/B!(R![Integers()!x:x in z[2]]):z in row]:row in s];
end function;
function C242SerializeDivisor(D,K,A)
 ps,ms:=Support(D);rows:=[* *];
 for i in [1..#ps] do
  aa,bb:=TwoGenerators(ps[i]);sa:=C242SerializeFunction(aa,K,A);sb:=C242SerializeFunction(bb,K,A);
  assert C242DeserializeFunction(sa,K,A) eq aa and C242DeserializeFunction(sb,K,A) eq bb;
  Append(~rows,[*Degree(ps[i]),IsFinite(ps[i]),ms[i],sa,sb*]);
 end for;
 return rows;
end function;
