/* Explicit independent lower bound and two finite-field upper bounds. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(10^9);t0:=Cputime();
K<z>:=CyclotomicField(20);assert Degree(K) eq 8;i:=z^5;
E:=EllipticCurve([1-i,-i,-i,K!0,K!0]);
assert Discriminant(E) eq 11-2*i and jInvariant(E) eq 1728;
P:=E![-1+z^2+z^3-z^4,-1+z^4-z^7];Q:=E![1+z^5,-2];
assert P ne E!0 and 5*P eq E!0;
assert 10*Q eq E!0 and 5*Q ne E!0 and 2*Q ne E!0;
assert WeilPairing(P,2*Q,5) eq -z^6;
assert (-z^6)^5 eq 1 and -z^6 ne 1;
for pr in [<41,2>,<61,8>] do
 p,zz:=Explode(pr);F:=GF(p);zz:=F!zz;
 assert Evaluate(PolynomialRing(F)!CyclotomicPolynomial(20),zz) eq 0;
 ii:=zz^5;Ef:=EllipticCurve([1-ii,-ii,-ii,F!0,F!0]);
 assert Discriminant(Ef) ne 0 and #Ef eq 50;
 print "CM_5_10_DIRECT_REDUCTION",p,Integers()!zz,#Ef;
end for;
print "CM_5_10_DIRECT_DONE subgroup=[5,10] upper_order=50 exact=true complete=true",Cputime(t0);quit;
