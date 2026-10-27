/* CCRS arXiv1307.6174, Section4.8: exact new degree-eight verification. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(2*10^9);t0:=Cputime();
K<z>:=CyclotomicField(20);assert Degree(K) eq 8;i:=z^5;assert i^2 eq -1;
E:=EllipticCurve([K!1-i,-i,-i,K!0,K!0]);assert Discriminant(E) ne 0;
assert jInvariant(E) eq 1728;
G,mp:=TorsionSubgroup(E);assert Invariants(G) eq [5,10];
P:=mp(G.1);Q:=mp(G.2);assert 5*P eq E!0 and P ne E!0;
assert 10*Q eq E!0 and 5*Q ne E!0 and 2*Q ne E!0;
w:=WeilPairing(P,2*Q,5);assert w^5 eq 1 and w ne 1;
print "CM_5_10_EXPLICIT_POINTS",[[Eltseq(a):a in Eltseq(T)]:T in [P,Q]];
print "CM_5_10_PAIRING",Eltseq(w);
for p in [41,61] do
 F:=GF(p);R:=PolynomialRing(F);roots:=Roots(R!CyclotomicPolynomial(20));
 assert #roots eq 8;zz:=roots[1][1];ii:=zz^5;assert ii^2 eq -1;
 Ef:=EllipticCurve([1-ii,-ii,-ii,F!0,F!0]);assert Discriminant(Ef) ne 0;
 assert #Ef eq 50;
 print "CM_5_10_GOOD_REDUCTION",p,"zeta20",Integers()!zz,"order",#Ef;
end for;
print "CM_5_10_VERIFIED degree=8 invariants=[5,10] exact=true complete=true",Cputime(t0);quit;
