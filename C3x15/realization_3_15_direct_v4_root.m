/* Standalone exact octic example. No modular-curve or RR dependency. */
SetColumns(0);SetQuitOnError(true);tt:=Cputime();
F<z>:=CyclotomicField(3);R<t>:=PolynomialRing(F);f:=t^4-5*t^3+8*t^2-6*t+3;
assert IsIrreducible(f);K<w>:=ext<F|f>;assert Degree(K) eq 4 and Degree(F) eq 2;
a1:=((-7*z-8)*w^3+(23*z+34)*w^2+(-8*z-43)*w)/6+(-3*z+10)/2;
a3:=((17*z+13)*w^3+(-82*z-71)*w^2+(106*z+119)*w)/12-2*z-5;
E:=EllipticCurve([a1,K!0,a3,K!0,K!0]);P:=E![0,0];
xq:=((-z+1)*w^3+(8*z-2)*w^2+(-17*z-1)*w)/6+3*z/2;
yq:=((z+2)*w^3+(-4*z-8)*w^2+(3*z+6)*w+2*z+3)/4;
Q:=E![xq,yq];assert Discriminant(E) ne 0 and 3*P eq E!0 and P ne E!0;
assert 15*Q eq E!0 and 5*Q ne E!0 and 3*Q ne E!0;assert WeilPairing(P,5*Q,3) eq K!z^2;
tau:=-yq^2/xq^3;assert yq eq xq^2 and tau eq -xq;
assert Rank(Matrix(F,[Eltseq(tau^i):i in [0..3]])) eq 4;
print "M3N15_FINE_MODULI_FIELD_RELATIVE_DEGREE",4;
assert jInvariant(E) eq K!(1331/8);assert Discriminant(f) eq -507;
print "M3N15_DIRECT_EXACT_OCTIC_FIELD_AND_LEVEL",true,Cputime(tt);
function EvAt43(a)
 return &+[(&+[(GF(43)!Eltseq(Eltseq(K!a)[i])[j])*(GF(43)!6)^(j-1):j in [1..#Eltseq(Eltseq(K!a)[i])]])*(GF(43)!4)^(i-1):i in [1..#Eltseq(K!a)]];
end function;
k:=GF(43);assert k!6^2+6+1 eq 0 and (k!4)^4-5*(k!4)^3+8*(k!4)^2-6*(k!4)+3 eq 0;
assert k!Discriminant(f) ne 0 and k!3 ne 0;
Er:=EllipticCurve([EvAt43(a):a in aInvariants(E)]);assert Discriminant(Er) ne 0 and #Er eq 45;
assert Invariants(AbelianGroup(Er)) eq [3,15];
print "M3N15_DIRECT_GOOD_UNRAMIFIED_PRIME43",aInvariants(Er),Discriminant(Er),#Er,Invariants(AbelianGroup(Er)),Cputime(tt);
Write("C3x15/data/realization_3_15_direct_v4_root_arithmetic.txt",Sprint([*[*Eltseq(c):c in Coefficients(f)*],[*[*Eltseq(c):c in Eltseq(a)*]:a in aInvariants(E)*],[*[*Eltseq(c):c in Eltseq(a)*]:a in [P[1],P[2],Q[1],Q[2]]*],[Integers()!a:a in aInvariants(Er)],Integers()!Discriminant(Er),#Er*]):Overwrite:=true);
print "M3N15_DIRECT_FULL_TORSION_3_15_CERTIFIED",Cputime(tt);quit;
