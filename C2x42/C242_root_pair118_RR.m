/* Forget the sign of the independent3-point over X1(2,14). */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Y:=MDX11(2,14,GF(5));C:=Curve(Y);A,am:=AlgorithmicFunctionField(FunctionField(C));
B:=am(Y`_coordinates[1]);CC:=am(Y`_coordinates[2]);D:=B*(1-B-CC);
E14:=EllipticCurve([A|0,CC,0,D,0]);P2:=E14![0,0];Q14:=E14![B,B];
assert 14*Q14 eq E14!0 and 7*Q14 ne E14!0 and 2*Q14 ne E14!0 and 7*Q14 ne P2;
R<U>:=PolynomialRing(A);f:=DivisionPolynomial(E14,3);assert Degree(f) eq 4 and IsIrreducible(f);
K<x3>:=FunctionField(f);assert Genus(A) eq 4 and Genus(K) eq 31;
print "C242_FORGET3SIGN_GENUS31_MODEL",Genus(K),Cputime(t0);

load "C2x42/serialize_genus31_divisors.m";
data:=eval Read("C2x42/data/C242_g31_pair118_refinement_arithmetic.txt");
assert data[2] eq [* [*118,[*2,4,3,4,85*],13,true*],[*118,[*2,4,3,4,85*],17,false*] *];
assert #data[3] eq 1;DD:=DivisorGroup(K)!0;
for row in data[3][1] do
 d,isfin,m,sa,sb:=Explode(row);aa:=C242DeserializeFunction(sa,K,A);bb:=C242DeserializeFunction(sb,K,A);
 O:=isfin select MaximalOrderFinite(K) else MaximalOrderInfinite(K);I:=ideal<O|aa,bb>;fac:=Factorization(I);
 assert #fac eq 1 and fac[1][2] eq 1;q:=Place(fac[1][1]);assert Degree(q) eq d and IsFinite(q) eq isfin;
 DD+:=m*Divisor(q);
end for;
assert Degree(DD) eq 0 and not IsPrincipal(DD);
assert #Basis(DD:Reduction:=false,Simplification:="None") eq 0;
print "C242_ROOT_PAIR118_A17_NATIVE_RR_ZERO_PASS",Cputime(t0);quit;
