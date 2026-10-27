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

if not assigned Shard then Shard:=0;end if;
if Type(Shard) eq MonStgElt then Shard:=StringToInteger(Shard);end if;assert Shard in [0,1];
load "C2x42/serialize_genus31_divisors.m";
data:=eval Read(Sprintf("C2x42/data/C242_g31_remaining_shard%o_arithmetic.txt",Shard));
assert data[3] and #data[2] eq 13 and #data[4] eq 13;
for i in [1..#data[4]] do
 idx,desc,ell,rows:=Explode(data[4][i]);assert data[2][i][1] eq idx and data[2][i][2] eq desc and data[2][i][3] eq ell;
 DD:=DivisorGroup(K)!0;
 for row in rows do
  d,isfin,m,sa,sb:=Explode(row);aa:=C242DeserializeFunction(sa,K,A);bb:=C242DeserializeFunction(sb,K,A);
  O:=isfin select MaximalOrderFinite(K) else MaximalOrderInfinite(K);I:=ideal<O|aa,bb>;fac:=Factorization(I);
  assert #fac eq 1 and fac[1][2] eq 1;q:=Place(fac[1][1]);assert Degree(q) eq d and IsFinite(q) eq isfin;
  DD+:=m*Divisor(q);
 end for;
 assert Degree(DD) eq 0;yes:=IsPrincipal(DD);assert yes eq data[2][i][4];
 if not yes then assert #Basis(DD:Reduction:=false,Simplification:="None") eq 0;end if;
 print "C242_G31_SERIALIZED_RR_REPLAY",Shard,idx,ell,yes,Cputime(t0);
end for;
print "C242_G31_SERIALIZED_RR_REPLAY_DONE",Shard,#data[4],"complete=true",Cputime(t0);quit;
