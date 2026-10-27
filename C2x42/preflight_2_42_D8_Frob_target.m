/* Exact lightweight source counts and genus4 target setup. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/class_rational_reduction.m");rows:=[* *];
for tup in [[84,43],[28,15]] do
 J:=JH(tup[1],[tup[2]]);g:=Dimension(J);f:=FrobeniusPolynomial(J,5);assert Degree(f) eq 2*g;tr:=[];cl:=[];
 for n in [1..8] do
  a:=-n*Coefficient(f,2*g-n);for j in [1..n-1] do a-:=Coefficient(f,2*g-j)*tr[n-j];end for;
  Append(~tr,Integers()!a);num:=5^n+1-Integers()!a;
  for e in Divisors(n) do if e lt n then num-:=e*cl[e];end if;end for;
  assert num mod n eq 0;Append(~cl,num div n);
 end for;
 Append(~rows,[*tup,g,Coefficients(f),cl,Integers()!Evaluate(f,1)*]);
 print "M2N42_D8_FROB_PREFLIGHT",tup,g,cl,Evaluate(f,1),Cputime(t0);
end for;
Y:=MDX11(2,14,GF(5));C:=Curve(Y);cs:=Cusps(Y);c0:=Rep({c:c in cs|Degree(c) eq 1});
assert Genus(Y) eq 4;G,m1,m2:=ClassGroupRationalReduction(C,c0:Proof:=true);
hg:=[m2(Divisor(c)-Degree(c)*Divisor(c0)):c in cs];H:=sub<G|hg>;
print "M2N42_D8_TARGET_GROUP",Invariants(G),Invariants(H),[Degree(c):c in cs],Cputime(t0);
Append(~rows,[*Invariants(G),[Eltseq(r):r in Rows(RelationMatrix(G))],[Degree(c):c in cs],[Eltseq(a):a in hg],Invariants(H)*]);
Write("C2x42/data/preflight_2_42_D8_Frob_target_arithmetic.txt",Sprint(rows):Overwrite:=true);
print "M2N42_D8_FROB_TARGET_DONE complete=true",Cputime(t0);quit;
