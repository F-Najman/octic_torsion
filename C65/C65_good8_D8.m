// Fresh octic pure NC8 necessary test; no cusp complement.
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);
AttachSpec("mdmagma/v2/mdmagma.spec");t0:=Cputime();
X13:=MDX1(13,GF(2));C13:=Curve(X13);F13:=FunctionField(C13);
EE:=EllipticCurve(X13`_E);PP:=EE!X13`_P;
xy3:=TateNormalForm_xy(EE,3*PP);xy9:=TateNormalForm_xy(EE,9*PP);
phi3:=iso<C13->C13|[xy3[1],xy3[2],F13!1],[xy9[1],xy9[2],F13!1]:Check:=false>;
GG:=AutomorphismGroup(C13,[phi3]);assert #GG eq 3;
Q,pi:=CurveQuotient(GG);assert Genus(Q) eq 0 and Degree(pi) eq 3;
FQ:=FunctionField(Q);zfun:=Pullback(pi,FQ.2/FQ.1)+1;
function J13(z)
 return (z^16+z^4+1)/(z^28+z^25+z^24+z^21+z^20+z^17+z^16+z^13);
end function;
assert J13(zfun) eq jInvariant(EE);
function Pair65(E,P5,P13)
 b,c:=Explode(TateNormalForm_bc(E,P5));assert b eq c;
 K:=BaseRing(E);xy:=TateNormalForm_xy(E,P13);pt:=C13(K)![xy[1],xy[2],K!1];
 z:=Evaluate(zfun,pt);assert J13(z) eq jInvariant(E);
 return z,b;
end function;
F<a>:=GF(256);pairset:={};
for j in F do
 E0:=EllipticCurveFromjInvariant(j);
 for E in Twists(E0) do
  if #E mod 65 ne 0 then continue;end if;
  assert #E eq 260;
  A,mp:=AbelianGroup(E);assert Invariants(A) eq [260];P:=4*mp(A.1);
  for k in [k:k in [1..64]|GCD(k,65) eq 1] do
   R5:=13*k*P;R13:=5*k*P;
   z,t5:=Pair65(E,R5,R13);
   z3,t3:=Pair65(E,3*R5,3*R13);assert z3 eq z and t3 eq 1/t5;
   z6,t6:=Pair65(E,6*R5,6*R13);assert z6 eq z+1 and t6 eq t5;
   Include(~pairset,<z,t5>);
  end for;
 end for;
end for;
assert #pairset eq 64;
function Key(z,t5)
 return Minimum([Eltseq(z^(2^i),GF(2)) cat Eltseq(t5^(2^i),GF(2)):i in [0..7]]);
end function;
lookup:=AssociativeArray();for r in pairset do lookup[Key(r[1],r[2])]:=r;end for;
assert #Keys(lookup) eq 8;left:=Keys(lookup);reps:=[];
while #left gt 0 do
 k:=Minimum(Setseq(left));r:=lookup[k];Append(~reps,r);
 orb:={Key(r[1],r[2]),Key(r[1]+1,r[2]),Key(r[1],1/r[2]),Key(r[1]+1,1/r[2])};
 assert #orb eq 4 and orb subset Keys(lookup);left diff:=orb;
end while;
assert #reps eq 2;
printf "OCTIC65_GOOD8_EXHAUSTIVE pairs=%o closed_places=%o diamond_orbits=%o time=%o\n",#pairset,#Keys(lookup),#reps,Cputime(t0);
print "OCTIC65_GOOD8_SOURCE_DIAMONDS_VERIFIED a3=t_inverse a6=z_plus1";
for i in [1..#reps] do
 printf "OCTIC65_GOOD8_REP index=%o fieldpoly=%o z=%o t=%o\n",i,Eltseq(DefiningPolynomial(F,GF(2))),Eltseq(reps[i][1],GF(2)),Eltseq(reps[i][2],GF(2));
end for;
F0<z>:=FunctionField(GF(2));j13:=J13(z);R<t>:=PolynomialRing(F0);
f:=t^12+t^8+t^4+1+j13*(t^7+t^6+t^5);assert IsIrreducible(f);
AF<t5>:=FunctionField(f);assert Genus(AF) eq 21;
O:=MaximalOrderFinite(AF);den:=AF!Denominator(j13);
printf "OCTIC65_GOOD8_MODEL genus=21 degree=12 time=%o\n",Cputime(t0);
function PairPlace(z0,t0)
 K:=Parent(z0);mp:=MinimalPolynomial(z0,GF(2));d:=Degree(mp);
 MB:=Matrix(GF(2),[Eltseq(z0^i,GF(2)):i in [0..d-1]]);assert Rank(MB) eq d;
 xv:=Vector(GF(2),Eltseq(t0,GF(2)));assert Rank(VerticalJoin(MB,Matrix(xv))) eq d;
 co:=Solution(MB,xv);g:=&+[AF!co[i+1]*(AF!z)^i:i in [0..d-1]];
 I:=ideal<O|O!Evaluate(mp,AF!z),O!(den*(t5-g))>;assert IsPrime(I);
 P:=Place(I);assert Degree(P) eq d;return P;
end function;
// Reconstruct a representative pair by the original exhaustive curve scan.
for index in [1..#reps] do
 r:=reps[index];found:=false;
 for E in Twists(EllipticCurveFromjInvariant(J13(r[1]))) do
  if #E ne 260 then continue;end if;
  A,mp:=AbelianGroup(E);P65:=4*mp(A.1);
  for k in [k:k in [1..64]|GCD(k,65) eq 1] do
   P5:=13*k*P65;P13:=5*k*P65;zz,tt:=Pair65(E,P5,P13);
   if zz eq r[1] and tt eq r[2] then Erep:=E;R5rep:=P5;R13rep:=P13;found:=true;break;end if;
  end for;
  if found then break;end if;
 end for;
 assert found;E:=Erep;P5:=R5rep;P13:=R13rep;
 P:=PairPlace(r[1],r[2]);PT:=PairPlace(r[1],1/r[2]);
 isogs:=MDIsogenies(E,3);assert #isogs eq 4;out:=[];outT:=[];
 for phi in isogs do
  Ephi:=Codomain(phi);R5phi:=phi(P5);R13phi:=phi(P13);
  zz,tt:=Pair65(Ephi,R5phi,R13phi);
  z3,t3:=Pair65(Ephi,3*R5phi,3*R13phi);assert z3 eq zz and t3 eq 1/tt;
  Append(~out,PairPlace(zz,tt));Append(~outT,PairPlace(zz,1/tt));
 end for;
 T:=DivisorGroup(AF)!0;TT:=DivisorGroup(AF)!0;
 for Q in Setseq(Seqset(out)) do
  mult:=#[q:q in out|q eq Q];assert (8*mult) mod Degree(Q) eq 0;T+:=((8*mult) div Degree(Q))*Divisor(Q);
 end for;
 for Q in Setseq(Seqset(outT)) do
  mult:=#[q:q in outT|q eq Q];assert (8*mult) mod Degree(Q) eq 0;TT+:=((8*mult) div Degree(Q))*Divisor(Q);
 end for;
 assert Degree(T) eq 32 and Degree(TT) eq 32;
 B:=TT-T+2*Divisor(PT)-2*Divisor(P);assert Degree(B) eq 0;
 printf "OCTIC65_GOOD8_TEST_START index=%o image_degrees=%o time=%o\n",index,[Degree(Q):Q in Seqset(out)],Cputime(t0);
 prin:=IsPrincipal(B);
 printf "OCTIC65_GOOD8_A3_DIFFERENCE index=%o principal=%o time=%o\n",index,prin,Cputime(t0);
end for;
print "OCTIC65_GOOD8_A3_DONE";quit;
