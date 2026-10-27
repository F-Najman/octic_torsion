SetColumns(0);SetQuitOnError(true);SetMemoryLimit(24*10^9);tt:=Cputime();
CQ:=SmallModularCurve(11);jraw:=jInvariant(CQ,11);
CC:=ChangeRing(CQ,GF(2));assert IsNonsingular(CC) and Genus(CC) eq 1;
RC:=CoordinateRing(Ambient(CC));jn:=RC!Numerator(jraw);jd:=RC!Denominator(jraw);
K,kmap:=AlgorithmicFunctionField(FunctionField(CC));jK:=kmap(FunctionField(CC)!(jn/jd));
// The classical X0(7) parameter has j=(v^2+v+1)^4/v^7 in characteristic2.
Fv<v>:=FunctionField(GF(2));j7:=(v^2+v+1)^4/v^7;j7dual:=(v^2+v+1)^4/v;
Phi7:=ClassicalModularPolynomial(7);
assert Evaluate(Phi7,[j7,j7dual]) eq 0;
R<x>:=PolynomialRing(K);f:=x^8+jK*x^7+x^4+1;
assert IsIrreducible(f);AF<u>:=FunctionField(f);assert Genus(AF) eq 7;
jaf:=AF!jK;
cs:=Poles(jaf);widths:=[-Valuation(jaf,c):c in cs];
assert #cs eq 4 and &and[Degree(c) eq 1:c in cs] and Seqset(widths) eq {1,7,11,77};
printf "OCTIC77_X0_MODEL genus7=true modular_polynomial_verified=true cusp_widths=%o time=%o\n",widths,Cputime(tt);
DG:=DivisorGroup(AF);Dzero:=DG!0;G,gm:=RayClassGroup(Dzero);
inv:=Invariants(G);assert &*[a:a in inv|a ne 0] eq 1800 and #[a:a in inv|a eq 0] eq 1;
printf "OCTIC77_PICARD invariants=%o time=%o\n",inv,Cputime(tt);
C6:=AbelianGroup([6]);chars:=[];
cols:=[];
for i in [1..4] do
 Append(~cols,Eltseq(Divisor(cs[i])@@gm));
 printf "OCTIC77_CUSP_CLASS index=%o class=%o time=%o\n",i,cols[i],Cputime(tt);
end for;
expect:=[w eq 77 select 0 else (w eq 11 select 2 else (w eq 7 select 3 else 5)):w in widths];
for r in CartesianPower({0..5},Ngens(G)) do
 rr:=[r[i]:i in [1..Ngens(G)]];
 if not &and[inv[i]*rr[i] mod 6 eq 0:i in [1..#inv]] then continue;end if;
 if &and[(&+[cols[i][j]*rr[j]:j in [1..#inv]]) mod 6 eq expect[i]:i in [1..4]] then Append(~chars,rr);end if;
end for;
printf "OCTIC77_C6_CUSP_CHARACTERS count=%o images=%o time=%o\n",#chars,chars,Cputime(tt);
assert #chars eq 6;
wraw:=DefiningPolynomials(AtkinLehnerInvolution(CQ,11,11));ww:=[RC!e:e in wraw];
jndual:=Evaluate(jn,ww);jddual:=Evaluate(jd,ww);
O:=MaximalOrderFinite(AF);OK:=MaximalOrderFinite(K);
numY,denY:=IntegralSplit(K.1,OK);numJ,denJ:=IntegralSplit(jK,OK);
IX:=AF!BaseField(K).1;IY:=AF!K!numY;IT:=(AF!K!denJ)*AF.1;
assert &and[v in O:v in [IX,IY,IT]];
fx:=(K!BaseField(K).1)@@kmap;fy:=(K!numY)@@kmap;fden:=(K!denJ)@@kmap;

function CoordPlace(r)
 L:=Parent(r[1]);degrees:=[Degree(MinimalPolynomial(s,GF(2))):s in r];d:=LCM(degrees);
 candidates:=[<IX,r[1]>,<IY,r[2]>,<IT,r[3]>,<IX+IY,r[1]+r[2]>,<IX+IT,r[1]+r[3]>,<IY+IT,r[2]+r[3]>,<IX+IY+IT,r[1]+r[2]+r[3]>];
 found:=false;
 for q in candidates do
  mp:=MinimalPolynomial(q[2],GF(2));if Degree(mp) eq d then gen:=q[1];val:=q[2];found:=true;break;end if;
 end for;
 assert found;
 MB:=Matrix(GF(2),[Eltseq(val^i,GF(2)):i in [0..d-1]]);assert Rank(MB) eq d;
 gens:=[O!Evaluate(mp,gen)];
 for i in [1..3] do
  xv:=Vector(GF(2),Eltseq(r[i],GF(2)));assert Rank(VerticalJoin(MB,Matrix(xv))) eq d;
  co:=Solution(MB,xv);poly:=&+[AF!co[k+1]*gen^k:k in [0..d-1]];
  Append(~gens,O!([IX,IY,IT][i]-poly));
 end for;
 II:=ideal<O|gens>;assert IsPrime(II);P:=Place(II);assert Degree(P) eq d;return P;
end function;

tested:=0;
for e in [10,11] do
 F<a>:=GF(2^e);lookup:=AssociativeArray();CCF:=ChangeRing(CC,F);
 for p0 in Points(CCF) do
  p:=CC(F)!Eltseq(p0);co:=Eltseq(p);
  if Evaluate(jd,co) eq 0 or Evaluate(jddual,co) eq 0 then continue;end if;
  jj:=Evaluate(jn,co)/Evaluate(jd,co);jj11:=Evaluate(jndual,co)/Evaluate(jddual,co);
  ky:=<jj,jj11>;
  if not IsDefined(lookup,ky) then lookup[ky]:=[];end if;
  Append(~lookup[ky],p);
 end for;
 printf "OCTIC77_X011_FINITE_LOOKUP exponent=%o pairs=%o time=%o\n",e,#Keys(lookup),Cputime(tt);
 PR<X>:=PolynomialRing(F);seen:={};used:=0;
 for jj in F do
  jkey:=Minimum([Eltseq(jj^(2^i),GF(2)):i in [0..e-1]]);
  if jkey in seen then continue;end if;
  E0:=EllipticCurveFromjInvariant(jj);
  for E in Twists(E0) do
   if #E mod 77 ne 0 then continue;end if;
   A,mp:=AbelianGroup(E);
   if Invariants(A)[#Invariants(A)] mod 77 ne 0 then continue;end if;
   assert #Invariants(A) eq 1;
   PN:=(#E div 77)*mp(A.1);assert Order(PN) eq 77;
   R7:=11*PN;R11:=7*PN;
   ker7:=&*[X-(k*R7)[1]:k in [1..3]];ker11:=&*[X-(k*R11)[1]:k in [1..5]];
   E7,ph7:=IsogenyFromKernel(E,ker7);E11,ph11:=IsogenyFromKernel(E,ker11);
   assert Order(ph7(R11)) eq 11 and Order(ph11(R7)) eq 7;
   ky:=<jInvariant(E),jInvariant(E11)>;
   if not IsDefined(lookup,ky) or #lookup[ky] ne 1 then continue;end if;
   rs:=[r[1]:r in Roots(X^8+jInvariant(E)*X^7+X^4+1)|r[1] ne 0 and (r[1]^8+r[1]^4+1)/r[1] eq jInvariant(E7)];
   if #rs ne 1 then continue;end if;
   p:=lookup[ky][1];v7:=rs[1];r:=<Evaluate(fx,p),Evaluate(fy,p),Evaluate(fden,p)*v7>;
   P:=CoordPlace(r);assert Degree(P) eq e;fullsourcecoordinates:=r;
   cl:=Eltseq(Divisor(P)@@gm);
   vals:=[(&+[cl[i]*r[i]:i in [1..#cl]]) mod 6:r in chars];
   chars:=[r:r in chars|(&+[cl[i]*r[i]:i in [1..#cl]]) mod 6 eq 0];assert #chars gt 0;
   tested+:=1;used+:=1;Include(~seen,jkey);
   printf "OCTIC77_FULL_SOURCE_SPLIT exponent=%o Eorder=%o invariants=%o j=%o class=%o previous_values=%o remaining_characters=%o time=%o\n",e,#E,aInvariants(E),Eltseq(jInvariant(E),GF(2)),cl,vals,#chars,Cputime(tt);
   if #chars eq 1 or used ge 8 then break;end if;
  end for;
  if #chars eq 1 or used ge 8 then break;end if;
 end for;
 if #chars eq 1 then break;end if;
end for;
printf "OCTIC77_COVER_IDENTIFICATION_DONE source_places=%o remaining=%o characters=%o complete=%o time=%o\n",tested,#chars,chars,#chars eq 1,Cputime(tt);

assert #chars eq 1;h6:=hom<G->C6|[C6![a]:a in chars[1]]>;
C2:=AbelianGroup([2]);h2:=hom<G->C2|[C2![a mod 2]:a in chars[1]]>;
load "C77/C77_AS_component.m";
print "OCTIC77_FIXED_AS_COEFFICIENT_LOADED";
// Independent Artin-Schreier character verification on Picard/2 generators.
G2,q2:=quo<G|[2*g:g in Generators(G)]>;span2:=sub<G2|>;pe:=PlaceEnumInit(AF);tried:=0;samples:=0;
repeat
 P:=PlaceEnumNext(pe);tried+:=1;assert tried le 150;
 if Valuation(a2,P) lt 0 then continue;end if;
 F:=ResidueClassField(P);actual:=Integers()!Trace(F!Evaluate(a2,P),GF(2));cl:=Divisor(P)@@gm;
 assert actual eq Eltseq(h2(cl))[1];span2:=sub<G2|span2,q2(cl)>;samples+:=1;
 printf "OCTIC77_AS_FROBENIUS degree=%o class=%o value=%o generated=%o total=%o time=%o\n",Degree(P),Eltseq(cl),actual,#span2,#G2,Cputime(tt);
until span2 eq G2;
printf "OCTIC77_AS_CHARACTER_VERIFIED samples=%o all_Picard_mod2=true time=%o\n",samples,Cputime(tt);
print "OCTIC77_FIXED_AS_CHARACTER_VALIDATED";
quit;
