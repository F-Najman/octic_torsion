AttachSpec("mdmagma/v2/mdmagma.spec");
SetColumns(0);SetQuitOnError(true);
try
SetMemoryLimit(16*10^9);
if assigned N then N:=StringToInteger(N);else N:=57;end if;
assert N eq 57;m:=N eq 55 select 5 else 3;ell:=N div m;h:=(ell-1) div 2;tstart:=Cputime();
CQ:=SmallModularCurve(ell);jraw:=jInvariant(CQ,ell);wraw:=DefiningPolynomials(AtkinLehnerInvolution(CQ,ell,ell));
CC:=ChangeRing(CQ,GF(2));assert IsNonsingular(CC) and Genus(CC) eq (ell eq 23 select 2 else 1);
RC:=CoordinateRing(Ambient(CC));jn:=RC!Numerator(jraw);jd:=RC!Denominator(jraw);ww:=[RC!e:e in wraw];
jell:=FunctionField(CC)!(jn/jd);jNell:=FunctionField(CC)!(Evaluate(jn,ww)/Evaluate(jd,ww));
K,mK:=AlgorithmicFunctionField(FunctionField(CC));jK:=mK(jell);
jdegree:=-&+[Valuation(jK,P)*Degree(P):P in Poles(jK)];assert jdegree eq ell+1;
printf "DATABASE_ELL_MODEL N=%o ell=%o genus=%o good_reduction=true j_degree=%o time=%o\n",N,ell,Genus(K),jdegree,Cputime(tstart);
R<v>:=PolynomialRing(K);f:=N eq 55 select v^6+jK*v+jK else v^4+v^3+1/jK;assert IsIrreducible(f);
AF<t3>:=FunctionField(f);gg:=Genus(AF);assert gg eq (N eq 69 select 7 else 5);
O:=MaximalOrderFinite(AF);OK:=MaximalOrderFinite(K);
numY,denY:=IntegralSplit(K.1,OK);numJ,denJ:=IntegralSplit(jK,OK);
IX:=AF!BaseField(K).1;IY:=AF!K!numY;multT:=N eq 55 select K!denJ else K!numJ;IT:=(AF!multT)*t3;
integral_coordinates:=[O!IX,O!IY,O!IT];
fx:=(K!BaseField(K).1)@@mK;fy:=(K!numY)@@mK;fJ:=multT@@mK;
printf "X0_3ELL_MODEL N=%o genus=%o degree_over_X0ell=4 integral_coordinates_ready=true time=%o\n",N,gg,Cputime(tstart);
function ThreeParameter(E,R3)
 a1,a2,a3,a4,a6:=Explode(aInvariants(E));x:=R3[1];y:=R3[2];
 aa3:=2*y+a3+a1*x;aa2:=3*x+a2;aa4:=3*x^2+2*x*a2+a4-a1*y;
 s:=aa4/aa3;aone:=a1+2*s;assert aa2-s^2-a1*s eq 0;
 tt:=aa3/aone^3;assert 1/(tt^3*(tt+1)) eq jInvariant(E);return tt;
end function;
function SmallParameter(E,R)
 if N eq 55 then
  b,c:=Explode(TateNormalForm_bc(E,R));assert b eq c;u:=b+1/b;assert u^6/(u+1) eq jInvariant(E);return u;
 else return ThreeParameter(E,R);end if;
end function;
function EllPoint(E,Rl)
 L:=BaseRing(E);PL<X>:=PolynomialRing(L);polyker:=&*[X-(k*Rl)[1]:k in [1..(ell-1) div 2]];
 Et,phi:=IsogenyFromKernel(E,polyker);jj:=jInvariant(E);jjN:=jInvariant(Et);
 A:=Evaluate(jn,[X,PL!0,PL!1]);BB:=Evaluate(jn,[X,PL!1,PL!1])-A;DD:=Evaluate(jd,[X,PL!0,PL!1]);
 assert Evaluate(jd,[X,PL!1,PL!1]) eq DD;
 ce:=DefiningEquation(CC);c0:=Evaluate(ce,[X,PL!0,PL!1]);c2:=PL!1;
 if ell eq 19 then c1:=PL!1;else c1:=X^3+X+1;end if;
 assert Evaluate(ce,[X,PL!1,PL!1]) eq c0+c1+c2;
 YY:=jj*DD-A;res:=c2*YY^2+c1*YY*BB+c0*BB^2;
 out:=[];
 for rt in Roots(res) do
  xx:=rt[1];bb:=Evaluate(BB,xx);
  if bb ne 0 then ys:=[Evaluate(YY,xx)/bb];else ys:=[zz[1]:zz in Roots(c2*X^2+Evaluate(c1,xx)*X+Evaluate(c0,xx))];end if;
  for yy in ys do
   pt:=CC(L)![xx,yy,L!1];
   if Evaluate(jd,Eltseq(pt)) ne 0 and Evaluate(jell,pt) eq jj and Evaluate(jNell,pt) eq jjN then Append(~out,pt);end if;
  end for;
 end for;
 assert #out eq 1;return out[1];
end function;
function Coordinates0(E,R3,Rl)
 pt:=EllPoint(E,Rl);return <Evaluate(fx,pt),Evaluate(fy,pt),Evaluate(fJ,pt)*SmallParameter(E,R3)>;
end function;
lookup:=AssociativeArray();models:=AssociativeArray();
function Key0(r,e)
 return Minimum([&cat[Eltseq(s^(2^i),GF(2)):s in r]:i in [0..e-1]]);
end function;
for e in [7,8] do
 F<a>:=GF(2^e);count:=0;before:=#Keys(lookup);
 for j in F do
  E0:=EllipticCurveFromjInvariant(j);
  for E in Twists(E0) do
   if #E mod 57 ne 0 then continue;end if;
   assert #E eq (e eq 7 select 114 else 228);count+:=1;
   A,mp:=AbelianGroup(E);assert Invariants(A) eq [#E];PN:=(#E div 57)*mp(A.1);assert Order(PN) eq 57;
   R3:=19*PN;Rl:=3*PN;r:=Coordinates0(E,R3,Rl);key:=<e,Key0(r,e)>;
   lookup[key]:=r;models[key]:=<E,R3,Rl>;
  end for;
 end for;
 expected:=e eq 7 select 2 else 1;
 assert count eq (e eq 7 select 14 else 8) and #Keys(lookup)-before eq expected;
 printf "C57_BASE_NC_COVERAGE degree=%o curves=%o closed_image_points=%o expected=%o time=%o\n",e,count,#Keys(lookup)-before,expected,Cputime(tstart);
end for;
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

cs:=Poles(AF!jK);assert #cs eq 4 and &and[Degree(c) eq 1:c in cs];
Dzero:=DivisorGroup(AF)!0;RG,rm:=RayClassGroup(Dzero);Q,q:=quo<RG|[2*g:g in Generators(RG)]>;
assert Invariants(RG) eq [15,30,0] and Invariants(Q) eq [2,2];
cols:=[Vector(GF(2),Eltseq(q(Divisor(c)@@rm))):c in cs];
expected:=[GF(2)|(-Valuation(AF!jK,c) in [3,19]) select 1 else 0:c in cs];
chars:=[v:v in VectorSpace(GF(2),Ngens(Q))|&and[InnerProduct(v,cols[i]) eq expected[i]:i in [1..4]]];
print "C57_COVER_CUSP_CHARACTERS",#chars,chars;assert #chars eq 1;
ps:=[CoordPlace(lookup[key]):key in Sort(Setseq(Keys(lookup)))];
pcs:=[Vector(GF(2),Eltseq(q(Divisor(P)@@rm))):P in ps];
chars:=[v:v in chars|&and[InnerProduct(v,cl) eq 0:cl in pcs]];
assert #chars eq 1;chi:=chars[1];
print "C57_COVER_IDENTIFIED",chi,"cusp_widths",[-Valuation(AF!jK,c):c in cs],"cusp_values",expected,"nc7_classes",pcs,Cputime(tstart);
C2:=AbelianGroup([2]);hh:=hom<RG->C2|[C2![Integers()!InnerProduct(chi,Vector(GF(2),Eltseq(q(RG.i))))]:i in [1..Ngens(RG)]]>;
assert #Image(hh) eq 2; assert &and[hh(Divisor(cs[i])@@rm) eq C2![Integers()!expected[i]]:i in [1..#cs]]; assert &and[hh(Divisor(P)@@rm) eq C2!0:P in ps]; print "COVER_HOM_VALUES_VERIFIED",[hh(RG.i):i in [1..Ngens(RG)]]; U:=Kernel(hh);AA:=AbelianExtension(Dzero,U);
print "C57_COVER_CLASS_FIELD_START",Cputime(tstart);
CYns:=FunctionField(AA);polys:=DefiningPolynomials(CYns);assert #polys eq 1 and Degree(polys[1]) eq 2;CY:=FunctionField(polys[1]);assert Genus(CY) eq 9;
print "C57_CLASS_FIELD_READY",Genus(CY),Cputime(tstart);

function RelativeDecomposition(F,P)
 assert BaseField(F) cmpeq FunctionField(P);
 f,g:=TwoGenerators(P);u:=UniformizingElement(P);
 assert Valuation(f,P) gt 0 and Valuation(g,P) gt 0 and Valuation(u,P) eq 1;
 raw:=Decomposition(F,P);
 out:=[Q:Q in raw|Valuation(F!f,Q) gt 0 and Valuation(F!g,Q) gt 0];
 assert #out gt 0 and &+[Degree(Q):Q in out] eq 2*Degree(P);
 assert &and[Valuation(F!u,Q) eq 1:Q in out];
 return out;
end function;
// The genus equality, degree two and separability imply zero different
// by Riemann--Hurwitz; a split rational cusp checks the exact constants.
cN:=[c:c in cs|-Valuation(AF!jK,c) eq N][1];dcN:=RelativeDecomposition(CY,cN);
assert #dcN eq 2 and &and[Degree(a) eq 1:a in dcN];
for i in [1..Ngens(RG)] do
 Dg:=rm(RG.i);fr:=0;
 for B in Support(Dg) do
  dd:=RelativeDecomposition(CY,B);
  if #dd eq 1 then assert Degree(dd[1]) eq 2*Degree(B);fr+:=Valuation(Dg,B);
  else assert #dd eq 2 and &and[Degree(a) eq Degree(B):a in dd];end if;
 end for;
 assert hh(RG.i) eq C2![fr mod 2];
 print "COVER_GENERATOR_FROBENIUS_VERIFIED",i,fr mod 2;
end for;
print "COVER_CHARACTER_VERIFIED",N,Genus(CY),Cputime(tstart);
tcg:=Cputime();cN:=[c:c in cs|-Valuation(AF!jK,c) eq N][1];cyc:=RelativeDecomposition(CY,cN);assert #cyc eq 2 and Degree(cyc[1]) eq 1;
GY,_,phY:=ClassGroup(CY:ReductionDivisor:=Divisor(cyc[1]),Proof:=true);assert Invariants(GY) eq [105,210,0];
print "C57_CLASS_GROUP_READY",Invariants(GY),Cputime(tcg);

assert assigned GY and assigned phY;
placecache:=AssociativeArray();
procedure CachedPlaceClass(P,~cache,~value)
 if not IsDefined(cache,P) then cache[P]:=phY(Divisor(P));print "C57_POINT_CLASS",Degree(P),Eltseq(cache[P]),Cputime(tstart);end if;value:=cache[P];
end procedure;

assert KroneckerSymbol(57,5) eq -1;
idx:=0;bad:=0;total:=0;arithmetic_records:=[**];
for key in Sort(Setseq(Keys(lookup))) do
 idx+:=1;e:=key[1];r:=lookup[key];md:=models[key];E:=md[1];R3:=md[2];Rl:=md[3];P:=CoordPlace(r);assert Degree(P) eq e;
 dp:=RelativeDecomposition(CY,P);assert #dp eq 2 and &and[Degree(t) eq e:t in dp];Py:=dp[1];
 pclass:=phY(Divisor(dp[1]));spclass:=phY(Divisor(dp[2]));target_records:=[**];test_records:=[];
 isogs:=MDIsogenies(E,5);assert #isogs eq 6;out:=[];
 for ph in isogs do Append(~out,CoordPlace(Coordinates0(Codomain(ph),ph(R3),ph(Rl))));end for;
 choices:=[<DivisorGroup(CY)!0,DivisorGroup(CY)!0>];base_degree:=0;
 for B in Seqset(out) do
  mult:=#[z:z in out|z eq B];assert (e*mult) mod Degree(B) eq 0;mult:=(e*mult) div Degree(B);base_degree+:=mult*Degree(B);
  dd:=RelativeDecomposition(CY,B);
  assert #dd in [1,2];
  if #dd eq 2 then
   assert &and[Degree(t) eq Degree(B):t in dd];
   Append(~target_records,[*Degree(B),mult,1,Eltseq(phY(Divisor(dd[1]))),Eltseq(phY(Divisor(dd[2])))*]);
   options:=[<a*Divisor(dd[1])+(mult-a)*Divisor(dd[2]),a*Divisor(dd[2])+(mult-a)*Divisor(dd[1])>:a in [0..mult]];
  else
   assert Degree(dd[1]) eq 2*Degree(B) and IsEven(mult);
   Append(~target_records,[*Degree(B),mult,2,Eltseq(phY(Divisor(dd[1]))),Eltseq(phY(Divisor(dd[1])))*]);
   options:=[<(mult div 2)*Divisor(dd[1]),(mult div 2)*Divisor(dd[1])>];
  end if;
  choices:=[<T[1]+D[1],T[2]+D[2]>:T in choices,D in options];
 end for;
 assert base_degree eq 6*e;
 survivors:=0;printf "C57_GOOD78_CASE_START index=%o degree=%o lift_choices=%o time=%o\n",idx,e,#choices,Cputime(tstart);
 for ii in [1..#choices] do
  T:=choices[ii][1];Td:=choices[ii][2];assert Degree(T) eq 6*e and Degree(Td) eq 6*e;
  A:=Td-T+4*(Divisor(dp[2])-Divisor(Py));assert Degree(A) eq 0;
  value:=GY!0;
  for B in Support(A) do cv:=GY!0;CachedPlaceClass(B,~placecache,~cv);value+:=Valuation(A,B)*cv;end for;
  pr:=value eq GY!0;
  if total eq 0 or pr then assert IsPrincipal(A) eq pr;end if;
  printf "C57_GOOD78_OPERATOR_CLASS index=%o degree=%o choice=%o value=%o principal=%o\n",idx,e,ii,Eltseq(value),pr;
  Append(~test_records,Eltseq(value));total+:=1;if pr then survivors+:=1;end if;
 end for;
 Append(~arithmetic_records,[*idx,e,key,Eltseq(pclass),Eltseq(spclass),target_records,test_records*]);
 printf "C57_GOOD78_CASE_DONE index=%o degree=%o survivors=%o time=%o\n",idx,e,survivors,Cputime(tstart);
 if survivors gt 0 then bad+:=1;end if;
end for;
Write("C57/data/C57_good78_class_arithmetic.txt",Sprint([*Invariants(GY),[Eltseq(r):r in Rows(RelationMatrix(GY))],arithmetic_records*]):Overwrite:=true);
printf "C57_GOOD78_DONE cases=%o tests=%o remaining=%o complete=true time=%o\n",idx,total,bad,Cputime(tstart);

catch err
 print "C57_GOOD78_FAILED";print err`Object;
end try;
