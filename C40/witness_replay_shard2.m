/* Standalone direct-principal-identity replay.  Explicit cases [3, 7, 11, 15, 19, 23, 27, 31, 35, 39, 43, 47, 51, 55, 59, 63, 67, 71, 75, 79, 83, 87, 91, 95]. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");
Attach("common/fast_hecke_11.m");
load "C40/data/witness_replay_shard2_data.m";
X:=MDX1(40,GF(3));C:=Curve(X);
F0:=AlgorithmicFunctionField(FunctionField(C));F:=RationalExtensionRepresentation(F0);
B:=BaseField(F);kk:=ConstantField(F);PP:=PolynomialRing(kk);
function RatData(a)
 return [[Integers()!x:x in Coefficients(Numerator(a))],[Integers()!x:x in Coefficients(Denominator(a))]];
end function;
assert [RatData(a):a in Coefficients(DefiningPolynomial(F))] eq [[[Integers()!u:u in z]:z in v]:v in C40ReplayFieldEquation];
function GetFunction(data)
 return F![B!(PP!v[1])/B!(PP!v[2]):v in data];
end function;
OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);
fb:=[**];
for j in [1..#C40ReplayIdealBases] do
 OO:=C40ReplayFiniteFlags[j] eq 1 select OF else OI;
 I:=ideal<OO|[GetFunction(a):a in C40ReplayIdealBases[j]]>;
 assert Degree(Place(I)) eq C40ReplayDegrees[j];Append(~fb,I);
end for;
assert #fb eq 616;
function DivisorIdeals(sparse)
 IF:=ideal<OF|1>;II:=ideal<OI|1>;degree:=0;seen:={Integers()|};
 for row in sparse do
  j,e:=Explode(row);assert j in [1..616] and e gt 0 and j notin seen;Include(~seen,j);
  degree+:=e*C40ReplayDegrees[j];
  if C40ReplayFiniteFlags[j] eq 1 then IF*:=fb[j]^e;else II*:=fb[j]^e;end if;
 end for;
 return IF,II,degree;
end function;
vA:=C40ReplayCases[1][2][6];AF,AI,degree:=DivisorIdeals(vA);assert degree eq 1;
basecusps:=[];
for cusp in Cusps(X) do
 if Degree(cusp) ne 1 then continue;end if;
 I,J:=Ideals(DivisorGroup(F)!FunctionFieldDivisor(Divisor(cusp)));
 if I eq AF and J eq AI then Append(~basecusps,cusp);end if;
end for;
assert #basecusps eq 1;c0:=basecusps[1];
K:=ext<GF(3)|PolynomialRing(GF(3))!C40ReplayFieldPolynomial>;
for item in C40ReplayCases do
 index:=item[1];w:=item[2];assert w[2] eq 0 and w[6] eq vA;
 EF,EI,ed:=DivisorIdeals(w[4]);RF,RI,rd:=DivisorIdeals(w[5]);assert ed eq 25 and rd eq 25;
 b:=GetFunction(w[7]);assert b ne 0;
 bc:=[Evaluate(PolynomialRing(K)!v,K.1):v in C40ReplayTateBC[index]];
 P:=MDPlaceCommonZeros(X`_coordinates,bc);assert Degree(P) eq 8;
 DD:=DivisorGroup(F)!FunctionFieldDivisor(Divisor(P)-8*Divisor(c0));
 ID,JD:=Ideals(DD);
 assert EF/(RF*(b*OF)) eq ID and EI/(RI*(b*OI)) eq JD;
 printf "C40_INDEPENDENT_PRINCIPAL_WITNESS index=%o both_orders_exact=true time=%o\n",index,Cputime(tt);
end for;
printf "C40_INDEPENDENT_PRINCIPAL_REPLAY_DONE cases=%o classgroup_search=false time=%o\n",#C40ReplayCases,Cputime(tt);
quit;
