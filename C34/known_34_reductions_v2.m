/* Actual maximal-order prime reductions. Requires existing X=X1(34)/F3. */
load "realizations/data/vanhoeij_octic_input.m";
R34<zz34>:=PolynomialRing(Rationals());
known34divisors:=[];known34orbits:=[];known34sources:=[];
RF34:=recformat<P>;
for ix34 in [2..4] do
 n34,fc34,xc34,yc34:=Explode(OcticExamples[ix34]);assert n34 eq 34;
 f34:=R34!fc34;assert Degree(f34) eq 8 and IsIrreducible(f34);
 K34<a34>:=NumberField(f34);O34:=MaximalOrder(K34);
 xx34:=Evaluate(R34!xc34,a34);yy34:=Evaluate(R34!yc34,a34);
 r34:=(xx34^2*yy34-xx34*yy34+yy34-1)/(xx34^2*yy34-xx34);s34:=(xx34*yy34-yy34+1)/(xx34*yy34);
 b34:=r34*s34*(r34-1);c34:=s34*(r34-1);E34:=EllipticCurve([1-c34,-b34,-b34,0,0]);P34:=E34![0,0];
 assert 34*P34 eq E34!0 and 17*P34 ne E34!0 and 2*P34 ne E34!0;
 dr34:=DivisorGroup(C)!0;primeinfo34:=[];
 for item34 in Decomposition(O34,3) do
  prime34:=item34[1];e34:=item34[2];F34,red34:=ResidueClassField(prime34);
  assert Valuation(Discriminant(E34),prime34) eq 0;
  assert Valuation(b34,prime34) ge 0 and Valuation(c34,prime34) ge 0;
  Ebar34:=EllipticCurve([1-red34(c34),-red34(b34),-red34(b34),0,0]);Pbar34:=Ebar34![0,0];
  assert 34*Pbar34 eq Ebar34!0 and 17*Pbar34 ne Ebar34!0 and 2*Pbar34 ne Ebar34!0;
  pp34:=Places(C(F34)![red34(xx34),red34(yy34)]);assert #pp34 eq 1;
  point34:=pp34[1];assert Degree(point34) eq Degree(F34) and not IsCusp(X,point34);
  assert ModuliPoint(X,Ebar34,rec<RF34|P:=Pbar34>) eq point34;
  dr34+:=e34*Divisor(point34);Append(~primeinfo34,[e34,Degree(F34)]);
 end for;
 assert Degree(dr34) eq 8;
 orb34:=Seqset([DiamondOperator(X,k34,dr34):k34 in [1..17]|GCD(k34,34) eq 1]);
 assert &and[Degree(dd34) eq 8:dd34 in orb34];
 known34divisors cat:=Setseq(orb34);Append(~known34orbits,Setseq(orb34));Append(~known34sources,primeinfo34);
 printf "C34_KNOWN_REDUCTION index=%o primes_ef=%o diamond_divisors=%o\n",ix34,primeinfo34,#orb34;
end for;
known34divisors:=Setseq(Seqset(known34divisors));
printf "C34_KNOWN_REDUCTIONS_DONE classes_not_yet_mapped=true reduced_divisors=%o\n",#known34divisors;
