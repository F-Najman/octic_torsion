/* Exact maximal-order reduction of the known octic point. Requires X1(40)/F3. */
load "realizations/data/vanhoeij_octic_input.m";
R40<zz40>:=PolynomialRing(Rationals());
known40divisors:=[];known40orbits:=[];known40sources:=[];
RF40:=recformat<P>;
for ix40 in [7] do
 n40,fc40,xc40,yc40:=Explode(OcticExamples[ix40]);assert n40 eq 40;
 f40:=R40!fc40;assert Degree(f40) eq 8 and IsIrreducible(f40);
 K40<a40>:=NumberField(f40);O40:=MaximalOrder(K40);
 xx40:=Evaluate(R40!xc40,a40);yy40:=Evaluate(R40!yc40,a40);
 r40:=(xx40^2*yy40-xx40*yy40+yy40-1)/(xx40^2*yy40-xx40);s40:=(xx40*yy40-yy40+1)/(xx40*yy40);
 b40:=r40*s40*(r40-1);c40:=s40*(r40-1);E40:=EllipticCurve([1-c40,-b40,-b40,0,0]);P40:=E40![0,0];
 assert 40*P40 eq E40!0 and 20*P40 ne E40!0 and 8*P40 ne E40!0;
 dr40:=DivisorGroup(C)!0;primeinfo40:=[];
 for item40 in Decomposition(O40,3) do
  prime40:=item40[1];e40:=item40[2];F40,red40:=ResidueClassField(prime40);
  assert Valuation(Discriminant(E40),prime40) eq 0;
  assert Valuation(b40,prime40) ge 0 and Valuation(c40,prime40) ge 0;
  Ebar40:=EllipticCurve([1-red40(c40),-red40(b40),-red40(b40),0,0]);Pbar40:=Ebar40![0,0];
  assert 40*Pbar40 eq Ebar40!0 and 20*Pbar40 ne Ebar40!0 and 8*Pbar40 ne Ebar40!0;
  pp40:=Places(C(F40)![red40(xx40),red40(yy40)]);assert #pp40 eq 1;
  point40:=pp40[1];assert Degree(point40) eq Degree(F40) and not IsCusp(X,point40);
  assert ModuliPoint(X,Ebar40,rec<RF40|P:=Pbar40>) eq point40;
  dr40+:=e40*Divisor(point40);Append(~primeinfo40,[e40,Degree(F40)]);
 end for;
 assert Degree(dr40) eq 8;
 orb40:=Seqset([DiamondOperator(X,k40,dr40):k40 in [1..20]|GCD(k40,40) eq 1]);
 assert &and[Degree(dd40) eq 8:dd40 in orb40];
 known40divisors cat:=Setseq(orb40);Append(~known40orbits,Setseq(orb40));Append(~known40sources,primeinfo40);
 printf "C40_KNOWN_REDUCTION index=%o primes_ef=%o diamond_divisors=%o\n",ix40,primeinfo40,#orb40;
end for;
known40divisors:=Setseq(Seqset(known40divisors));
printf "C40_KNOWN_REDUCTIONS_DONE classes_not_yet_mapped=true reduced_divisors=%o\n",#known40divisors;
