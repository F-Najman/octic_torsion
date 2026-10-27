Attach("C54/C54_coordinate_place.m");
tv54:=Cputime();data7:=eval Read("C54/data/F5_degree7_pointdata.txt");def7,records7:=Explode(data7);
assert #records7 eq 828;R7:=PolynomialRing(GF(5));f7:=R7!def7;assert IsIrreducible(f7) and Degree(f7) eq 7;
K7<a7>:=ext<GF(5)|f7>;FF5:=FunctionField(C5);good7fresh54:=[];weighted7fresh54:=0;jkeys7fresh54:={};
for r in records7 do
 xycoeff,n,r2,jcoeff:=Explode(r);assert n mod 54 eq 0 and r2 in [1,3];
 xy:=[Evaluate(R7!co,a7):co in xycoeff];jj:=Evaluate(R7!jcoeff,a7);assert Degree(MinimalPolynomial(jj,GF(5))) eq 7;
 jkey:=Eltseq(MinimalPolynomial(jj,GF(5)));assert jkey notin jkeys7fresh54;Include(~jkeys7fresh54,jkey);
 P:=MDPlaceCoordinateFullDegree([FF5.1,FF5.2],xy);EP:=EllipticCurve(X5,P);PP:=LevelStructure(X5,P)`P;
 assert 27*PP eq EP!0 and 9*PP ne EP!0;
 b2e,b4e,b6e,b8e:=Explode(bInvariants(EP));RR<TT>:=PolynomialRing(BaseRing(EP));assert #Roots(4*TT^3+b2e*TT^2+2*b4e*TT+b6e) eq r2;
 assert MinimalPolynomial(jInvariant(EP),GF(5)) eq MinimalPolynomial(jj,GF(5));
 Append(~good7fresh54,P);weighted7fresh54+:=9*r2;
end for;
assert #good7fresh54 eq 828 and #Seqset(good7fresh54) eq 828 and #jkeys7fresh54 eq 828 and weighted7fresh54 eq 11412;
for i in [1..#good7fresh54] do
 TestClass54(b5(Divisor(good7fresh54[i])-7*Divisor(c05)),7,<7,i>,~checked54,~members54,~passing54,~hits54);
end for;
printf "OCTIC54_DEGREE7_REVALIDATED representatives=828 weighted_source_places=11412 checked=%o members=%o passing=%o complete=true time=%o\n",checked54,members54,passing54,Cputime(tv54);
