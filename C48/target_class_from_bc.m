/* Uses the current proven target presentation and exact RR map. */
function C48TargetFromBC(bc,dd)
 K:=Parent(bc[1]);assert #K eq 5^dd;
 b,c:=Explode(bc);E0:=EllipticCurve([1-c,-b,-b,K!0,K!0]);P0:=E0![0,0];
 E2,L2:=C48TwoIsogenousModuli(E0,P0);
 assert aInvariants(E2)[1] eq 0 and aInvariants(E2)[3] eq 0;
 targetbc:=DerickxNormalForm_bc(E2,L2`P,L2`Q);
 yy:=C48PairPlaceLinear(Y`_coordinates,targetbc);
 assert dd mod Degree(yy) eq 0 and dd div Degree(yy) in [1,2];
 return yy;
end function;
targetClassCache48:=AssociativeArray();
procedure C48TargetClass(bc,dd,~cache,~answer)
 yy:=C48TargetFromBC(bc,dd);ee:=Degree(yy);
 if not IsDefined(cache,yy) then
  cache[yy]:=phiRR48(FunctionFieldDivisor(Divisor(yy)-ee*Divisor(c0)));
 end if;
 answer:=(dd div ee)*cache[yy];
end procedure;
