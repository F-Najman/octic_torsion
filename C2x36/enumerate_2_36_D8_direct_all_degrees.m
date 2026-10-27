/* Exact torsion-first source enumeration, without point counts or large place lists.
   Limit=0 means full; a positive Limit is explicitly an incomplete timing probe. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
if not assigned Limit then Limit:=1000;end if;if Type(Limit) eq MonStgElt then Limit:=StringToInteger(Limit);end if;
if not assigned Tag then Tag:="profile1000";end if;
if not assigned D then D:=8;end if;if Type(D) eq MonStgElt then D:=StringToInteger(D);end if;assert D in [3..8];
K:=GF(5^D);Fp:=GF(5);R<x>:=PolynomialRing(K);
function BC(E,P,Q)
 ai:=aInvariants(E);assert ai[1] eq 0 and ai[3] eq 0 and P[2] eq 0;
 xx:=Q[1]-P[1];yy:=Q[2];return [xx^3/yy^2,(ai[2]+3*P[1])*xx^2/yy^2];
end function;
function Key(v)
 keys:=[];
 for i in [1..D] do Append(~keys,&cat[[Integers()!a:a in Eltseq(z,Fp)]:z in v]);v:=[z^5:z in v];end for;
 return Sort(keys)[1];
end function;
function Orb(E,P,Q)
 return [Key(BC(E,P+18*b*Q,c*P+d*Q)):b,c in [0,1],d in [1,5,7,11,13,17]];
end function;
function Deg(v)
 for d in Divisors(D) do if &and[z^(5^d) eq z:z in v] then return d;end if;end for;error "fielddegree";
end function;
seen:={K|};orbs:=AssociativeArray();points:=AssociativeArray();records:=[* *];jcount:=0;curves:=0;good4:=0;good9:=0;divtime:=0.;hashTime:=0.;complete:=true;
for j in K do
 if j in seen then continue;end if;
 if Limit gt 0 and jcount ge Limit then complete:=false;break;end if;
 for i in [0..D-1] do Include(~seen,j^(5^i));end for;jcount+:=1;
 es:=Twists(EllipticCurveWithjInvariant(j));assert #es in [2,4,6];
 for ie in [1..#es] do
  E:=es[ie];curves+:=1;aa:=aInvariants(E);
  E:=EllipticCurve([K!0,aa[2]+aa[1]^2/4,K!0,aa[4]+aa[1]*aa[3]/2,aa[5]+aa[3]^2/4]);
  ai:=aInvariants(E);assert ai[1] eq 0 and ai[3] eq 0;
  cubic:=x^3+ai[2]*x^2+ai[4]*x+ai[5];roots:=Roots(cubic);
  if #roots ne 3 then if #es eq 2 then break;else continue;end if;end if;
  e2:=[r[1]:r in roots];fours:=[* *];
  for i in [1..3] do
   other:=[e2[k]:k in [1..3]|k ne i];a:=e2[i]-other[1];b:=e2[i]-other[2];
   if not IsSquare(a) or not IsSquare(b) then continue;end if;
   u:=Sqrt(a);v:=Sqrt(b);Q4:=E![e2[i]+u*v,u*v*(u+v)];P2:=E![other[1],0];
   assert 4*Q4 eq E!0 and 2*Q4 eq E![e2[i],0] and P2 ne 2*Q4;Append(~fours,[*P2,Q4*]);
  end for;
  if #fours eq 0 then continue;end if;good4+:=1;t:=Cputime();psi:=DivisionPolynomial(E,9);rr:=Roots(psi);divtime+:=Cputime(t);
  nines:=[* *];sg:={};
  for root in rr do
   xx:=root[1];val:=Evaluate(cubic,xx);if not IsSquare(val) then continue;end if;
   Q9:=E![xx,Sqrt(val)];assert 9*Q9 eq E!0;if 3*Q9 eq E!0 then continue;end if;
   ky:=Sort([[Integers()!a:a in Eltseq((k*Q9)[1],Fp)]:k in [1,2,4]]);
   if ky in sg then continue;end if;Include(~sg,ky);Append(~nines,Q9);
  end for;
  if #nines eq 0 then continue;end if;good9+:=1;assert #nines in [1,3,12];t:=Cputime();
  for four in fours do for Q9 in nines do
   P,Q4:=Explode(four);Q:=Q4+7*Q9;assert 36*Q eq E!0 and 18*Q ne E!0 and 12*Q ne E!0 and P ne 18*Q;
   ims:=Orb(E,P,Q);orbkey:=Sort(ims)[1];
   if IsDefined(orbs,orbkey) then continue;end if;
   orbs[orbkey]:=#records+1;qt:=[Evaluate(R!orbkey[1..D],K.1),Evaluate(R!orbkey[D+1..2*D],K.1)];d:=Deg(qt);
   for im in Seqset(ims) do if IsDefined(points,im) then assert points[im] eq d;else points[im]:=d;end if;end for;
   Append(~records,[*d,orbkey*]);
  end for;end for;hashTime+:=Cputime(t);
 end for;
 if jcount mod 500 eq 0 then print "M2N36_D8_DIRECT_PROGRESS",jcount,curves,good4,good9,#records,#Keys(points),divtime,hashTime,Cputime(t0);end if;
end for;
counts:=[#[ky:ky in Keys(points)|points[ky] eq d]:d in Divisors(D)];
expected:=[0,0,48,120,720,2904,10416,50334];
if complete then assert #seen eq #K;assert counts eq [expected[d]:d in Divisors(D)];end if;
Write(Sprintf("C2x36/data/M2N36_D8_direct_%o_arithmetic.txt",Tag),Sprint([*D,Coefficients(DefiningPolynomial(K,Fp)),complete,jcount,curves,good4,good9,counts,records*]):Overwrite:=true);
print "M2N36_D8_DIRECT_DONE",D,complete,jcount,curves,good4,good9,counts,#records,divtime,hashTime,Cputime(t0);quit;
