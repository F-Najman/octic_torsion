/* Forget the sign of the independent3-point over X1(2,14). */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Y:=MDX11(2,14,GF(5));C:=Curve(Y);A,am:=AlgorithmicFunctionField(FunctionField(C));
B:=am(Y`_coordinates[1]);CC:=am(Y`_coordinates[2]);D:=B*(1-B-CC);
E14:=EllipticCurve([A|0,CC,0,D,0]);P2:=E14![0,0];Q14:=E14![B,B];
assert 14*Q14 eq E14!0 and 7*Q14 ne E14!0 and 2*Q14 ne E14!0 and 7*Q14 ne P2;
R<U>:=PolynomialRing(A);f:=DivisionPolynomial(E14,3);assert Degree(f) eq 4 and IsIrreducible(f);
K<x3>:=FunctionField(f);assert Genus(A) eq 4 and Genus(K) eq 31;
print "C242_FORGET3SIGN_GENUS31_MODEL",Genus(K),Cputime(t0);

Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
load "C2x42/af_place_three_via_base.m";
if not assigned Limit then Limit:=4;end if;
if Type(Limit) eq MonStgElt then Limit:=StringToInteger(Limit);end if;
data:=eval Read("C2x42/data/M2N42_41_remaining_supports.m");
R5:=PolynomialRing(GF(5));fieldobjects:=AssociativeArray();
for row in data[1] do fieldobjects[row[1]]:=ext<GF(5)|R5![Integers()!a:a in row[2]]>;end for;
function Source(D,key)
 key:=[Integers()!a:a in key];
 L:=fieldobjects[D];b:=Evaluate(R5!key[1..D],L.1);c:=Evaluate(R5!key[D+1..2*D],L.1);
 E:=EllipticCurve([L|0,c,0,b*(1-b-c),0]);return E,E![0,0],E![b,b];
end function;
function ZValues(E,P,Q)
 aa:=aInvariants(E);assert aa[1] eq 0 and aa[3] eq 0 and P[2] eq 0;
 Q14:=3*Q;Q3:=14*Q;xx:=Q14[1]-P[1];yy:=Q14[2];uu:=yy/xx;
 bc:=[xx^3/yy^2,(aa[2]+3*P[1])*xx^2/yy^2];zz:=(Q3[1]-P[1])/uu^2;
 DD:=bc[1]*(1-bc[1]-bc[2]);assert 3*zz^4+4*bc[2]*zz^3+6*DD*zz^2-DD^2 eq 0;
 return bc cat [zz];
end function;
coords:=[K!B,K!CC,x3];cache:=AssociativeArray();export:=[* *];results:=[* *];
procedure ImageA(~answer,~cache,~export,D,key,ell)
 key:=[Integers()!a:a in key];
 ck:=[ell,D] cat key;if IsDefined(cache,ck) then answer:=cache[ck];return;end if;
 E,P,Q:=Source(D,key);vv:=ZValues(E,P,Q);p:=MDPlaceAFThreeViaBase(coords,vv,[B,CC]);ee:=Degree(p);assert D mod ee eq 0;
 N:=(D div ee)*Divisor(p);pp:=MDPlaceAFThreeViaBase(coords,ZValues(E,P,ell*Q),[B,CC]);ef:=Degree(pp);assert D mod ef eq 0;
 delta:=(D div ef)*Divisor(pp);isogs:=MDIsogeniesFastPoints(E,[P,Q],ell);T:=DivisorGroup(K)!0;isogdata:=[* *];
 print "C242_G31_HECKE_BEGIN",D,ell,#isogs,Cputime(t0);
 for item in isogs do
  EG,pts,k:=Explode(item);v:=ZValues(EG,pts[1],pts[2]);z:=MDPlaceAFThreeViaBase(coords,v,[B,CC]);e:=Degree(z);assert (D*k) mod e eq 0;
  mult:=(D*k) div e;T+:=mult*Divisor(z);Append(~isogdata,[*k,e,mult*]);
  print "C242_G31_HECKE_IMAGE",D,ell,k,e,mult,Cputime(t0);
 end for;
 assert Degree(T) eq D*(ell+1);answer:=T-N-ell*delta;assert Degree(answer) eq 0;cache[ck]:=answer;
 Append(~export,[*ell,D,key,ee,D div ee,ef,D div ef,isogdata*]);
end procedure;

row:=[rr:rr in data[2]|rr[1] eq 193][1];assert row[2] eq [*2,4,5,4,38*];
D:=row[3][1][1];key:=row[3][1][2];assert D eq 4;
E,P,Q:=Source(D,key);ell:=11;T:=DivisorGroup(K)!0;
for item in MDIsogeniesFastPoints(E,[P,Q],ell) do
 EG,pts,k:=Explode(item);z:=MDPlaceAFThreeViaBase(coords,ZValues(EG,pts[1],pts[2]),[B,CC]);e:=Degree(z);assert (D*k) mod e eq 0;
 T+:=((D*k) div e)*Divisor(z);
end for;
ND:=DivisorGroup(K)!0;counts:=AssociativeArray();native:=MDIsogenies(E,ell);assert #native eq ell+1;
for phi in native do
 EE:=Domain(phi);pp:=phi(EE![BaseRing(EE)!a:a in Eltseq(P)]);qq:=phi(EE![BaseRing(EE)!a:a in Eltseq(Q)]);
 z:=MDPlaceAFThreeViaBase(coords,ZValues(Codomain(phi),pp,qq),[B,CC]);e:=Degree(z);
 if IsDefined(counts,z) then counts[z]+:=1;else counts[z]:=1;end if;
end for;
for z in Keys(counts) do e:=Degree(z);assert (D*counts[z]) mod e eq 0;ND+:=((D*counts[z]) div e)*Divisor(z);end for;
assert ND eq T and Degree(T) eq 48;
print "C242_G31_NATIVE_MD_ISOGENIES_REFERENCE_EQUAL",ell,D,Degree(T),Cputime(t0);
print "C242_G31_NATIVE_REFERENCE_DONE complete=true",Cputime(t0);quit;
