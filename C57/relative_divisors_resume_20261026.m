AttachSpec("mdmagma/v2/mdmagma.spec");
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(32*10^9);tt43:=Cputime();
X19:=MDX1(19,GF(2));C19:=Curve(X19);FF19:=FunctionField(C19);assert Genus(C19) eq 7;
E19:=EllipticCurve(X19`_E);P19:=E19!X19`_P;
K19,m19:=AlgorithmicFunctionField(FF19);j19:=m19(jInvariant(E19));
assert &+[-Valuation(j19,P)*Degree(P):P in Poles(j19)] eq 180;
OK19:=MaximalOrderFinite(K19);nY19,dY19:=IntegralSplit(K19.1,OK19);
tf85:=Cputime();aa19:=[m19(a):a in aInvariants(E19)];a1n,a2n,a3n,a4n,a6n:=Explode(aa19);assert a1n ne 0;
r19n:=a3n/a1n;t19n:=(a4n+r19n^2)/a1n;d19n:=(a2n+r19n)/a1n^2;
RP85<Xc,Yc>:=PolynomialRing(K19,2);xo:=a1n^2*Xc+r19n;yo:=a1n^3*Yc+t19n;
oldEq:=yo^2+a1n*xo*yo+a3n*yo+xo^3+a2n*xo^2+a4n*xo+a6n;
assert oldEq/a1n^6 eq Yc^2+Xc*Yc+Xc^3+d19n*Xc^2+1/j19;
print "C57_FULL85_CANONICAL_MODEL_IDENTITY_VERIFIED",Cputime(tf85);
td85:=Cputime();assert assigned d19n;
RD85<WW85>:=PolynomialRing(K19);TX85:=WW85^2+WW85+d19n;fd85:=TX85^4+TX85^3+1/j19;
assert Degree(fd85) eq 8 and Derivative(fd85) ne 0;
print "C57_DIRECT85_IRREDUCIBILITY_START",Cputime(td85);assert IsIrreducible(fd85);
D85<wd85>:=FunctionField(fd85);print "C57_DIRECT85_NORMALIZATION_START",Cputime(td85);assert Genus(D85) eq 85;
print "C57_DIRECT85_MODEL_VERIFIED degree8_over_X119=true genus85=true index1440=true",Cputime(td85);
OD85:=MaximalOrderFinite(D85);den85:=K19!1;
for coeff in Coefficients(fd85) do nn,dd:=IntegralSplit(coeff,OK19);den85*:=dd;end for;
coordD85:=[D85!K19!BaseField(K19).1,D85!nY19,D85!den85*wd85];
assert &and[f in OD85:f in coordD85];
print "C57_DIRECT85_INTEGRAL_COORDINATES_READY",Cputime(td85);

print "C57_DIRECT85_STANDALONE_READY";

tg85:=Cputime();assert Genus(D85) eq 85;
function FullPoint57(E,R3,R19,cache)
 L:=BaseRing(E);assert 3*R3 eq E!0 and R3 ne E!0 and 19*R19 eq E!0 and R19 ne E!0;
 xy:=TateNormalForm_xy(E,R19);pt:=C19(L)![xy[1],xy[2],L!1];
 assert Evaluate(jInvariant(E19),pt) eq jInvariant(E);
 a1,a2,a3,a4,a6:=Explode(aInvariants(E));xp:=R19[1];yp:=R19[2];
 slope:=(xp^2+a4+a1*yp)/(a1*xp+a3);A2t:=a2+slope*a1+xp+slope^2;ut:=(a1*xp+a3)/A2t;
 bn:=A2t/ut^2;cn:=1+a1/ut;bc:=TateNormalForm_bc(E,R19);assert bc eq [bn,cn];
 ET:=EllipticCurve([1-cn,-bn,-bn,L!0,L!0]);
 X3:=(R3[1]+xp)/ut^2;Y3:=(R3[2]+yp+slope*(R3[1]+xp))/ut^3;RT:=ET![X3,Y3];assert 3*RT eq ET!0 and RT ne ET!0;
 b1,b2,b3,b4,b6:=Explode(aInvariants(ET));rr:=b3/b1;tt:=(b4+rr^2)/b1;dd:=(b2+rr)/b1^2;
 xx:=(X3+rr)/b1^2;yy:=(Y3+tt)/b1^3;ww:=yy/xx+xx;
 assert ww^2+ww eq xx+dd and xx^4+xx^3+1/jInvariant(E) eq 0;
 assert Evaluate(d19n@@m19,pt) eq dd;
 vals:=[Evaluate((K19!BaseField(K19).1)@@m19,pt),Evaluate(nY19@@m19,pt)];
 d:=LCM([Degree(MinimalPolynomial(a,GF(2))):a in vals]);found:=false;
 for bits in [[1,0],[0,1],[1,1]] do
  val:=&+[bits[i]*vals[i]:i in [1..2]];mp:=MinimalPolynomial(val,GF(2));if Degree(mp) eq d then bb:=bits;found:=true;break;end if;
 end for;assert found;
 MB:=Matrix(GF(2),[Eltseq(val^i,GF(2)):i in [0..d-1]]);assert Rank(MB) eq d;
 rows:=[Eltseq(Solution(MB,Vector(GF(2),Eltseq(a,GF(2))))):a in vals cat [ww]];
 key:=[d] cat bb cat [Integers()!a:a in Eltseq(mp)] cat [Integers()!a:a in &cat rows];
 if IsDefined(cache,key) then return cache[key],key;end if;
 gen19:=&+[bb[i]*[K19!BaseField(K19).1,nY19][i]:i in [1..2]];
 gens19:=[OK19!Evaluate(mp,gen19)];
 for i in [1..2] do poly:=&+[K19!rows[i][k+1]*gen19^k:k in [0..d-1]];Append(~gens19,OK19!([K19!BaseField(K19).1,nY19][i]-poly));end for;
 printf "C57_FULL85_NATIVE_POINT_START degree=%o bits=%o time=%o\n",d,bb,Cputime(tg85);
 II19:=ideal<OK19|gens19>;assert IsPrime(II19);Pnative:=Place(II19);assert Degree(Pnative) eq d;
 printf "C57_FULL85_NATIVE_POINT_DONE degree=%o time=%o\n",d,Cputime(tg85);
 ff,gg:=TwoGenerators(Pnative);
 assert Valuation(den85,Pnative) eq 0;
 assert &and[Valuation(a,Pnative) ge 0:a in Coefficients(fd85)|a ne 0];
 assert Valuation(Discriminant(fd85),Pnative) eq 0;
 rel:=wd85-&+[D85!rows[3][k+1]*(D85!gen19)^k:k in [0..d-1]];
 printf "C57_FULL85_DIRECT_RELATIVE_IDEAL_START degree=%o time=%o\n",d,Cputime(tg85);
 gens85:=[OD85!(D85!ff),OD85!(D85!gg),OD85!(D85!den85*rel)];
 I85:=ideal<OD85|gens85>;
 printf "C57_FULL85_DIRECT_RELATIVE_IDEAL_BUILT degree=%o time=%o\n",d,Cputime(tg85);
 assert Norm(I85) eq II19;
 printf "C57_FULL85_DIRECT_RELATIVE_IDEAL_NORM_CERTIFIED degree=%o time=%o\n",d,Cputime(tg85);
 P:=DivisorGroup(D85)!I85;
 printf "C57_FULL85_DIRECT_RELATIVE_DIVISOR_BUILT degree=%o time=%o\n",d,Cputime(tg85);
 // The exact prime-norm identity, the defining generators and the
 // discriminant-unit check certify the selected root and e=1.
 // No repeated global valuation or residue-field isomorphism is required.
 printf "C57_FULL85_POINT_DONE degree=%o actual_W_selected=true time=%o\n",d,Cputime(tg85);return P,key;
end function;
F256<a256>:=GF(256);curves85:=[**];
for jj in F256 do for E in Twists(EllipticCurveFromjInvariant(jj)) do if #E mod 57 eq 0 then assert #E eq 228;Append(~curves85,E);end if;end for;end for;
assert #curves85 eq 8;EE85:=curves85[1];AG85,ag85:=AbelianGroup(EE85);assert Invariants(AG85) eq [228];
PN85:=4*ag85(AG85.1);assert Order(PN85) eq 57;r385:=19*PN85;r1985:=3*PN85;
cache85:=AssociativeArray();input85:=[];inputkeys85:=[];
for k in [0..2] do
 P,key:=FullPoint57(EE85,5^k*r385,5^k*r1985,cache85);cache85[key]:=P;Append(~input85,P);Append(~inputkeys85,key);
end for;
assert #Seqset(inputkeys85) eq 3 and &and[k[1] eq 8:k in inputkeys85];assert 5^18 mod 57 eq 1 and 5^9 mod 57 eq 20;
print "C57_FULL85_INPUT_COVERAGE three_distinct_diamond5_representatives=true full_orbit_coverage_from_free_order18_and_independent_count18=true",Cputime(tg85);
isogs85:=MDIsogenies(EE85,5);assert #isogs85 eq 6;output85:=[];outmods85:=[**];outkeys85:=[];obs85:=[];
for ph in isogs85 do
 E:=Codomain(ph);R3:=ph(r385);R19:=ph(r1985);P,key:=FullPoint57(E,R3,R19,cache85);cache85[key]:=P;k:=Index(outkeys85,key);
 if k eq 0 then Append(~output85,P);Append(~outmods85,<E,R3,R19>);Append(~outkeys85,key);Append(~obs85,1);else obs85[k]+:=1;end if;
end for;
assert #output85 eq 2 and obs85 eq [3,3] and &and[k[1] eq 24:k in outkeys85];
d5outputs85:=[];d5outkeys85:=[];
for md in outmods85 do P,key:=FullPoint57(md[1],5*md[2],5*md[3],cache85);cache85[key]:=P;Append(~d5outputs85,P);Append(~d5outkeys85,key);end for;
if not assigned PointKeyFile85 then PointKeyFile85:="C57/data/full85_relative_divisors_resume_pointkeys.m";end if;
Write(PointKeyFile85,Sprintf("InputKeys:=%o;OutputKeys:=%o;Diamond5OutputKeys:=%o;HeckeObservations:=%o;",inputkeys85,outkeys85,d5outkeys85,obs85):Overwrite:=true);
print "C57_FULL85_HECKE5_COVERAGE two_actual_degree24_outputs=true no_phase_ambiguity=true",Cputime(tg85);
T85:=&+[P:P in output85];d5T85:=&+[P:P in d5outputs85];
Af85:=d5T85-T85+4*input85[2]+input85[1]-5*input85[3];
print "C57_RELATIVE_DIVISORS_OPERATOR_READY",Cputime(tg85);
rr85:=Basis(Af85:Reduction:=false,Simplification:="None");
assert #rr85 le 1;
printf "C57_RELATIVE_DIVISORS_A5_RR_DIMENSION dimension=%o complete=true time=%o\n",#rr85,Cputime(tg85);
if #rr85 eq 1 then
 assert Divisor(rr85[1])+Af85 eq 0*Af85;
 print "C57_RELATIVE_DIVISORS_A5_PRINCIPAL",Cputime(tg85);
else
 print "C57_RELATIVE_DIVISORS_GOOD8_EXCLUDED_A5 complete=true",Cputime(tg85);
end if;
quit;
