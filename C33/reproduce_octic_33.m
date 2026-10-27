/* Fresh complete degree8 Picard-class enumeration on X1(33)/F2.
   No earlier-degree exclusion is used. Global septic places are merely
   freshly reconstructed effective degree7+cusp witnesses. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(16*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/class_rational_reduction.m");
if not assigned C33Log then C33Log:="C33/logs/root_fresh.log";end if;
out33:=Open(C33Log,"w");
procedure Note33(s)
 fprintf out33,"%o\n",s;Flush(out33);print s;
end procedure;
fp:=FrobeniusPolynomial(JOne(33),2);jf2:=Integers()!Evaluate(fp,1);
Note33(Sprintf("C33_D8_FROBENIUS_ORDER order=%o time=%o",jf2,Cputime(t0)));
X:=MDX1(33,GF(2));C:=Curve(X);cs:=Cusps(X);os:=CuspOrbitsQ(X);
assert Genus(X) eq 21;
assert Sort([Degree(o):o in os]) eq [1:i in [1..10]] cat [2:i in [1..5]] cat [10,10];
assert {*Degree(c):c in cs*} eq {*1^^10,2^^5,5^^2,10*};
assert #Seqset([CuspSignature(X,Support(o)[1]):o in os]) eq 17;
c0:=Support([o:o in os|Degree(o) eq 1][1])[1];
Note33(Sprintf("C33_D8_MODEL_READY genus=21 Q_cusp_orbits=17 finite_cusps=%o time=%o",#cs,Cputime(t0)));
G,m1,m2:=ClassGroupRationalReduction(C,c0:Proof:=true);
assert &*[n:n in Invariants(G)|n ne 0] eq jf2;
qc:=[m2(o-Degree(o)*Divisor(c0)):o in os];HQ:=sub<G|qc>;
assert Invariants(G) eq [15,42373650,0] and Invariants(HQ) eq [5,42373650];
Note33(Sprintf("C33_D8_GROUP_READY invariants=%o rational_image=%o time=%o",Invariants(G),Invariants(HQ),Cputime(t0)));
places:=[NoncuspidalPlaces(X,d):d in [1..8]];
assert [#p:p in places[1..7]] eq [0,0,0,0,2,20,30];
assert Degree(fp) eq 42;traces:=[];closed:=[];
for n in [1..8] do
 sn:=-n*Coefficient(fp,42-n);
 for j in [1..n-1] do sn-:=Coefficient(fp,42-j)*traces[n-j];end for;
 Append(~traces,Integers()!sn);count:=2^n+1-Integers()!sn;
 for d in Divisors(n) do if d lt n then count-:=d*closed[d];end if;end for;
 assert count mod n eq 0;Append(~closed,count div n);
 assert closed[n] eq #places[n]+#[c:c in cs|Degree(c) eq n];
end for;
Note33(Sprintf("C33_D8_COUNTS_VERIFIED noncuspidal=%o all_places=%o time=%o",[#p:p in places],closed,Cputime(t0)));

/* Fresh global witness of degree7; used only after adding a rational cusp. */
R<z>:=PolynomialRing(Rationals());
f7:=z^7-6*z^6+9*z^5-5*z^3-z+1;
y7:=(-3*z^6+16*z^5-17*z^4-9*z^3+11*z^2+3*z+4)/3;
assert IsIrreducible(f7);K7<a7>:=NumberField(f7);
function C33E(xx,yy)
 r:=(xx^2*yy-xx*yy+yy-1)/(xx^2*yy-xx);s:=(xx*yy-yy+1)/(xx*yy);
 b:=r*s*(r-1);c:=s*(r-1);return EllipticCurve([1-c,-b,-b,0,0]);
end function;
E7:=C33E(a7,Evaluate(y7,a7));P7:=E7![0,0];
assert 33*P7 eq E7!0 and 11*P7 ne E7!0 and 3*P7 ne E7!0;
assert Degree(MinimalPolynomial(jInvariant(E7))) eq 7;
assert IsIrreducible(Polynomial(GF(2),f7));
L7<aa7>:=ext<GF(2)|Polynomial(GF(2),f7)>;
yy7:=Evaluate(Polynomial(GF(2),y7),aa7);
p7:=Places(C(L7)![aa7,yy7]);assert #p7 eq 1 and Degree(p7[1]) eq 7 and not IsCusp(X,p7[1]);
known7:=Seqset(DiamondOrbit(X,p7[1]));assert #known7 eq 10;
classes7:={m2(Divisor(p)-7*Divisor(c0)):p in known7};assert #classes7 eq 10;
rcs:=[c:c in cs|Degree(c) eq 1];
rcc:=[m2(Divisor(c)-Divisor(c0)):c in rcs];
known71:={a+b:a in classes7,b in rcc};assert #known71 eq 100;

/* The octic model is independently certified by the degree8 realization
   package. We verify its global order and exact reduction again here. */
load "realizations/data/vanhoeij_octic_input.m";
N8,fc,xc,yc:=Explode(OcticExamples[1]);assert N8 eq 33;
f8:=R!fc;xp8:=R!xc;yp8:=R!yc;assert IsIrreducible(f8) and Degree(f8) eq 8;
K8<a8>:=NumberField(f8);E8:=C33E(Evaluate(xp8,a8),Evaluate(yp8,a8));P8:=E8![0,0];
assert 33*P8 eq E8!0 and 11*P8 ne E8!0 and 3*P8 ne E8!0;
assert Degree(MinimalPolynomial(jInvariant(E8))) eq 4;
assert IsIrreducible(Polynomial(GF(2),f8));
L8<aa8>:=ext<GF(2)|Polynomial(GF(2),f8)>;
p8:=Places(C(L8)![Evaluate(Polynomial(GF(2),xp8),aa8),Evaluate(Polynomial(GF(2),yp8),aa8)]);
assert #p8 eq 1 and Degree(p8[1]) eq 8 and not IsCusp(X,p8[1]);
known8:=Seqset(DiamondOrbit(X,p8[1]));
classes8:={m2(Divisor(p)-8*Divisor(c0)):p in known8};
assert #classes8 eq #known8;
Note33(Sprintf("C33_D8_KNOWN_CLASSES septic_plus_cusp=%o octic_orbit=%o octic_classes=%o time=%o",#known71,#known8,#classes8,Cputime(t0)));

cc:=[m2(Divisor(c)-Degree(c)*Divisor(c0)):c in cs];
procedure C33Effective(weights,classes,remaining,start,current,~out,~number)
 if remaining eq 0 then Include(~out,current);number+:=1;return;end if;
 for i in [start..#weights] do
  if weights[i] le remaining then
   C33Effective(weights,classes,remaining-weights[i],i,current+classes[i],~out,~number);
  end if;
 end for;
end procedure;
qeffective:={G|};number:=0;
C33Effective([Degree(o):o in os],qc,8,1,G!0,~qeffective,~number);
assert #qeffective eq number;
Note33(Sprintf("C33_D8_Q_EFFECTIVE divisors=%o classes=%o time=%o",number,#qeffective,Cputime(t0)));
allcusps:={G|};number:=0;
C33Effective([Degree(c):c in cs],cc,8,1,G!0,~allcusps,~number);
Note33(Sprintf("C33_D8_ALL_CUSPS divisors=%o classes=%o time=%o",number,#allcusps,Cputime(t0)));
complements:=AssociativeArray();
for k in [0..3] do
 out:={G|};number:=0;C33Effective([Degree(c):c in cs],cc,k,1,G!0,~out,~number);
 complements[k]:=out;
 Note33(Sprintf("C33_D8_COMPLEMENTS degree=%o divisors=%o classes=%o",k,number,#out));
end for;

if not assigned C33Arithmetic then C33Arithmetic:="C33/data/class_arithmetic.py";end if;
arith:=Open(C33Arithmetic,"w");
fprintf arith,"relations = %o\n",[Eltseq(r):r in Rows(RelationMatrix(G))];
fprintf arith,"rational_cusp_generators = %o\n",[Eltseq(a):a in qc];
fprintf arith,"rational_effective_cusps = %o\n",[Eltseq(a):a in qeffective];
fprintf arith,"known_septic_plus_cusp = %o\n",[Eltseq(a):a in known71];
fprintf arith,"known_octic_classes = %o\n",[Eltseq(a):a in classes8];
fprintf arith,"tests = [\n";
candidates:={G|};typecounts:=[];
procedure Test33(a,degree,~total,~passes,~qpasses,~candidates)
 inH:=a in HQ;inQ:=a in qeffective;
 total+:=1;if inH then passes+:=1;end if;
 if inH and inQ then qpasses+:=1;end if;
 if inH and not inQ then Include(~candidates,a);end if;
 fprintf arith,"[%o,%o,%o,%o],\n",degree,Eltseq(a),inH select 1 else 0,inQ select 1 else 0;
end procedure;
total:=0;passes:=0;qpasses:=0;
for a in allcusps do Test33(a,0,~total,~passes,~qpasses,~candidates);end for;
Append(~typecounts,<0,total,passes,qpasses>);
Note33(Sprintf("C33_D8_TEST_TYPE noncuspidal_degree=0 total=%o in_rational_image=%o rational_effective=%o candidates=%o time=%o",total,passes,qpasses,#candidates,Cputime(t0)));
for d in [5..8] do
 total:=0;passes:=0;qpasses:=0;
 for i in [1..#places[d]] do
  a:=m2(Divisor(places[d][i])-d*Divisor(c0));
  for b in complements[8-d] do Test33(a+b,d,~total,~passes,~qpasses,~candidates);end for;
  if i mod 20 eq 0 then Note33(Sprintf("C33_D8_TYPE_PROGRESS degree=%o place=%o/%o time=%o",d,i,#places[d],Cputime(t0)));end if;
 end for;
 Append(~typecounts,<d,total,passes,qpasses>);
 Note33(Sprintf("C33_D8_TEST_TYPE noncuspidal_degree=%o total=%o in_rational_image=%o rational_effective=%o candidates=%o time=%o",d,total,passes,qpasses,#candidates,Cputime(t0)));
end for;
fprintf arith,"]\ncomplete = True\n";delete arith;
unknown:=candidates diff (known71 join classes8);
assert known71 subset candidates and classes8 subset candidates;
Note33(Sprintf("C33_OCTIC_CLASSIFICATION_DONE candidates=%o septic_plus_cusp=%o octic_classes=%o unknown=%o complete=%o time=%o",#candidates,#known71,#classes8,#unknown,#unknown eq 0,Cputime(t0)));
Note33(Sprintf("C33_D8_TYPE_COUNTS %o",typecounts));
if #unknown eq 0 then delete out33;quit;end if;
