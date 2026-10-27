/* Fresh octic X1(35)/F2 filter. No all-cusp enumeration is needed:
   at most two of the four non-Q cubic cusps fit in degree8, and their
   membership in the rational image is tested individually and in pairs. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(16*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/class_rational_reduction.m");
out35:=Open("C35/logs/classify_octic_35.log","w");
procedure Note35(s)
 fprintf out35,"%o\n",s;Flush(out35);print s;
end procedure;
Note35("C35_D8_START");
fp:=FrobeniusPolynomial(JOne(35),2);jf2:=Integers()!Evaluate(fp,1);
X:=MDX1(35,GF(2));C:=Curve(X);assert Genus(X) eq 25;
Note35(Sprintf("C35_D8_MODEL_NORMALIZED genus=25 time=%o",Cputime(t0)));
cs:=Cusps(X);osraw:=CuspOrbitsQ(X);
os:=[Type(o) eq PlcCrvElt select Divisor(o) else o:o in osraw];
assert {*Degree(c):c in cs*} eq {*1^^12,3^^4,4^^3,12*};
assert {*Degree(o):o in os*} eq {*1^^12,4^^3,6^^2,12*};
assert #Seqset([CuspSignature(X,Support(o)[1]):o in os]) eq 18;
assert &and[#Support(o) eq 2 and &and[Degree(c) eq 3:c in Support(o)]:o in os|Degree(o) eq 6];
c0:=Support([o:o in os|Degree(o) eq 1][1])[1];
Note35(Sprintf("C35_D8_CUSPS_READY finite=%o rational_orbits=%o JF2_order=%o time=%o",#cs,#os,jf2,Cputime(t0)));
G,m1,m2:=ClassGroupRationalReduction(C,c0:Proof:=true);
assert &*[n:n in Invariants(G)|n ne 0] eq jf2;
qc:=[m2(o-Degree(o)*Divisor(c0)):o in os];HQ:=sub<G|qc>;
assert Invariants(HQ) eq [13,109148520] and 13*#HQ eq jf2;
Note35(Sprintf("C35_D8_GROUP_READY invariants=%o rational_image=%o time=%o",Invariants(G),Invariants(HQ),Cputime(t0)));
places:=[NoncuspidalPlaces(X,d):d in [1..8]];
assert [#p:p in places] eq [0,0,0,0,0,24,24,30];
assert Degree(fp) eq 50;traces:=[];closed:=[];
for n in [1..8] do
 sn:=-n*Coefficient(fp,50-n);
 for j in [1..n-1] do sn-:=Coefficient(fp,50-j)*traces[n-j];end for;
 Append(~traces,Integers()!sn);count:=2^n+1-Integers()!sn;
 for d in Divisors(n) do if d lt n then count-:=d*closed[d];end if;end for;
 assert count mod n eq 0;Append(~closed,count div n);
 assert closed[n] eq #places[n]+#[c:c in cs|Degree(c) eq n];
end for;
Note35(Sprintf("C35_D8_COUNTS_VERIFIED noncuspidal=%o time=%o",[#p:p in places],Cputime(t0)));
bad:=[c:c in cs|Degree(c) eq 3];assert #bad eq 4;
badclasses:=[m2(Divisor(c)-3*Divisor(c0)):c in bad];
bad_pairs:=[];bad_complete:=true;
for i in [1..4] do
 inH:=badclasses[i] in HQ;if inH then bad_complete:=false;end if;
 Note35(Sprintf("C35_D8_CUSP_SINGLE index=%o in_rational_image=%o",i,inH));
 for j in [i..4] do
  expected:=i ne j and &or[Seqset(Support(o)) eq {bad[i],bad[j]}:o in os|Degree(o) eq 6];
  actual:=badclasses[i]+badclasses[j] in HQ;
  if actual ne expected then bad_complete:=false;end if;
  Append(~bad_pairs,[i,j,actual select 1 else 0,expected select 1 else 0]);
  Note35(Sprintf("C35_D8_CUSP_PAIR i=%o j=%o in_rational_image=%o rational_orbit_pair=%o",i,j,actual,expected));
 end for;
end for;
assert &and[m2(Divisor(c)-Degree(c)*Divisor(c0)) in HQ:c in cs|Degree(c) in [1,4]];
Note35(Sprintf("C35_D8_CUSP_CLASSIFICATION_DONE complete=%o",bad_complete));
load "realizations/data/vanhoeij_octic_input.m";
N8,fc,xc,yc:=Explode(OcticExamples[5]);assert N8 eq 35;
R<z>:=PolynomialRing(Rationals());f8:=R!fc;xp8:=R!xc;yp8:=R!yc;
assert IsIrreducible(f8) and Degree(f8) eq 8;
K8<a8>:=NumberField(f8);xx:=Evaluate(xp8,a8);yy:=Evaluate(yp8,a8);
r:=(xx^2*yy-xx*yy+yy-1)/(xx^2*yy-xx);s:=(xx*yy-yy+1)/(xx*yy);
b:=r*s*(r-1);c:=s*(r-1);E8:=EllipticCurve([1-c,-b,-b,0,0]);P8:=E8![0,0];
assert 35*P8 eq E8!0 and 7*P8 ne E8!0 and 5*P8 ne E8!0;
assert Degree(MinimalPolynomial(jInvariant(E8))) eq 4;
assert IsIrreducible(Polynomial(GF(2),f8));
L8<aa8>:=ext<GF(2)|Polynomial(GF(2),f8)>;
p8:=Places(C(L8)![Evaluate(Polynomial(GF(2),xp8),aa8),Evaluate(Polynomial(GF(2),yp8),aa8)]);
assert #p8 eq 1 and Degree(p8[1]) eq 8 and not IsCusp(X,p8[1]);
known8:=Seqset(DiamondOrbit(X,p8[1]));classes8:={m2(Divisor(p)-8*Divisor(c0)):p in known8};
assert #known8 eq #classes8;
Note35(Sprintf("C35_D8_KNOWN_OCTICS orbit=%o classes=%o time=%o",#known8,#classes8,Cputime(t0)));
arith:=Open("C35/data/class_arithmetic.py","w");
fprintf arith,"relations = %o\n",[Eltseq(r):r in Rows(RelationMatrix(G))];
fprintf arith,"rational_cusp_generators = %o\n",[Eltseq(a):a in qc];
fprintf arith,"bad_cubic_cusps = %o\n",[Eltseq(a):a in badclasses];
fprintf arith,"bad_cusp_pairs = %o\n",bad_pairs;
fprintf arith,"known_octic_classes = %o\n",[Eltseq(a):a in classes8];
fprintf arith,"frobenius_coefficients = %o\n",Coefficients(fp);
fprintf arith,"tests = [\n";
passes6:=0;passes7:=0;passes8:=0;candidates8:={G|};
for d in [6..8] do
 for i in [1..#places[d]] do
  a:=m2(Divisor(places[d][i])-d*Divisor(c0));inH:=a in HQ;
  if inH then
   if d eq 6 then passes6+:=1;elif d eq 7 then passes7+:=1;else passes8+:=1;Include(~candidates8,a);end if;
  end if;
  fprintf arith,"[%o,%o,%o],\n",d,Eltseq(a),inH select 1 else 0;
  Note35(Sprintf("C35_D8_CLASS_TEST degree=%o index=%o/%o in_rational_image=%o time=%o",d,i,#places[d],inH,Cputime(t0)));
 end for;
end for;
fprintf arith,"]\ncomplete = True\n";delete arith;
unknown8:=candidates8 diff classes8;assert classes8 subset candidates8;
complete35:=bad_complete and passes6 eq 0 and passes7 eq 0 and #unknown8 eq 0;
Note35(Sprintf("C35_OCTIC_CLASSIFICATION_DONE sextic_members=%o septic_members=%o octic_members=%o known_octic_classes=%o unknown_octic_classes=%o cusp_classification=%o complete=%o time=%o",passes6,passes7,passes8,#classes8,#unknown8,bad_complete,complete35,Cputime(t0)));
if complete35 then delete out35;quit;end if;
