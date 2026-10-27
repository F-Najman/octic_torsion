/* New canonical degree20 quotient model attempt; no exclusion claim. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(12*10^9);tt:=Cputime();
R<T,x,y>:=PolynomialRing(GF(2),3);B:=FieldOfFractions(R);
F25:=eval Read("mdmagma/models_X1_n/FFFc25.txt");
xx:=B!x;yy:=B!y;r:=(xx^2*yy-xx*yy+yy-1)/(xx*(xx*yy-1));s:=(xx*yy-yy+1)/(xx*yy);
b:=r*s*(r-1);c:=s*(r-1);
// The exact degree5 modular unit from the archived vanHoeij gonality table.
f3:=b;f4:=c;f5:=b-c;f6:=s-1;f7:=s-r;f8:=r*s-2*r+1;f9:=s^2-s-r+1;
f10:=xx-yy+1;f12:=xx-yy;f13:=xx^3*yy-xx^2*yy^2-xx^2*yy+xx*yy^2-yy+1;
h:=f3*f4*f6^3*f7^2*f10*f12^3/(f5*f8^2*f9^2*f13^3);
N:=R!Numerator(h);D:=R!Denominator(h);
print "N75_GON5_CANONICAL_RESULTANT_START",TotalDegree(N),TotalDegree(D),Cputime(tt);
res:=Resultant(F25,N-T*D,3);assert res ne 0;
K<t>:=RationalFunctionField(GF(2));S<X>:=PolynomialRing(K);proj:=hom<R->S|[S!t,X,S!0]>;
pol:=proj(res);factors:=Factorization(pol);
print "N75_GON5_CANONICAL_RESULTANT_FACTORS",[<Degree(z[1]),z[2]>:z in factors],Cputime(tt);
fs:=[z[1]:z in factors|Degree(z[1]) eq 5];assert #fs eq 1;
L<a>:=FunctionField(fs[1]);SY<Y>:=PolynomialRing(L);ev:=hom<R->SY|[SY!t,SY!a,Y]>;
gg:=GCD(ev(F25),ev(N-T*D));assert Degree(gg) eq 1;
ya:=-Coefficient(gg,0)/Coefficient(gg,1);assert Evaluate(ev(F25),ya) eq 0;
xa:=a;ra:=(xa^2*ya-xa*ya+ya-1)/(xa*(xa*ya-1));sa:=(xa*ya-ya+1)/(xa*ya);
ba:=ra*sa*(ra-1);ca:=sa*(ra-1);E:=EllipticCurve([1-ca,-ba,-ba,L!0,L!0]);P:=E![0,0];
assert 25*P eq E!0 and 5*P ne E!0 and Genus(L) eq 12;
print "N75_GON5_CANONICAL_BASE_MODEL_VERIFIED",Degree(L),Genus(L),Cputime(tt);
S3<W3>:=PolynomialRing(L);j25:=jInvariant(E);f:=W3^4+W3^3+1/j25;assert Degree(f) eq 4;
assert IsIrreducible(f);M<w>:=FunctionField(f);
print "N75_GON5_CANONICAL_DEGREE20_GENUS_START",Cputime(tt);
assert Genus(M) eq 73;
print "N75_GON5_CANONICAL_GENUS73_MODEL_VERIFIED",Cputime(tt);
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke_11.m");
coords:=[M!t,M!a,w];
Attach("common/fast_hecke.m");
load "common/pair_af_place_linear.m";
OL75:=MaximalOrderFinite(L);OM75:=MaximalOrderFinite(M);
n75,d75:=IntegralSplit(1/j25,OL75);assert M!d75*w in OM75;
function NativePoint75Relative(vals,wv)
 p0:=C48PairAFPlaceLinear([L!t,a],vals);d:=Degree(p0);
 assert d eq LCM([Degree(MinimalPolynomial(v,GF(2))):v in vals]);
 assert Valuation(d75,p0) eq 0 and Valuation(j25,p0) eq 0;
 assert &and[Valuation(c,p0) ge 0:c in Coefficients(f)|c ne 0];
 assert Valuation(Discriminant(f),p0) eq 0;
 F:=Parent(wv);found:=false;
 for bits in [[1,0],[0,1],[1,1]] do
  theta:=&+[bits[i]*vals[i]:i in [1..2]];
  if Degree(MinimalPolynomial(theta,GF(2))) eq d then chosen:=bits;found:=true;break;end if;
 end for;assert found;
 mat:=Matrix(GF(2),[Eltseq(theta^i,GF(2)):i in [0..d-1]]);
 cs:=Eltseq(Solution(mat,Vector(GF(2),Eltseq(wv,GF(2)))));
 thetaf:=chosen[1]*(L!t)+chosen[2]*a;
 rel:=w-&+[M!cs[i+1]*(M!thetaf)^i:i in [0..d-1]];
 ff,gg:=TwoGenerators(p0);
 print "N75_RELATIVE_POINT_START",d,Cputime(tt);
 I:=ideal<OM75|[OM75!(M!ff),OM75!(M!gg),OM75!(M!d75*rel)]>;
 print "N75_RELATIVE_IDEAL_READY",d,Cputime(tt);
 assert IsPrime(I);pp:=Place(I);assert Degree(pp) eq d;
 assert Valuation(M!ff,pp) gt 0 and Valuation(M!gg,pp) gt 0 and Valuation(rel,pp) gt 0;
 assert Valuation(M!UniformizingElement(p0),pp) eq 1;
 print "N75_RELATIVE_POINT_DONE",d,Cputime(tt);return pp;
end function;

function AbsoluteRecord(vals)
 F:=Universe(vals);d:=LCM([Degree(MinimalPolynomial(v,GF(2))):v in vals]);assert d in [7,14];
 found:=false;
 for v in vals do if Degree(MinimalPolynomial(v,GF(2))) eq d then theta:=v;found:=true;break;end if;end for;
 assert found;mp:=MinimalPolynomial(theta,GF(2));
 mat:=Matrix(GF(2),[Eltseq(theta^j,GF(2)):j in [0..d-1]]);assert Rank(mat) eq d;
 enc:=[];
 for v in vals do cc:=Eltseq(Solution(mat,Vector(GF(2),Eltseq(v,GF(2)))));Append(~enc,&+[Integers()!cc[j+1]*2^j:j in [0..d-1]]);end for;
 return [Integers()!c:c in Coefficients(mp)],enc;
end function;
rowsfile:="C75/data/C75_good7_moduli_records_root_v2.txt";
fo:=Open(rowsfile,"w");delete fo;
function ModuliY(E,P25,P3,K,coords)
 assert Order(P25) eq 25 and Order(P3) eq 3;
 a1,a2,a3,a4,a6:=Explode(aInvariants(E));x0:=P25[1];y0:=P25[2];
 aa3:=2*y0+a3+a1*x0;aa4:=3*x0^2+2*x0*a2+a4-a1*y0;sl:=aa4/aa3;
 aa2:=3*x0+a2-sl*a1-sl^2;sc:=aa2/aa3;
 b:=-aa2^3/aa3^2;c:=1-(a1+2*sl)*sc;
 EC:=EllipticCurve([1-c,-b,-b,BaseRing(E)!0,BaseRing(E)!0]);
 RR:=EC![sc^2*(P3[1]-x0),sc^3*(P3[2]-y0-sl*(P3[1]-x0))];assert Order(RR) eq 3;
 aa:=aInvariants(EC);a1,a2,a3,a4,a6:=Explode(aa);
 xc:=(RR[1]-a3/a1)/a1^2;assert xc^4+xc^3+1/jInvariant(E) eq 0;
 xy:=TateNormalForm_xy(E,P25);xx,yy:=Explode(xy);
 r:=(xx^2*yy-xx*yy+yy-1)/(xx*(xx*yy-1));s:=(xx*yy-yy+1)/(xx*yy);
 f3:=b;f4:=c;f5:=b-c;f6:=s-1;f7:=s-r;f8:=r*s-2*r+1;f9:=s^2-s-r+1;
 f10:=xx-yy+1;f12:=xx-yy;f13:=xx^3*yy-xx^2*yy^2-xx^2*yy+xx*yy^2-yy+1;
 hv:=f3*f4*f6^3*f7^2*f10*f12^3/(f5*f8^2*f9^2*f13^3);
 pp:=NativePoint75Relative([hv,xx],xc);
 poly,enc:=AbsoluteRecord(aInvariants(E) cat [P25[1],P25[2],P3[1],P3[2],b,c,xx,yy,hv,xc]);
 row:=[*poly,enc,Degree(pp)*];
 Write(rowsfile,Sprint(row));return pp;
end function;
F:=GF(128);curves:=[* *];
for j in F do
 for E in Twists(EllipticCurveFromjInvariant(j)) do
  inv:=Invariants(AbelianGroup(E));
  if #inv gt 0 and inv[#inv] mod 75 eq 0 then assert #E eq 150 and inv eq [150];Append(~curves,E);end if;
 end for;
end for;
assert #curves eq 7;E:=curves[1];AG,am:=AbelianGroup(E);P75:=2*am(AG.1);P25:=3*P75;P3:=25*P75;
assert Order(P75) eq 75;P:=ModuliY(E,P25,P3,M,coords);assert Degree(P) eq 7;
print "N75_GON5_CANONICAL_INPUT_POINT",Degree(P),Cputime(tt);
for q in [7] do
 print "N75_GON5_CANONICAL_HECKE_START",q,Cputime(tt);
 terms:=MDIsogeniesFastPoints(E,[P25,P3],q);H:=DivisorGroup(M)!0;
 for term in terms do
  EG,points,e:=Explode(term);pp:=ModuliY(EG,points[1],points[2],M,coords);
  m:=e*7/Degree(pp);assert Denominator(m) eq 1;H+:=(Integers()!m)*Divisor(pp);
  print "N75_GON5_CANONICAL_HECKE_TERM",q,e,Degree(pp),Cputime(tt);
 end for;
 assert Degree(H) eq (q+1)*7;delta:=ModuliY(E,q*P25,q*P3,M,coords);assert Degree(delta) eq 7;
 D:=H-Divisor(P)-q*Divisor(delta);assert Degree(D) eq 0;
 print "N75_GON5_CANONICAL_PRINCIPAL_START",q,Cputime(tt);
 yes:=IsPrincipal(D);print "N75_GON5_CANONICAL_PRINCIPAL_RESULT",q,yes,Cputime(tt);
 assert not yes;
 assert #Basis(D:Reduction:=false,Simplification:="None") eq 0;
 ND:=DivisorGroup(M)!0;native:=MDIsogenies(E,q);assert #native eq q+1;
 for phi in native do
  EE:=Domain(phi);p25:=phi(EE![BaseRing(EE)!a:a in Eltseq(P25)]);p3:=phi(EE![BaseRing(EE)!a:a in Eltseq(P3)]);
  yy:=ModuliY(Codomain(phi),p25,p3,M,coords);ND+:=Divisor(yy);
 end for;
 assert ND eq 2*H;
 print "C75_GOOD7_NATIVE_REFERENCE_AND_RR_PASS",Degree(ND),Cputime(tt);
 print "C75_GOOD7_COMPLETE_EXCLUSION",q,Cputime(tt);
end for;
print "C75_GOOD7_CERTIFICATE_DONE complete=true",Cputime(tt);quit;
