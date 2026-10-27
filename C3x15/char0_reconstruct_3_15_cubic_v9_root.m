/* Independent lower-coefficient model: cyclic cubic over X1(15).
   Characteristic-zero degree4 RR, exact signed global cusp classes. */
SetColumns(0);SetQuitOnError(true);SetLogFile("C3x15/logs/char0_reconstruct_3_15_cubic_v9_root.log");SetMemoryLimit(16*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
data:=eval Read("C3x15/data/M3N15_D4_exact_cusp_lifts.txt");
mm,nn,dd,sigs0,ds0,baseidx,coeffs,targets,pairs:=Explode(data);
F<z>:=CyclotomicField(3);X:=MDX1(15,F);FF:=FunctionField(Curve(X));A,afmap:=AlgorithmicFunctionField(FF);
E0:=EllipticCurve([afmap(FF!f):f in X`_E]);Q0:=E0![afmap(FF!f):f in X`_P];assert 15*Q0 eq E0!0 and 5*Q0 ne E0!0 and 3*Q0 ne E0!0;
print "M3N15_CUBIC_V9_BASE_READY",Genus(A),Cputime(tt);
base_cs:=Poles(jInvariant(E0));print "M3N15_CUBIC_V9_BASE_CUSPS",[Degree(c):c in base_cs],Cputime(tt);
R0:=5*Q0;aa1,aa2,aa3,aa4,aa6:=Explode(aInvariants(E0));
rx:=R0[1];ry:=R0[2];lm:=(3*rx^2+2*aa2*rx+aa4-aa1*ry)/(2*ry+aa1*rx+aa3);
assert aa2-lm*aa1+3*rx-lm^2 eq 0;
assert aa4-lm*aa3+2*aa2*rx-(ry+lm*rx)*aa1+3*rx^2-2*lm*ry eq 0;
Ebase:=EllipticCurve([aa1+2*lm,A!0,2*ry+aa1*rx+aa3,A!0,A!0]);
Qbase:=Ebase![Q0[1]-rx,Q0[2]-ry-lm*(Q0[1]-rx)];assert 15*Qbase eq Ebase!0;
assert jInvariant(Ebase) eq jInvariant(E0);
PZ:=Parent(DivisionPolynomial(Ebase,3));pol:=ExactQuotient(DivisionPolynomial(Ebase,3),PZ.1);pol/:=LeadingCoefficient(pol);
print "M3N15_CUBIC_V9_POLYNOMIAL_READY",Degree(pol),Cputime(tt);
assert Degree(pol) eq 3 and IsIrreducible(pol);
print "M3N15_CUBIC_V9_IRREDUCIBLE",Cputime(tt);
K<r>:=ext<A|pol>;E:=ChangeRing(Ebase,K);Q:=E!Qbase;
a1,a2,a3,a4,a6:=Explode(aInvariants(E));
yp:=-3*z^2*r^3/((1+2*z)*a1*r+3*z*a3);
assert yp^2+a1*r*yp+a3*yp eq r^3;
P:=E![r,yp];assert 3*P eq E!0 and P ne E!0;
e:=WeilPairing(P,5*Q,3);assert e in {K!z,K!z^2};if e ne z^2 then P:=-P;end if;assert WeilPairing(P,5*Q,3) eq z^2;
print "M3N15_CUBIC_V9_FIELD",Genus(K),Degree(K,A),Cputime(tt);assert Genus(K) eq 9;
x0:=P[1];y0:=P[2];sl:=(3*x0^2+2*a2*x0+a4-a1*y0)/(2*y0+a1*x0+a3);
assert a2-sl*a1+3*x0-sl^2 eq 0;
assert a4-sl*a3+2*a2*x0-(y0+sl*x0)*a1+3*x0^2-2*sl*y0 eq 0;
xq:=Q[1]-x0;yq:=Q[2]-y0-sl*xq;sc:=yq/xq^2;
A1:=(a1+2*sl)*sc;A3:=(2*y0+a1*x0+a3)*sc^3;
E:=EllipticCurve([A1,K!0,A3,K!0,K!0]);P:=E![0,0];Q:=E![sc^2*xq,sc^3*yq];
assert 15*Q eq E!0 and WeilPairing(P,5*Q,3) eq z^2;
t:=-Q[1];assert Q[2] eq t^2;q:=(A1-t-2)/(t-1);assert A3 eq q*t^2+(1-q)*t;
fs_short:=[A1^3-27*A3,(Q+2*P)[1],Q[1]-(Q+P)[1]];
print "M3N15_CUBIC_V9_UNITS_READY",#fs_short,Cputime(tt);
j:=K!jInvariant(E0);assert j eq jInvariant(E);cs:=[];
for ci in [1..#base_cs] do
 cp:=base_cs[ci];if Degree(cp) ne 1 or (-Valuation(jInvariant(E0),cp)) mod 3 ne 0 then continue;end if;fcp,gcp:=TwoGenerators(cp);ucp:=UniformizingElement(cp);
 print "M3N15_CUBIC_V9_CUSP_FIBRE_START",ci,Degree(cp),Cputime(tt);
 raw:=Decomposition(K,cp);
 above:=[qq:qq in raw|Valuation(K!fcp,qq) gt 0 and Valuation(K!gcp,qq) gt 0];
 assert #above gt 0 and &+[Valuation(K!ucp,qq)*Degree(qq):qq in above] eq 3*Degree(cp);
 assert &and[Valuation(j,qq) lt 0:qq in above];cs cat:=above;
 print "M3N15_CUBIC_V9_CUSP_FIBRE_DONE",ci,[Degree(qq):qq in above],Cputime(tt);
end for;
assert #cs eq 12 and &and[Degree(qq) eq 1:qq in cs];
print "M3N15_CUBIC_V9_CUSPS",[Degree(c):c in cs],Cputime(tt);
blocks:=AssociativeArray();
for c in cs do
 print "M3N15_CUBIC_V9_SIGNATURE_START",Index(cs,c),Cputime(tt);
 sel:=[2,7,9];reduced_sigs:=[[v[k]:k in sel]:v in sigs0];assert #Seqset(reduced_sigs) eq 20;
 sig:=[Valuation(f,c):f in fs_short];i:=Index(reduced_sigs,sig);assert i gt 0 and ds0[i] eq 1;assert not IsDefined(blocks,i);blocks[i]:=Divisor(c);
 print "M3N15_CUBIC_V9_SIGNATURE_MATCH",i,Cputime(tt);
end for;
assert &and[IsDefined(blocks,i):i in [1,7,9,12,17,18]];
c0:=blocks[baseidx];print "M3N15_CUBIC_V9_MATCHED_NEEDED_RATIONAL_BLOCKS",Cputime(tt);
// All12 classes form one F-defined automorphism orbit; test short lift4 first.
for i in [4] do
 v:=coeffs[i];assert &and[ds0[j] eq 1:j in [1..#ds0]|v[j] ne 0];
 L:=4*c0+&+[v[j]*(blocks[j]-ds0[j]*c0):j in [1..#ds0]|v[j] ne 0];assert Degree(L) eq 4;
 print "M3N15_CUBIC_V9_RR_START",i,Cputime(tt);
 rr_basis:=Basis(L:Reduction:=false,Simplification:="None");h:=#rr_basis;assert h le 1;print "M3N15_CUBIC_V9_RR_DIMENSION",i,h,Cputime(tt);
 assert h eq 1;
 function EncodeBase(g)
  return [* [Eltseq(c):c in Coefficients(Numerator(g))], [Eltseq(c):c in Coefficients(Denominator(g))] *];
 end function;
 function EncodeA(g)
  return [*EncodeBase(c):c in Eltseq(A!g)*];
 end function;
 function EncodeK(g)
  return [*EncodeA(c):c in Eltseq(K!g)*];
 end function;
 witness:=[* [*EncodeBase(c):c in Coefficients(DefiningPolynomial(A))*], [*EncodeA(c):c in Coefficients(pol)*], EncodeK(rr_basis[1]), [*EncodeK(c):c in aInvariants(E)*], [*EncodeK(P[j]):j in [1,2]*], [*EncodeK(Q[j]):j in [1,2]*], EncodeK(q), EncodeK(t) *];
 Write("C3x15/data/M3N15_char0_RR_encoded_witness_root.m",Sprint(witness):Overwrite:=true);
 print "M3N15_CUBIC_V9_RR_WITNESS_EXPORTED",Cputime(tt);
end for;
print "M3N15_CUBIC_V9_CHARZERO_OCCURRENCE_RR_COMPLETE",Cputime(tt);quit;
