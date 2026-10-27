/* Exact independent reconstruction from the characteristic-zero place's
   two generators. No cusp construction or Riemann--Roch is needed here. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
R<t>:=PolynomialRing(Rationals());
f:=t^8-39/8*t^7+149/16*t^6-59/8*t^5-11/8*t^4+33/4*t^3-15/2*t^2+3*t-1/2;
assert IsIrreducible(f);K<z>:=NumberField(f);assert Degree(K) eq 8;
alpha:=(8/1355*z^7+221/271*z^6-9021/2710*z^5+13389/2710*z^4-272/271*z^3-4714/1355*z^2+4693/1355*z-2312/1355)/(z^2-z);
Y:=MDX1(28,Rationals());FF:=FunctionField(Curve(Y));AF,am:=AlgorithmicFunctionField(FF);
function EvalBase(r,z)
 return Evaluate(Numerator(r),z)/Evaluate(Denominator(r),z);
end function;
function EvalAF(g,z,alpha)
 v:=Eltseq(g);return &+[EvalBase(v[i],z)*alpha^(i-1):i in [1..#v]];
end function;
df:=DefiningPolynomial(AF);assert &+[EvalBase(Coefficient(df,i),z)*alpha^i:i in [0..Degree(df)]] eq 0;
ai:=[EvalAF(am(FF!a),z,alpha):a in Y`_E];
E:=EllipticCurve(ai);P:=E![EvalAF(am(FF!a),z,alpha):a in Y`_P];
assert Discriminant(E) ne 0 and 28*P eq E!0 and 14*P ne E!0 and 4*P ne E!0;
S<x>:=PolynomialRing(K);a1,a2,a3,a4,a6:=Explode(ai);
psi2:=4*x^3+(a1^2+4*a2)*x^2+(2*a1*a3+4*a4)*x+a3^2+4*a6;
fac:=Factorization(psi2);roots:=Roots(psi2);assert &+[Degree(h[1])*h[2]:h in fac] eq 3;
pts:=[];for rr in roots do xx:=rr[1];yy:=-(a1*xx+a3)/2;PP:=E![xx,yy];assert PP ne E!0 and 2*PP eq E!0;Append(~pts,PP);end for;
print "N28_D8_QUOTIENT_EXACT_FIELD",Coefficients(f);
print "N28_D8_QUOTIENT_EXACT_ORDER28",true,Cputime(t0);
print "N28_D8_QUOTIENT_PSI2_FACTOR_DEGREES",[Degree(h[1]):h in fac],"rational_2_points",#pts,Cputime(t0);
print "N28_D8_QUOTIENT_EXACT_A_INVARIANTS",ai;
print "N28_D8_QUOTIENT_EXACT_P28",P;
Write("C2x28/data/N28_D8_quotient_direct_torsion_arithmetic.txt",Sprint([*Coefficients(f),Eltseq(alpha),[Eltseq(a):a in ai],[Eltseq(P[1]),Eltseq(P[2])],[Degree(h[1]):h in fac],#pts*]):Overwrite:=true);
print "N28_D8_QUOTIENT_DIRECT_TORSION_DONE",Cputime(t0);quit;
