SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tn57:=Cputime();
try
CQ:=SmallModularCurve(19);CC:=ChangeRing(CQ,GF(2));assert IsNonsingular(CC) and #CC eq 3;
RC:=CoordinateRing(Ambient(CC));jr:=jInvariant(CQ,19);jn:=RC!Numerator(jr);jd:=RC!Denominator(jr);
K,km:=AlgorithmicFunctionField(FunctionField(CC));jK:=km(FunctionField(CC)!(jn/jd));
assert km(FunctionField(CC)!(RC.1/RC.3)) eq K!BaseField(K).1 and km(FunctionField(CC)!(RC.2/RC.3)) eq K.1;
R<v>:=PolynomialRing(K);AF<t3>:=FunctionField(v^4+v^3+1/jK);assert Genus(AF) eq 5;
// Exact AS coefficient over AF.
x0:=AF!K!BaseField(K).1;y0:=AF!K.1;u0:=AF.1;
aAS:=((x0^11 + x0^10 + x0^9 + x0^2 + x0)/(x0^11 + x0^10 + x0^9 + x0^8 + x0^3 + x0^2 + x0 + 1)*y0 + (x0^12 + x0^11 + x0^8 + x0^6 + x0^5 + x0^3 + 1)/(x0^8 + 1))*u0^3 + ((x0^13 + x0^12 + x0^10 + x0^9 + x0^8 + x0^7 + x0^6 + x0^5 + x0^4 + x0 + 1)/(x0^11 + x0^10 + x0^9 + x0^8 + x0^3 + x0^2 + x0 + 1)*y0 + (x0^11 + x0^8 + x0^7 + x0^5 + x0^4 + x0^2 + 1)/(x0^8 + 1))*u0^2 + ((x0^8 + x0^6 + x0^5 + x0 + 1)/(x0^6 + x0^4 + x0^2 + 1)*y0 + (x0^5 + x0^4 + x0^2 + x0 + 1)/(x0^3 + x0^2 + x0 + 1))*u0 + (x0^5 + x0^3 + x0)/(x0^3 + x0^2 + x0 + 1)*y0 + (x0^6 + x0^4 + x0)/(x0^3 + x0^2 + x0 + 1);assert Parent(aAS) cmpeq AF;

RY<z>:=PolynomialRing(AF);CY9<zz>:=FunctionField(z^2+z+aAS);assert Genus(CY9) eq 9;
cs9:=[P:P in Poles(CY9!AF!jK)|-Valuation(CY9!AF!jK,P) eq 57];assert #cs9 eq 2 and &and[Degree(P) eq 1:P in cs9];
Gfresh,_,phiFresh:=ClassGroup(CY9:ReductionDivisor:=Divisor(cs9[1]),Proof:=true);
assert Invariants(Gfresh) eq [105,210,0];print "C57_NORM9_FRESH_CLASSGROUP_PROVED",Invariants(Gfresh),Cputime(tn57);
Norm9Invariants:=[ 105, 210, 0 ];Norm9Relations:=[
[ 105, 0, 0 ],
[ 0, 210, 0 ]
];Norm9Classes:=AssociativeArray();
Norm9Classes[[ 14, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1 ]]:=[ 95, 18, 14 ];
Norm9Classes[[ 14, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1 ]]:=[ 35, 189, 14 ];
Norm9Classes[[ 14, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0 ]]:=[ 56, 48, 14 ];
Norm9Classes[[ 14, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1 ]]:=[ 80, 108, 14 ];
Norm9Classes[[ 14, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1 ]]:=[ 35, 189, 14 ];
Norm9Classes[[ 14, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0 ]]:=[ 56, 78, 14 ];
Norm9Classes[[ 7, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0 ]]:=[ 88, 93, 7 ];
Norm9Classes[[ 7, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0 ]]:=[ 73, 138, 7 ];
Norm9Classes[[ 24, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1 ]]:=[ 73, 158, 24 ];
Norm9Classes[[ 24, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0 ]]:=[ 41, 28, 24 ];
Norm9Classes[[ 24, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1 ]]:=[ 88, 68, 24 ];
Norm9Classes[[ 24, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0 ]]:=[ 71, 58, 24 ];
Norm9Classes[[ 8, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0 ]]:=[ 90, 200, 8 ];
Norm9Classes[[ 8, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0 ]]:=[ 15, 20, 8 ];
Norm9Classes[[ 14, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1 ]]:=[ 15, 17, 14 ];
Norm9Classes[[ 14, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1 ]]:=[ 22, 43, 14 ];
Norm9Classes[[ 14, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1 ]]:=[ 38, 183, 14 ];
Norm9Classes[[ 14, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1 ]]:=[ 90, 137, 14 ];
Norm9Classes[[ 14, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1 ]]:=[ 97, 13, 14 ];
Norm9Classes[[ 14, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1 ]]:=[ 53, 153, 14 ];
Norm9Classes[[ 7, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1 ]]:=[ 21, 127, 7 ];
Norm9Classes[[ 7, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1 ]]:=[ 21, 202, 7 ];

Gold:=AbelianGroup(Norm9Invariants);assert [Eltseq(r):r in Rows(RelationMatrix(Gold))] eq Norm9Relations;
O9:=MaximalOrderFinite(CY9);OA:=MaximalOrderFinite(AF);OK:=MaximalOrderFinite(K);
nJ,dJ:=IntegralSplit(jK,OK);nA,dA:=IntegralSplit(aAS,OA);
funcs9:=[CY9!AF!K!BaseField(K).1,CY9!AF!K.1,(CY9!AF!nJ)*CY9!t3,(CY9!dA)*zz];assert &and[f in O9:f in funcs9];
function EvalRatNative57(a,xv)
 L:=Parent(xv);return Evaluate(PolynomialRing(L)!Numerator(a),xv)/Evaluate(PolynomialRing(L)!Denominator(a),xv);
end function;
function EvalKNative57(a,xv,yv)
 cc:=Eltseq(K!a);return &+[EvalRatNative57(cc[i],xv)*yv^(i-1):i in [1..#cc]];
end function;
function EvalAFNative57(a,xv,yv,tv)
 cc:=Eltseq(AF!a);return &+[EvalKNative57(cc[i],xv,yv)*tv^(i-1):i in [1..#cc]];
end function;
function PlaceFromKey57(key)
 d:=key[1];PP<X>:=PolynomialRing(GF(2));mp:=PP!key[7..7+d];L<a>:=ext<GF(2)|mp>;
 coords:=[&+[L!key[8+d+i*d+j]*a^j:j in [0..d-1]]:i in [0..4]];
 assert &+[L!key[1+i]*coords[i]:i in [1..5]] eq a;
 jv,xv,yv,tv,zv:=Explode(coords);assert EvalKNative57(jK,xv,yv) eq jv;
 assert tv^4+tv^3+1/jv eq 0 and zv^2+zv+EvalAFNative57(aAS,xv,yv,tv) eq 0;
 vals:=[xv,yv,EvalKNative57(K!nJ,xv,yv)*tv,EvalAFNative57(AF!dA,xv,yv,tv)*zv];found:=false;
 for i in [1..4] do
  pp:=MinimalPolynomial(vals[i],GF(2));if Degree(pp) eq d then gen:=funcs9[i];val:=vals[i];found:=true;break;end if;
 end for;
 assert found;MB:=Matrix(GF(2),[Eltseq(val^i,GF(2)):i in [0..d-1]]);assert Rank(MB) eq d;
 gens:=[O9!Evaluate(pp,gen)];
 for i in [1..4] do
  co:=Solution(MB,Vector(GF(2),Eltseq(vals[i],GF(2))));
  Append(~gens,O9!(funcs9[i]-&+[CY9!co[k+1]*gen^k:k in [0..d-1]]));
 end for;
 II:=ideal<O9|gens>;assert IsPrime(II);P:=Place(II);assert Degree(P) eq d;return P;
end function;
if not assigned NormFreshFile then NormFreshFile:="C57/data/C57_all_genus9_fresh_data.m";end if;
Write(NormFreshFile,Sprintf("FreshNorm9Invariants:=%o;FreshNorm9Relations:=%o;FreshNorm9Classes:=AssociativeArray();",Invariants(Gfresh),[Eltseq(r):r in Rows(RelationMatrix(Gfresh))]):Overwrite:=true);
FreshClasses:=AssociativeArray();allkeys:=Keys(Norm9Classes);
for key in allkeys do
 P:=PlaceFromKey57(key);cl:=phiFresh(Divisor(P));FreshClasses[key]:=cl;
 Write(NormFreshFile,Sprintf("FreshNorm9Classes[%o]:=%o;",key,Eltseq(cl)));
 printf "C57_NORM9_FRESH_POINT_CLASS degree=%o class=%o time=%o\n",Degree(P),Eltseq(cl),Cputime(tn57);
end for;
inputkeys:=[k:k in allkeys|k[1] eq 8];targetkeys:=[k:k in allkeys|k[1] eq 24];assert #inputkeys eq 2 and #targetkeys eq 4;
assert &and[k[2..6] eq [1,0,0,0,0]:k in inputkeys cat targetkeys];
firstkey:=targetkeys[1];targetA:=[k:k in targetkeys|k[1..#k-24] eq firstkey[1..#firstkey-24]];targetB:=[k:k in targetkeys|k notin targetA];assert #targetA eq 2 and #targetB eq 2;
checked57:=0;zeros57:=0;
for i,a,b in [1..2] do
 po:=Gold!Norm9Classes[inputkeys[i]];spo:=Gold!Norm9Classes[inputkeys[3-i]];
 vo:=Gold!Norm9Classes[targetA[3-a]]+Gold!Norm9Classes[targetB[3-b]]-Gold!Norm9Classes[targetA[a]]-Gold!Norm9Classes[targetB[b]]+4*(spo-po);
 vf:=FreshClasses[targetA[3-a]]+FreshClasses[targetB[3-b]]-FreshClasses[targetA[a]]-FreshClasses[targetB[b]]+4*(FreshClasses[inputkeys[3-i]]-FreshClasses[inputkeys[i]]);
 assert (vo eq Gold!0) eq (vf eq Gfresh!0);checked57+:=1;if vf eq Gfresh!0 then zeros57+:=1;end if;
 printf "C57_NORM9_FRESH_ZERO_STATUS input=%o A=%o B=%o principal=%o\n",i,a,b,vf eq Gfresh!0;
end for;
assert checked57 eq 8 and zeros57 eq 2;

OldGood78Arithmetic:=[*
[ 105, 210, 0 ],
[
[ 105, 0, 0 ],
[ 0, 210, 0 ]
],
[* [* 1, 7,
<7, [ 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0 ]>,
[ 21, 202, 7 ],
[ 21, 127, 7 ],
[* [* 14, 1, 1,
[ 15, 17, 14 ],
[ 90, 137, 14 ]
*], [* 14, 1, 1,
[ 22, 43, 14 ],
[ 97, 13, 14 ]
*], [* 14, 1, 1,
[ 38, 183, 14 ],
[ 53, 153, 14 ]
*] *],
[
[ 45, 60, 0 ],
[ 90, 90, 0 ],
[ 90, 0, 0 ],
[ 30, 30, 0 ],
[ 75, 0, 0 ],
[ 15, 30, 0 ],
[ 15, 150, 0 ],
[ 60, 180, 0 ]
]
*], [* 2, 7,
<7, [ 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0 ]>,
[ 73, 138, 7 ],
[ 88, 93, 7 ],
[* [* 14, 1, 1,
[ 95, 18, 14 ],
[ 80, 108, 14 ]
*], [* 14, 1, 1,
[ 35, 189, 14 ],
[ 35, 189, 14 ]
*], [* 14, 1, 1,
[ 56, 48, 14 ],
[ 56, 78, 14 ]
*] *],
[
[ 75, 120, 0 ],
[ 45, 90, 0 ],
[ 75, 120, 0 ],
[ 45, 90, 0 ],
[ 75, 180, 0 ],
[ 45, 150, 0 ],
[ 75, 180, 0 ],
[ 45, 150, 0 ]
]
*], [* 3, 8,
<8, [ 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1 ]>,
[ 15, 20, 8 ],
[ 90, 200, 8 ],
[* [* 24, 1, 1,
[ 73, 158, 24 ],
[ 88, 68, 24 ]
*], [* 24, 1, 1,
[ 41, 28, 24 ],
[ 71, 58, 24 ]
*] *],
[
[ 45, 150, 0 ],
[ 75, 180, 0 ],
[ 0, 0, 0 ],
[ 30, 30, 0 ]
]
*] *]
*];
FreshByOldVector:=AssociativeArray();
for key in allkeys do
 vec:=Norm9Classes[key];
 if IsDefined(FreshByOldVector,vec) then assert FreshByOldVector[vec] eq FreshClasses[key];else FreshByOldVector[vec]:=FreshClasses[key];end if;
end for;
allgoodtests57:=0;allgoodzeros57:=0;good7zeros57:=0;
for rec57 in OldGood78Arithmetic[3] do
 pc:=Gold!rec57[4];sc:=Gold!rec57[5];pf:=FreshByOldVector[rec57[4]];sf:=FreshByOldVector[rec57[5]];
 choices57:=[<Gold!0,Gfresh!0>];
 for target57 in rec57[6] do
  assert target57[3] eq 1;mm:=target57[2];
  olddiff:=Gold!target57[5]-Gold!target57[4];newdiff:=FreshByOldVector[target57[5]]-FreshByOldVector[target57[4]];
  opts57:=[<(2*a-mm)*olddiff,(2*a-mm)*newdiff>:a in [0..mm]];
  choices57:=[<c[1]+v[1],c[2]+v[2]>:c in choices57,v in opts57];
 end for;
 choices57:=[<c[1]+4*(sc-pc),c[2]+4*(sf-pf)>:c in choices57];
 assert Sort([Sprint(Eltseq(c[1])):c in choices57]) eq Sort([Sprint(v):v in rec57[7]]);
 zeroscase57:=0;
 for cc57 in choices57 do
  assert (cc57[1] eq Gold!0) eq (cc57[2] eq Gfresh!0);allgoodtests57+:=1;
  if cc57[2] eq Gfresh!0 then allgoodzeros57+:=1;zeroscase57+:=1;if rec57[2] eq 7 then good7zeros57+:=1;end if;end if;
 end for;
 printf "C57_GENUS9_FRESH_GOOD_CASE degree=%o choices=%o zeros=%o\n",rec57[2],#choices57,zeroscase57;
end for;
assert #allkeys eq 22 and allgoodtests57 eq 20 and allgoodzeros57 eq 1 and good7zeros57 eq 0;
printf "C57_GENUS9_FRESH_GOOD_DONE points=22 good7_tests=16 good7_zeros=0 good8_norm_tests=4 good8_norm_zeros=1 complete=true time=%o\n",Cputime(tn57);
printf "C57_NORM9_FRESH_DONE points=6 zero_tests=%o zeros=%o independent_class_basis=true complete=true time=%o\n",checked57,zeros57,Cputime(tn57);
catch err
 print "C57_NORM9_FRESH_FAILED";print err`Object;
end try;
