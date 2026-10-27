AttachSpec("mdmagma/v2/mdmagma.spec");
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(2*10^9);tm:=Cputime();
X19:=MDX1(19,GF(2));C19:=Curve(X19);FF19:=FunctionField(C19);
E19:=EllipticCurve(X19`_E);K19,m19:=AlgorithmicFunctionField(FF19);OK19:=MaximalOrderFinite(K19);nY19,dY19:=IntegralSplit(K19.1,OK19);
AA:=[m19(a):a in aInvariants(E19)];a1n,a2n,a3n,a4n,a6n:=Explode(AA);r19n:=a3n/a1n;t19n:=(a4n+r19n^2)/a1n;d19n:=(a2n+r19n)/a1n^2;
function EvalRat(f,t)
 return Evaluate(Numerator(f),t)/Evaluate(Denominator(f),t);
end function;
function EvalNative(f,t,y)
 return &+[EvalRat(Eltseq(f)[i],t)*y^(i-1):i in [1..#Eltseq(f)]];
end function;
function DecodeKey57(key)
 d:=key[1];bb:=key[2..3];mp:=Polynomial(GF(2),key[4..d+4]);assert Degree(mp) eq d and IsIrreducible(mp);
 L<a>:=ext<GF(2)|mp>;rows:=[key[d+5+i*d..d+4+(i+1)*d]:i in [0..2]];v:=[Evaluate(Polynomial(GF(2),rr),a):rr in rows];
 t:=v[1];nv:=v[2];w:=v[3];assert bb[1]*t+bb[2]*nv eq a;
 dy:=EvalRat(BaseField(K19)!dY19,t);assert dy ne 0;y:=nv/dy;
 assert &+[EvalRat(Coefficient(DefiningPolynomial(K19),i),t)*y^i:i in [0..Degree(DefiningPolynomial(K19))]] eq 0;
 aa:=[EvalNative(f,t,y):f in AA];E:=EllipticCurve(aa);R19:=E![L!0,L!0];assert Order(R19) eq 19;
 r:=EvalNative(r19n,t,y);tt:=EvalNative(t19n,t,y);dd:=EvalNative(d19n,t,y);xx:=w^2+w+dd;
 assert xx^4+xx^3+1/jInvariant(E) eq 0;yy:=xx*(w+xx);
 R3:=E![aa[1]^2*xx+r,aa[1]^3*yy+tt];assert Order(R3) eq 3 and Order(R3+13*R19) eq 57;
 return E,R3,R19;
end function;
function ModuliKey57(E,R3,R19)
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
 return key;
end function;
load "C57/data/full85_relative_divisors_resume_pointkeys.m";
assert #InputKeys eq 3 and #OutputKeys eq 2 and #Diamond5OutputKeys eq 2 and HeckeObservations eq [3,3];
E,R3,R19:=DecodeKey57(InputKeys[1]);assert #E eq 228 and Degree(MinimalPolynomial(jInvariant(E),GF(2))) eq 8;
allkeys:=[ModuliKey57(E,5^k*R3,5^k*R19):k in [0..17]];
assert #Seqset(allkeys) eq 18 and &and[k[1] eq 8:k in allkeys] and allkeys[1..3] eq InputKeys;
for key in OutputKeys cat Diamond5OutputKeys do EC,R3C,R19C:=DecodeKey57(key);assert ModuliKey57(EC,R3C,R19C) eq key;end for;
phs:=MDIsogenies(E,5);assert #phs eq 6;fresh:=[];freshd5:=[];
for ph in phs do
 EC:=Codomain(ph);rr3:=ph(R3);rr19:=ph(R19);Append(~fresh,ModuliKey57(EC,rr3,rr19));Append(~freshd5,ModuliKey57(EC,5*rr3,5*rr19));
end for;
assert Seqset(fresh) eq Seqset(OutputKeys) and &and[#[k:k in fresh|k eq key] eq 3:key in OutputKeys];
for i in [1..#fresh] do assert freshd5[i] eq Diamond5OutputKeys[Index(OutputKeys,fresh[i])];end for;
printf "C57_GOOD8_COMPLETION_MODULI_REPLAY_PASS all18inputorbit=true actual6fiveisogenies=true two24outputs=true signed3pointretained=true time=%o\n",Cputime(tm);
quit;
