AttachSpec("mdmagma/v2/mdmagma.spec");
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(32*10^9);tt:=Cputime();
X57:=MDX1(57,GF(2));C57:=Curve(X57);FF57:=FunctionField(C57);K57,m57:=AlgorithmicFunctionField(FF57);
assert Degree(K57) eq 37;print "C57_NATIVE37_GENUS_START",Cputime(tt);assert Genus(K57) eq 85;
O57:=MaximalOrderFinite(K57);nY57,dY57:=IntegralSplit(K57.1,O57);co57:=[K57!BaseField(K57).1,nY57];
assert &and[f in O57:f in co57];
print "C57_NATIVE37_MODEL_VERIFIED genus85=true degree37=true",Cputime(tt);
// The native plane model has x=t,y=generator. Check this once; then
// evaluate every AF function coefficientwise, avoiding generic inverse maps.
print "C57_NATIVE37_RAW_MAP_IDENTIFICATION_START",Cputime(tt);assert m57(FF57.1) eq K57!BaseField(K57).1 and m57(FF57.2) eq K57.1;print "C57_NATIVE37_RAW_MAP_IDENTIFICATION_PASS",Cputime(tt);
function EvalNative57(f,xx,yy)
 cs:=Eltseq(K57!f);L:=Parent(xx);v:=L!0;
 for i in [1..#cs] do
  num:=Evaluate(Numerator(cs[i]),xx);den:=Evaluate(Denominator(cs[i]),xx);assert den ne 0;
  v+:=(num/den)*yy^(i-1);
 end for;return v;
end function;
function NativePoint57(E,PN,cache)
 assert Order(PN) eq 57;L:=BaseRing(E);xy:=TateNormalForm_xy(E,PN);
 xx,yy:=Explode(xy);
 mp57:=DefiningPolynomial(K57);cs57:=Coefficients(mp57);
 assert &+[Evaluate(Numerator(cs57[i]),xx)/Evaluate(Denominator(cs57[i]),xx)*yy^(i-1):i in [1..#cs57]] eq 0;
 r:=(xx^2*yy-xx*yy+yy-1)/(xx*(xx*yy-1));ss:=(xx*yy-yy+1)/(xx*yy);
 assert TateNormalForm_bc(E,PN) eq [r*ss*(r-1),ss*(r-1)];
 assert EvalNative57(dY57,xx,yy) ne 0;
 vals:=[EvalNative57(f,xx,yy):f in co57];print "C57_NATIVE37_RAW_AF_VALUES_DONE",Cputime(tt);d:=LCM([Degree(MinimalPolynomial(v,GF(2))):v in vals]);found:=false;
 for bits in [[1,0],[0,1],[1,1]] do val:=&+[bits[i]*vals[i]:i in [1..2]];mp:=MinimalPolynomial(val,GF(2));if Degree(mp) eq d then bb:=bits;found:=true;break;end if;end for;assert found;
 MB:=Matrix(GF(2),[Eltseq(val^i,GF(2)):i in [0..d-1]]);assert Rank(MB) eq d;
 rows:=[Eltseq(Solution(MB,Vector(GF(2),Eltseq(v,GF(2))))):v in vals];
 key:=[d] cat bb cat [Integers()!a:a in Eltseq(mp)] cat [Integers()!a:a in &cat rows];
 if IsDefined(cache,key) then return cache[key],key;end if;
 gen:=&+[bb[i]*co57[i]:i in [1..2]];gens:=[O57!Evaluate(mp,gen)];
 for i in [1..2] do Append(~gens,O57!(co57[i]-&+[K57!rows[i][k+1]*gen^k:k in [0..d-1]]));end for;
 print "C57_NATIVE37_GENERATORS_READY",d,Cputime(tt);I:=ideal<O57|gens>;print "C57_NATIVE37_IDEAL_READY",d,Cputime(tt);assert Degree(Norm(I)) eq d;
 // The coordinate evaluation surjects onto GF(2^d); equality of the
 // ideal norm degree certifies its unique normalization prime.
 printf "C57_NATIVE37_EXACT_POINT_IDEAL degree=%o time=%o\n",d,Cputime(tt);return I,key;
end function;
F<a>:=GF(256);ces:=[**];
for j in F do for E in Twists(EllipticCurveFromjInvariant(j)) do if #E mod 57 eq 0 then assert #E eq 228;Append(~ces,E);end if;end for;end for;assert #ces eq 8;
E:=ces[1];G,g:=AbelianGroup(E);assert Invariants(G) eq [228];PN:=4*g(G.1);assert Order(PN) eq 57 and Degree(MinimalPolynomial(jInvariant(E),GF(2))) eq 8;
cache:=AssociativeArray();ins:=[];inkeys:=[];
for k in [0..2] do I,key:=NativePoint57(E,5^k*PN,cache);cache[key]:=I;Append(~ins,I);Append(~inkeys,key);end for;
assert &and[ins[i] ne ins[j]:i,j in [1..3]|i lt j] and &and[k[1] eq 8:k in inkeys];
assert 5^18 mod 57 eq 1 and 5^9 mod 57 eq 20;
print "C57_NATIVE37_INPUT_COVERAGE jdegree8=true diamond5order18=true independent_total18=true",Cputime(tt);
phs:=MDIsogenies(E,5);assert #phs eq 6;outs:=[];outmods:=[**];outkeys:=[];obs:=[];
for ph in phs do EC:=Codomain(ph);PC:=ph(PN);I,key:=NativePoint57(EC,PC,cache);cache[key]:=I;i:=Index(outs,I);if i eq 0 then Append(~outs,I);Append(~outmods,<EC,PC>);Append(~outkeys,key);Append(~obs,1);else obs[i]+:=1;end if;end for;
assert #outs eq 2 and obs eq [3,3] and &and[k[1] eq 24:k in outkeys];
douts:=[];doutkeys:=[];for md in outmods do I,key:=NativePoint57(md[1],5*md[2],cache);cache[key]:=I;Append(~douts,I);Append(~doutkeys,key);end for;
Write("C57/data/native57_raw_af_q5_pointkeys_20261026.m",Sprintf("InputKeys57:=%o;OutputKeys57:=%o;Diamond5OutputKeys57:=%o;Observations57:=%o;",inkeys,outkeys,doutkeys,obs):Overwrite:=true);
print "C57_NATIVE37_Q5_COVERAGE actual6kernels=true two24outputs=true weights1=true",Cputime(tt);
TI:=&*outs;dTI:=&*douts;
AI:=dTI/TI*ins[2]^4*ins[1]/ins[3]^5;
print "C57_NATIVE37_OPERATOR_IDEAL_READY positive_bound88",Cputime(tt);
DD:=DivisorGroup(K57)!AI;print "C57_NATIVE37_OPERATOR_DIVISOR_READY",Cputime(tt);
pr:=IsPrincipal(DD);printf "C57_NATIVE37_DIAMOND5_A5_DONE principal=%o complete=true time=%o\n",pr,Cputime(tt);
if pr then
 BI:=TI/ins[1]/ins[2]^5;BI5:=BI^5;print "C57_NATIVE37_SCALAR5_IDEAL_READY",Cputime(tt);
 BD5:=DivisorGroup(K57)!BI5;pr5:=IsPrincipal(BD5);printf "C57_NATIVE37_SCALAR5_A5_DONE principal=%o complete=true time=%o\n",pr5,Cputime(tt);
 if pr5 then
  for k in [3..4] do I,key:=NativePoint57(E,5^k*PN,cache);cache[key]:=I;Append(~ins,I);Append(~inkeys,key);end for;
  phs11:=MDIsogenies(E,11);assert #phs11 eq 12;outs11:=[];mods11:=[**];keys11:=[];obs11:=[];
  for ph in phs11 do EC:=Codomain(ph);PC:=ph(PN);I,key:=NativePoint57(EC,PC,cache);cache[key]:=I;i:=Index(outs11,I);if i eq 0 then Append(~outs11,I);Append(~mods11,<EC,PC>);Append(~keys11,key);Append(~obs11,1);else obs11[i]+:=1;end if;end for;
  weights:=[8*obs11[i] div keys11[i][1]:i in [1..#keys11]];assert &and[8*obs11[i] mod keys11[i][1] eq 0:i in [1..#keys11]];assert &+[weights[i]*keys11[i][1]:i in [1..#keys11]] eq 96;
  TI11:=ideal<O57|1>;dTI11:=TI11;dk11:=[];
  for i in [1..#outs11] do TI11*:=outs11[i]^weights[i];md:=mods11[i];I,key:=NativePoint57(md[1],5*md[2],cache);cache[key]:=I;dTI11*:=I^weights[i];Append(~dk11,key);end for;
  Write("C57/data/native57_raw_af_q11_pointkeys_20261026.m",Sprintf("InputKeys57:=%o;OutputKeys57:=%o;Diamond5OutputKeys57:=%o;Observations57:=%o;Weights57:=%o;",inkeys,keys11,dk11,obs11,weights):Overwrite:=true);
  print "C57_NATIVE37_Q11_COVERAGE",[k[1]:k in keys11],weights,Cputime(tt);
  AI11:=dTI11/TI11/ins[2]*ins[1]/ins[5]^11*ins[4]^11;DD11:=DivisorGroup(K57)!AI11;
  pr11:=IsPrincipal(DD11);printf "C57_NATIVE37_DIAMOND5_A11_DONE principal=%o complete=true time=%o\n",pr11,Cputime(tt);
 end if;
end if;
print "C57_NATIVE37_RUN_COMPLETED",Cputime(tt);quit;
