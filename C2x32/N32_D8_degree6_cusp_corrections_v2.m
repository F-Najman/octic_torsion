/* Fresh sole degree6-support case, with all17 exact overapproximating
   quadratic cusp corrections for q5 and a generic T5 reference check. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
X:=MDX11(2,32,GF(3));C:=Curve(X);FF:=FunctionField(C);cs:=Cusps(X);cub:=[c:c in cs|Degree(c) eq 2];assert #cub eq 8;
E0:=EllipticCurve([FF!f:f in X`_E]);Q0:=E0![FF!f:f in X`_Q];ai:=aInvariants(E0);b:=Q0[1];cc:=ai[2];dd:=ai[4];
assert ai[1] eq 0 and ai[3] eq 0 and ai[5] eq 0 and dd eq b*(1-b-cc);
fs:=[b,dd,cc^2-4*dd];chosen:=[];
for k in [0,1,2,4,8,16,3,5,7,9,11,13,15] do
 if k ne 0 then Append(~fs,(k*Q0)[1]);Append(~chosen,k);end if;
 sigs:=[[Valuation(f,c):f in fs]:c in cub];us:=Sort(Setseq(Seqset(sigs)));sizes:=[#[t:t in sigs|t eq ss]:ss in us];
 print "N32_D8_QUADRATIC_UNIT_BLOCKS",chosen,sizes,Cputime(t0);
 if Maximum(sizes) le 2 then break;end if;
end for;
assert Maximum(sizes) le 2;
pairs:=[[i:i in [1..8]|sigs[i] eq ss]:ss in us|#[t:t in sigs|t eq ss] eq 2];assert #pairs gt 0;
vs:=[];
for p in pairs do
 for sign in [-1,1] do v:=[0:i in [1..8]];v[p[1]]:=-sign;v[p[2]]:=sign;Append(~vs,v);end for;
end for;assert #vs eq 2*#pairs;
print "N32_D8_DEGREE6_CUSP_PAIRS",pairs,Cputime(t0);
input:=eval Read("C2x32/data/N32_D8_fullaut_Hecke_input.txt");records:=input[1];assert #records eq 50;
R:=PolynomialRing(GF(3));tests:=[* *];remaining:=[];
for n in [2] do
 degree,idx,def,vals:=Explode(records[n]);assert degree eq 6 and idx eq n;
 f:=R!def;assert IsIrreducible(f) and Degree(f) eq 6;K<z>:=ext<GF(3)|f>;
 bc:=[Evaluate(R!v,z):v in vals];x:=MDPlaceCommonZeros(X`_coordinates,bc);assert Degree(x) eq 6;killed:=false;
 for q in [5] do
  T:=HeckeOperatorFast(X,q,x);A:=T-Divisor(x)-q*Divisor(DiamondOperator(X,q,x));assert Degree(A) eq 0;
  if n eq 2 and q eq 5 then
   E:=EllipticCurve(X,x);Ls:=LevelStructure(X,x);ys:=[];ms:=[];
   for phi in MDIsogenies(E,5) do
    La:=rec<recformat<P:PtEll,Q:PtEll>|P:=phi(Ls`P),Q:=phi(Ls`Q)>;
    y:=ModuliPointFast(X,Codomain(phi),La);ii:=Index(ys,y);
    if ii eq 0 then Append(~ys,y);Append(~ms,1);else ms[ii]+:=1;end if;
   end for;
   ref:=DivisorGroup(C)!0;
   for i in [1..#ys] do mu:=ms[i]*6/Degree(ys[i]);assert Denominator(mu) eq 1;ref+:=Integers()!mu*Divisor(ys[i]);end for;
   assert ref eq T;print "N32_D8_FRESH_T5_GENERIC_REFERENCE_EQUAL",Degree(ref),Cputime(t0);
  end if;
  corrections:=q eq 19 select [[0,0,0,0,0,0,0,0]] else [[0,0,0,0,0,0,0,0]] cat [ [v[i]+q*w[i]:i in [1..8]]:v,w in vs];
  assert #Seqset(corrections) eq (q eq 19 select 1 else 1+#vs^2);flags:=[];
  for coeff in corrections do
   D:=A+&+[coeff[i]*Divisor(cub[i]):i in [1..8]];assert Degree(D) eq 0;
   pr:=IsPrincipal(D);Append(~flags,pr select 1 else 0);
  end for;
  Append(~tests,[*n,q,corrections,flags*]);print "N32_D8_DEGREE6_CORRECTION_TEST",n,q,&+flags,Cputime(t0);
  if &+flags eq 0 then killed:=true;break;end if;
 end for;
 if not killed then Append(~remaining,n);end if;
end for;
Write("C2x32/data/N32_D8_degree6_correction_v2_arithmetic.txt",Sprint([*sigs,pairs,vs,tests,remaining*]):Overwrite:=true);
print "N32_D8_DEGREE6_CORRECTIONS_DONE",remaining,Cputime(t0);quit;
