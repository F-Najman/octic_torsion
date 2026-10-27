/* Fresh degree5-support cases, with all17 exact overapproximating cubic
   cusp corrections for q5/11; q19 literal-annihilation fallback. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
X:=MDX11(2,28,GF(3));C:=Curve(X);FF:=FunctionField(C);cs:=Cusps(X);cub:=[c:c in cs|Degree(c) eq 3];assert #cub eq 4;
E0:=EllipticCurve([FF!f:f in X`_E]);Q0:=E0![FF!f:f in X`_Q];ai:=aInvariants(E0);b:=Q0[1];cc:=ai[2];dd:=ai[4];
assert ai[1] eq 0 and ai[3] eq 0 and ai[5] eq 0 and dd eq b*(1-b-cc);
fs:=[b,dd,cc^2-4*dd];for k in [1..14] do Append(~fs,(k*Q0)[1]);end for;
sigs:=[[Valuation(f,c):f in fs]:c in cub];us:=Sort(Setseq(Seqset(sigs)));assert #us eq 2;
pairs:=[[i:i in [1..4]|sigs[i] eq s]:s in us];assert &and[#p eq 2:p in pairs];
vs:=[];
for p in pairs do
 for sign in [-1,1] do v:=[0:i in [1..4]];v[p[1]]:=-sign;v[p[2]]:=sign;Append(~vs,v);end for;
end for;assert #vs eq 4;
print "N28_D8_QUINTIC_CUSP_PAIRS",pairs,Cputime(t0);
input:=eval Read("C2x28/data/N28_D8_fullaut_Hecke_input.txt");records:=input[1];assert #records eq 49;
R:=PolynomialRing(GF(3));tests:=[* *];remaining:=[];
for n in [1,2] do
 degree,idx,def,vals:=Explode(records[n]);assert degree eq 5 and idx eq n;
 f:=R!def;assert IsIrreducible(f) and Degree(f) eq 5;K<z>:=ext<GF(3)|f>;
 bc:=[Evaluate(R!v,z):v in vals];x:=MDPlaceCommonZeros(X`_coordinates,bc);assert Degree(x) eq 5;killed:=false;
 for q in [5,11,19] do
  T:=HeckeOperatorFast(X,q,x);A:=T-Divisor(x)-q*Divisor(DiamondOperator(X,q,x));assert Degree(A) eq 0;
  if n eq 1 and q eq 5 then
   E:=EllipticCurve(X,x);Ls:=LevelStructure(X,x);ys:=[];ms:=[];
   for phi in MDIsogenies(E,5) do
    La:=rec<recformat<P:PtEll,Q:PtEll>|P:=phi(Ls`P),Q:=phi(Ls`Q)>;
    y:=ModuliPointFast(X,Codomain(phi),La);ii:=Index(ys,y);
    if ii eq 0 then Append(~ys,y);Append(~ms,1);else ms[ii]+:=1;end if;
   end for;
   ref:=DivisorGroup(C)!0;
   for i in [1..#ys] do mu:=ms[i]*5/Degree(ys[i]);assert Denominator(mu) eq 1;ref+:=Integers()!mu*Divisor(ys[i]);end for;
   assert ref eq T;print "N28_D8_FRESH_T5_GENERIC_REFERENCE_EQUAL",Degree(ref),Cputime(t0);
  end if;
  corrections:=q eq 19 select [[0,0,0,0]] else [[0,0,0,0]] cat [ [v[i]+q*w[i]:i in [1..4]]:v,w in vs];
  assert #Seqset(corrections) eq (q eq 19 select 1 else 17);flags:=[];
  for coeff in corrections do
   D:=A+&+[coeff[i]*Divisor(cub[i]):i in [1..4]];assert Degree(D) eq 0;
   pr:=IsPrincipal(D);Append(~flags,pr select 1 else 0);
  end for;
  Append(~tests,[*n,q,corrections,flags*]);print "N28_D8_QUINTIC_CORRECTION_TEST",n,q,&+flags,Cputime(t0);
  if &+flags eq 0 then killed:=true;break;end if;
 end for;
 if not killed then Append(~remaining,n);end if;
end for;
Write("C2x28/data/N28_D8_quintic_correction_arithmetic.txt",Sprint([*sigs,pairs,vs,tests,remaining*]):Overwrite:=true);
print "N28_D8_QUINTIC_CORRECTIONS_DONE",remaining,Cputime(t0);quit;
