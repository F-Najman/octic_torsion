SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
function CuspUnitsM4(X)
 FF:=FunctionField(Curve(X));E:=EllipticCurve([FF!f:f in X`_E]);P:=E![FF!f:f in X`_P];Q:=E![FF!f:f in X`_Q];
 ai:=aInvariants(E);B:=ai[2];assert ai eq [FF!1,B,B,FF!0,FF!0] and P eq E![0,0];
 assert Discriminant(E) eq B^4*(1-16*B);fs:=[B,1-16*B];
 for k in [1..8] do for j in [0..3] do
  A:=k*Q+j*P;assert A ne E!0 and A[1] ne 0 and A[2] ne 0 and A[1]+B ne 0;
  fs cat:=[A[1],A[2],A[1]+B,A[2]+A[1]+B];
 end for;end for;
 for k in [1..8] do A:=k*Q;for j in [1..3] do R:=A+j*P;assert A ne R and A ne -R;Append(~fs,A[1]-R[1]);end for;end for;
 rawcs:=Cusps(X);groups:=AssociativeArray();
 for c in rawcs do sig:=[Valuation(f,c):f in fs];
  if IsDefined(groups,sig) then Append(~groups[sig],c);else groups[sig]:=[c];end if;
 end for;
 sigs:=Sort(Setseq(Keys(groups)));assert #sigs eq 24;
 return [&+[Divisor(c):c in groups[sig]]:sig in sigs],sigs;
end function;
a:=eval Read("C4x16/data/M4N16_D4_doubled_cusp_lift.txt");tests:=[* *];
for p in [13] do
 X:=MDX11(4,16,GF(p));blocks,sigs:=CuspUnitsM4(X);ds:=[Degree(b):b in blocks];assert sigs eq a[4] and ds eq a[5];
 c0:=Support(blocks[a[6]])[1];v:=a[7];D:=8*Divisor(c0)+&+[v[i]*(blocks[i]-ds[i]*Divisor(c0)):i in [1..24]];
 V,mp:=RiemannRochSpace(D);dim:=Dimension(V);degs:=[];ms:=[];
 if dim eq 1 then T:=D+Divisor(mp(V.1));assert IsEffective(T) and Degree(T) eq 8;sp,ms:=Support(T);degs:=[Degree(c):c in sp];print "M4N16_D4_FINAL_RR_FUNCTION_WITNESS",mp(V.1);print "M4N16_D4_FINAL_EFFECTIVE_DIVISOR",T;assert degs eq [4,4] and ms eq [1,1];end if;
 Append(~tests,[*p,dim,degs,ms*]);print "M4N16_D4_DOUBLED_CUSP_FINITE_RR",tests[#tests],Cputime(t0);
end for;
Write("C4x16/data/M4N16_D4_final_RR13_arithmetic.txt",Sprint(tests):Overwrite:=true);
print "M4N16_D4_DOUBLED_CUSP_FINITE_RR_DONE",Cputime(t0);quit;
