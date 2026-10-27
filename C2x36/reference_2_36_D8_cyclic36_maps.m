SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
Attach("common/fast_hecke_cyclic_commonzero.m");
Y:=MDX1(36,GF(5));R:=PolynomialRing(GF(5));
a:=eval Read("C2x36/data/M2N36_D8_direct_full3_arithmetic.txt");
K<z>:=ext<GF(5)|R!a[2]>;key:=a[9][1][2];b:=Evaluate(R!key[1..3],z);c:=Evaluate(R!key[4..6],z);
E:=EllipticCurve([K|0,c,0,b*(1-b-c),0]);P:=E![0,0];Q:=E![b,b];
assert 36*Q eq E!0 and 18*Q ne P and 18*Q ne E!0 and 12*Q ne E!0;
for mapno in [1..3] do
 if mapno eq 1 then EG:=E;QG:=Q;
 else T:=mapno eq 2 select P else P+18*Q;EG,phi:=IsogenyFromKernel(E,PolynomialRing(K).1-T[1]);QG:=phi(Q);end if;
 assert 36*QG eq EG!0 and 18*QG ne EG!0 and 12*QG ne EG!0;
 L:=rec<recformat<P:PtEll>|P:=QG>;
 fast:=ModuliPointCyclicCommonZeros(Y,EG,L);ref:=ModuliPoint(Y,EG,L);
 assert fast eq ref and Degree(fast) eq 3;
 print "M2N36_CYCLIC36_NATIVE_REFERENCE",mapno,Degree(fast),Cputime(t0);
end for;
print "M2N36_CYCLIC36_NATIVE_REFERENCE_DONE";quit;
