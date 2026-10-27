/* Bounded characteristic-zero RR for the first short quotient cusp lift.
   All earlier exact finite certificates remain unchanged. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
a:=eval Read("C2x28/data/N28_D8_quotient_cusp_RR_arithmetic.txt");
sigs,ds,base,lifts:=Explode(a[1..4]);
Y:=MDX1(28,Rationals());C:=Curve(Y);assert Genus(Y) eq 10;
print "N28_D8_QUOTIENT_CHARZERO_MODEL",Cputime(t0);
raw:=CuspOrbitsQ(Y);os:=[Type(o) eq PlcCrvElt select Divisor(o) else o:o in raw];
ss:=[CuspSignature(Y,Support(o)[1]):o in os];order:=Sort([<ss[i],i>:i in [1..#os]]);
os:=[os[x[2]]:x in order];ss:=[x[1]:x in order];assert ss eq sigs and [Degree(o):o in os] eq ds;
c0:=Support(os[base])[1];print "N28_D8_QUOTIENT_CHARZERO_CUSPS",#os,Cputime(t0);
L:=8*Divisor(c0)+&+[lifts[1][i]*(os[i]-ds[i]*Divisor(c0)):i in [1..#os]];
assert Degree(L) eq 8;DD:=FunctionFieldDivisor(L);
print "N28_D8_QUOTIENT_CHARZERO_DIRECT_BASIS_START",Cputime(t0);
bb:=Basis(DD:Reduction:=false,Simplification:="None");
print "N28_D8_QUOTIENT_CHARZERO_RR_DIMENSION",#bb,Cputime(t0);
if #bb eq 1 then
 T:=DD+Divisor(bb[1]);assert IsEffective(T) and Degree(T) eq 8;sp,ms:=Support(T);
 print "N28_D8_QUOTIENT_CHARZERO_EFFECTIVE",[Degree(p):p in sp],ms,Cputime(t0);
 print "N28_D8_QUOTIENT_CHARZERO_FUNCTION",bb[1];
 print "N28_D8_QUOTIENT_CHARZERO_PLACES",sp;
end if;
print "N28_D8_QUOTIENT_CHARZERO_RR_DONE",Cputime(t0);quit;
