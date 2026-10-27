SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
J:=JH(64,[33]);g:=Dimension(J);print "OCTIC_M2N32_ACTUAL_GENUS",g;
f:=FrobeniusPolynomial(J,3);assert Degree(f) eq 2*g;tr:=[];cl:=[];
for n in [1..8] do
 a:=-n*Coefficient(f,2*g-n);
 for j in [1..n-1] do a-:=Coefficient(f,2*g-j)*tr[n-j];end for;
 Append(~tr,Integers()!a);num:=3^n+1-Integers()!a;
 for e in Divisors(n) do if e lt n then num-:=e*cl[e];end if;end for;
 assert num mod n eq 0;Append(~cl,num div n);
end for;
print "OCTIC_M2N32_FROB_COUNTS",g,cl,Cputime(t0);
Write("C2x32/data/preflight_2_32_Frob_arithmetic.txt",Sprint([*g,Coefficients(f),cl*]):Overwrite:=true);
print "OCTIC_M2N32_FROB_DONE complete=true",Cputime(t0);quit;
