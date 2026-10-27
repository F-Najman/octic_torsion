SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
J:=JH(63,[22]);g:=Dimension(J);print "M3N21_ACTUAL_GENUS",g;
f:=FrobeniusPolynomial(J,2);assert Degree(f) eq 2*g;tr:=[];
for n in [1..8] do
 a:=-n*Coefficient(f,2*g-n);
 for j in [1..n-1] do a-:=Coefficient(f,2*g-j)*tr[n-j];end for;
 Append(~tr,Integers()!a);
end for;
cl:=[];
for n in [1..4] do
 num:=4^n+1-tr[2*n];
 for e in Divisors(n) do if e lt n then num-:=e*cl[e];end if;end for;
 assert num mod n eq 0;Append(~cl,num div n);
end for;
Write("C3x21/data/preflight_3_21_Frob_arithmetic.txt",Sprint([*g,Coefficients(f),cl*]):Overwrite:=true);
print "M3N21_F4_CLOSED_COUNTS",g,cl,Cputime(t0);
print "M3N21_F4_FROB_DONE complete=true",Cputime(t0);quit;
