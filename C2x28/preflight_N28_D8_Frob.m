SetColumns(0);SetQuitOnError(true);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
J:=JH(56,[29]);assert Dimension(J) eq 25;f:=FrobeniusPolynomial(J,3);assert Degree(f) eq 50;
tr:=[];cs:=[];
for n in [1..8] do
 s:=-n*Coefficient(f,50-n);for j in [1..n-1] do s-:=Coefficient(f,50-j)*tr[n-j];end for;
 Append(~tr,Integers()!s);a:=3^n+1-s;
 for d in Divisors(n) do if d lt n then a-:=d*cs[d];end if;end for;
 assert a mod n eq 0;Append(~cs,Integers()!(a div n));
end for;
assert cs eq [12,6,4,0,48,66,312,606];
Write("C2x28/data/N28_D8_Frob_arithmetic.txt",Sprint([*25,Coefficients(f),cs*]):Overwrite:=true);
print "N28_D8_INDEPENDENT_FROB_DONE",cs,Cputime(tt);quit;
