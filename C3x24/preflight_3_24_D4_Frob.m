SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
J:=JH(72,[25]);g:=Dimension(J);assert g eq 25;rows:=[* *];
for p in [7,13,19] do
 f:=FrobeniusPolynomial(J,p);assert Degree(f) eq 2*g;tr:=[];cl:=[];
 for n in [1..4] do
  a:=-n*Coefficient(f,2*g-n);
  for j in [1..n-1] do a-:=Coefficient(f,2*g-j)*tr[n-j];end for;
  Append(~tr,Integers()!a);num:=p^n+1-Integers()!a;
  for e in Divisors(n) do if e lt n then num-:=e*cl[e];end if;end for;
  assert num mod n eq 0;Append(~cl,num div n);
 end for;
 Append(~rows,[*p,Coefficients(f),cl,Integers()!Evaluate(f,1)*]);print "M3N24_D4_FROB_ROW",p,cl,Evaluate(f,1),Cputime(t0);
end for;
Write("C3x24/data/preflight_3_24_D4_Frob_arithmetic.txt",Sprint([*g,rows*]):Overwrite:=true);
print "M3N24_D4_FROB_DONE complete=true",Cputime(t0);quit;
