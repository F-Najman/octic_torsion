SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");rows:=[* *];
for tup in [[72,37],[36,19]] do
 J:=JH(tup[1],[tup[2]]);g:=Dimension(J);f:=FrobeniusPolynomial(J,5);assert Degree(f) eq 2*g;tr:=[];cl:=[];
 for n in [1..8] do
  a:=-n*Coefficient(f,2*g-n);for j in [1..n-1] do a-:=Coefficient(f,2*g-j)*tr[n-j];end for;
  Append(~tr,Integers()!a);num:=5^n+1-Integers()!a;
  for e in Divisors(n) do if e lt n then num-:=e*cl[e];end if;end for;
  assert num mod n eq 0;Append(~cl,num div n);
 end for;
 Append(~rows,[*tup,g,Coefficients(f),cl,Integers()!Evaluate(f,1)*]);print "M2N36_D8_FROB_PREFLIGHT",tup,g,cl,Evaluate(f,1),Cputime(t0);
end for;
Write("C2x36/data/preflight_2_36_D8_Frob_arithmetic.txt",Sprint(rows):Overwrite:=true);print "M2N36_D8_FROB_PREFLIGHT_DONE complete=true",Cputime(t0);quit;
