SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");rows:=[* *];
for tup in [[75,1]] do
 J:=JH(tup[1],[tup[2]]);g:=Dimension(J);f:=FrobeniusPolynomial(J,2);assert Degree(f) eq 2*g;tr:=[];cl:=[];
 for n in [1..8] do
  a:=-n*Coefficient(f,2*g-n);for j in [1..n-1] do a-:=Coefficient(f,2*g-j)*tr[n-j];end for;
  Append(~tr,Integers()!a);num:=2^n+1-Integers()!a;
  for e in Divisors(n) do if e lt n then num-:=e*cl[e];end if;end for;
  assert num mod n eq 0;Append(~cl,num div n);
 end for;
 Append(~rows,[*tup,g,Coefficients(f),cl,Integers()!Evaluate(f,1)*]);print "C75_SOURCE_FROB_PREFLIGHT",tup,g,cl,Evaluate(f,1),Cputime(t0);
end for;
Write("C75/data/C75_source_frob_root_arithmetic.txt",Sprint(rows):Overwrite:=true);print "C75_SOURCE_FROB_PREFLIGHT_DONE complete=true",Cputime(t0);quit;
