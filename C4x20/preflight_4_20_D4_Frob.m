SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
J:=JH(80,[21]);g:=Dimension(J);print "M4N20_D4_FROB_GENUS",g;
rows:=[* *];orders:=[];
for p in [3,13,17] do
 f:=FrobeniusPolynomial(J,p);assert Degree(f) eq 2*g;
 step:=p eq 3 select 2 else 1;tr:=[];cl:=[];
 for n in [1..4*step] do
  a:=-n*Coefficient(f,2*g-n);
  for j in [1..n-1] do a-:=Coefficient(f,2*g-j)*tr[n-j];end for;
  Append(~tr,Integers()!a);
 end for;
 for n in [1..4] do
  num:=p^(n*step)+1-tr[n*step];
  for e in Divisors(n) do if e lt n then num-:=e*cl[e];end if;end for;
  assert num mod n eq 0;Append(~cl,num div n);
 end for;
 order:=step eq 2 select Evaluate(f,1)*Evaluate(f,-1) else Evaluate(f,1);
 Append(~orders,Integers()!order);Append(~rows,[*p,step,Coefficients(f),cl,Integers()!order*]);
 print "M4N20_D4_FROB_ROW",p,step,cl,order,Cputime(t0);
end for;
Write("C4x20/data/preflight_4_20_D4_Frob_arithmetic.txt",Sprint([*g,rows*]):Overwrite:=true);
print "M4N20_D4_JF_TORSION_ORDER_BOUND",GCD(orders);
print "M4N20_D4_FROB_DONE complete=true",Cputime(t0);quit;
