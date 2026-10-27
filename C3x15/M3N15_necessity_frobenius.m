SetColumns(0);SetQuitOnError(true);t0:=Cputime();
A:=JH(45,[16]);assert Dimension(A) eq 9;
rows:=[* *];
for p in [2,7,13,19] do
 f:=FrobeniusPolynomial(A,p);assert Degree(f) eq 18;
 order:=p eq 2 select Evaluate(f,1)*Evaluate(f,-1) else Evaluate(f,1);
 Append(~rows,[*p,Coefficients(f),Integers()!order*]);
 print "M3N15_NECESSITY_FROBENIUS",p,order,Cputime(t0);
end for;
Write("C3x15/data/M3N15_necessity_frobenius_arithmetic.txt",Sprint(rows):Overwrite:=true);
print "M3N15_NECESSITY_FROBENIUS_DONE",Cputime(t0);quit;
