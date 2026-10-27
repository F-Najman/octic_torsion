/* Independent modular-symbol Frobenius polynomial for the octic count check. */
SetColumns(0);SetQuitOnError(true);
f:=FrobeniusPolynomial(JOne(33),2);assert Degree(f) eq 42;
oo:=Open("C33/data/frobenius_coefficients.txt","w");
fprintf oo,"%o\n",Coefficients(f);delete oo;
print "C33_D8_INDEPENDENT_FROBENIUS_EXPORTED degree=42 complete=true";quit;
