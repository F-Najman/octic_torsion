SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);
AttachSpec("mdmagma/v2/mdmagma.spec");
X:=MDX1(34,GF(3));C:=Curve(X);
load "C34/known_34_reductions_v2.m";
quit;
