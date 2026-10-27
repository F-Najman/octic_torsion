/* Birationality of the archived source model and the universal Tate model. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
X:=MDX1(48,GF(5));F:=FunctionField(Curve(X));b,c:=Explode([F!z:z in X`_coordinates]);
E:=EllipticCurve([1-c,-b,-b,F!0,F!0]);P:=E![0,0];
xy:=TateNormalForm_xy(E,P);assert xy eq [F.1,F.2];
printf "C48_GENERIC_SOURCE_TATE_INVERSE_DONE exact_function_field_identity=true time=%o\n",Cputime(tt);
quit;
