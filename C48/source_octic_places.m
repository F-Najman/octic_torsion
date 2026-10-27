/* Fresh source octic place enumeration through the normalized function field.
   Earlier CurveFF enumeration failed before completing any support branch.
   Exact cardinality is independently predicted by the Frobenius polynomial. */
te48:=Cputime();
oe48:=Open("C48/logs/source_octic_places.log","w");
fprintf oe48,"C48_OCTIC_ENUM_START expected=48992 method=function_field time=%o\n",Cputime(te48);Flush(oe48);
p8f48:=Places(F48,8);assert #p8f48 eq 48992;
fprintf oe48,"C48_OCTIC_ENUM_DONE count=%o complete=true time=%o\n",#p8f48,Cputime(te48);Flush(oe48);
printf "C48_OCTIC_ENUM_DONE count=%o complete=true time=%o\n",#p8f48,Cputime(te48);
