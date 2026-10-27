/* Bounded diagnostic in the retained target-class-group session. */
tp48:=Cputime();
op48:=Open("C48/logs/probe_source_places.log","w");
procedure Probe48(s)
 fprintf op48,"%o\n",s;Flush(op48);print s;
end procedure;
X48:=MDX1(48,GF(5));CX48:=Curve(X48);
F48:=AlgorithmicFunctionField(FunctionField(CX48));
Probe48(Sprintf("C48_SOURCE_FUNCTION_FIELD_READY genus=%o AFdegree=%o time=%o",Genus(F48),Degree(F48),Cputime(tp48)));
p3f48:=Places(F48,3);assert #p3f48 eq 16;
Probe48(Sprintf("C48_SOURCE_FUNCTION_FIELD_PLACES degree=3 number=%o time=%o",#p3f48,Cputime(tp48)));
p4f48:=Places(F48,4);assert #p4f48 eq 68;
Probe48(Sprintf("C48_SOURCE_FUNCTION_FIELD_PLACES degree=4 number=%o time=%o",#p4f48,Cputime(tp48)));
