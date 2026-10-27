/* Bounded extension after low-degree function-field enumeration succeeds. */
tp48b:=Cputime();
op48b:=Open("C48/logs/probe_source_places_5_7.log","w");
procedure Probe48B(s)
 fprintf op48b,"%o\n",s;Flush(op48b);print s;
end procedure;
p5f48:=Places(F48,5);assert #p5f48 eq 656;
Probe48B(Sprintf("C48_SOURCE_FUNCTION_FIELD_PLACES degree=5 number=%o time=%o",#p5f48,Cputime(tp48b)));
p6f48:=Places(F48,6);assert #p6f48 eq 2712;
Probe48B(Sprintf("C48_SOURCE_FUNCTION_FIELD_PLACES degree=6 number=%o time=%o",#p6f48,Cputime(tp48b)));
p7f48:=Places(F48,7);assert #p7f48 eq 10784;
Probe48B(Sprintf("C48_SOURCE_FUNCTION_FIELD_PLACES degree=7 number=%o time=%o",#p7f48,Cputime(tp48b)));
