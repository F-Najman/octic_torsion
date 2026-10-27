/* Completes the all-cuspidal reduction branch of the fresh quartic lemma.
   We do not infer finite-field gonality from rational gonality. */
if not assigned C230Out then C230Out:="C2x30/logs/certified_run_v2";end if;
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
X:=MDX1(30,GF(7));cs:=Cusps(X);C:=Curve(X);
assert Genus(X) eq 9;F,mp:=AlgorithmicFunctionField(FunctionField(C));
log:=Open(C230Out cat "/quartic_allcusps.log","w");
ar:=Open(C230Out cat "/quartic_allcusps.py","w");
fprintf ar,"cusp_degrees = %o\ntests = [\n",[Degree(c):c in cs];
count:=0;exceptions:=0;
procedure TestCusps(j,left,v,~count,~exceptions)
 if j gt #cs then
  if left ne 0 then return;end if;
  D:=&+[v[i]*Divisor(cs[i]):i in [1..#cs]];assert Degree(D) eq 4;
  DD:=FunctionFieldDivisor(D);B:=Basis(DD:Reduction:=false,Simplification:="None");
  count+:=1;if #B ne 1 then exceptions+:=1;end if;
  fprintf ar,"[%o,%o],\n",v,#B;
  return;
 end if;
 for m in [0..left div Degree(cs[j])] do TestCusps(j+1,left-m*Degree(cs[j]),v cat [m],~count,~exceptions);end for;
end procedure;
TestCusps(1,4,[],~count,~exceptions);
fprintf ar,"]\ncomplete = True\n";delete ar;
fprintf log,"X130_D4_ALLCUSPS_DONE cuspdegrees=%o tested=%o exceptions=%o complete=true time=%o\n",[Degree(c):c in cs],count,exceptions,Cputime(tt);delete log;
assert exceptions eq 0;print "X130_D4_ALLCUSPS_DONE",count,exceptions,Cputime(tt);quit;
