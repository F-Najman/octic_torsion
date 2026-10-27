/* Exact complete256-class computation for the octic exclusion. */
if not assigned C230Out then C230Out:="C2x30/logs/certified_run_v2";end if;
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
load "C2x30/src/a3_model.m";
G,_,phi:=ClassGroup(A3K:ReductionDivisor:=Divisor(A3c0),Proof:=true);
assert Invariants(G) eq [2,2,2,2,2,2,102,10608,0];
hs:=[phi(o-Degree(o)*Divisor(A3c0)):o in A3orbits];H:=sub<G|hs>;
assert Invariants(H) eq [2,2,2,2,2,8] and #H eq 256;
log:=Open(C230Out cat "/classes256.log","w");
ar:=Open(C230Out cat "/classes256.py","w");
fprintf log,"A3_D8_PROBE_GROUP G=%o H=%o order=%o time=%o\n",Invariants(G),Invariants(H),#H,Cputime(tt);Flush(log);
fprintf ar,"cusp_block_signatures = %o\n",A3signatures;
fprintf ar,"relations = %o\ncusp_block_degrees = %o\ncusp_block_classes = %o\n",[Eltseq(r):r in Rows(RelationMatrix(G))],[Degree(o):o in A3orbits],[Eltseq(h):h in hs];
// Track shortest signed cusp expressions for every class, with no native
// inverse divisor/function representation.
lift:=AssociativeArray();lift[G!0]:=[0:i in [1..#hs]];queue:=[G!0];pos:=1;
while pos le #queue do
 a:=queue[pos];v:=lift[a];
 for j in [1..#hs] do
  for eps in [-1,1] do
   aa:=a+eps*hs[j];
   if not IsDefined(lift,aa) then
    w:=v;w[j]+:=eps;lift[aa]:=w;Append(~queue,aa);
   end if;
  end for;
 end for;
 pos+:=1;
end while;
assert #queue eq 256;
effective:=AssociativeArray();effective[0]:=AssociativeArray();effective[0][G!0]:=[0:i in [1..#hs]];
for d in [1..8] do effective[d]:=AssociativeArray();end for;
for j in [1..#hs] do
 w:=Degree(A3orbits[j]);
 for d in [w..8] do
  for a in Keys(effective[d-w]) do
   aa:=a+hs[j];v:=effective[d-w][a];v[j]+:=1;
   if not IsDefined(effective[d],aa) then effective[d][aa]:=v;end if;
  end for;
 end for;
end for;
fprintf log,"A3_D8_EFFECTIVE_CUSP_CLASSES counts=%o time=%o\n",[#Keys(effective[d]):d in [0..8]],Cputime(tt);Flush(log);
fprintf ar,"cases = [\n";dims:=AssociativeArray();
for idx in [1..#queue] do
 a:=queue[idx];known:=IsDefined(effective[8],a);
 if known then
  v:=effective[8][a];D:=&+[v[j]*A3orbits[j]:j in [1..#hs]];assert Degree(D) eq 8;
 else
  v:=lift[a];D:=8*Divisor(A3c0)+&+[v[j]*(A3orbits[j]-Degree(A3orbits[j])*Divisor(A3c0)):j in [1..#hs]];
 end if;
 assert Degree(D) eq 8 and phi(D-8*Divisor(A3c0)) eq a;
 B:=Basis(D:Reduction:=false,Simplification:="None");h:=#B;
 if IsDefined(dims,h) then dims[h]+:=1;else dims[h]:=1;end if;
 fprintf ar,"[%o,%o,%o,%o,%o],\n",idx,Eltseq(a),known select 1 else 0,v,h;Flush(ar);
 fprintf log,"A3_D8_RR_PROBE index=%o known_effective_cusp=%o dimension=%o time=%o\n",idx,known,h,Cputime(tt);Flush(log);
end for;
fprintf ar,"]\ncomplete = True\nexclusion_claimed = False\n";delete ar;
fprintf log,"A3_D8_RR_PROBE_DONE classes=256 dimensions=%o counts=%o complete=true exclusion_claimed=false time=%o\n",Sort(Setseq(Keys(dims))),[dims[h]:h in Sort(Setseq(Keys(dims)))],Cputime(tt);delete log;
print "A3_D8_RR_PROBE_DONE",dims,Cputime(tt);quit;
