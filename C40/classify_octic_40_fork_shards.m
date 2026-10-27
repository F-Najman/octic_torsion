/* Run in the retained, PROVEN C40 state only after export_fixed_octic_inputs.m.
   Fork copies the same exact field, group basis, cusp data, and master places.
   No independent random representative enumeration is used in any shard.
   Each child exits after its 24 tests; the original parent remains at its
   prompt after its own 24 tests so that the controlling terminal stays open. */
load "C40/data/fixed_octic_inputs.m";
assert C40FixedFieldPolynomial eq fmaster40 and C40FixedTateBC eq fixedbc40;
assert C40FixedOriginalPlaceIndices eq fixedplaceindices40;
assert C40FixedAllOrbitIndices eq fixedorbits40;
assert assigned C40FixedMasterSHA and #C40FixedMasterSHA eq 64;
assert #reps40[8] eq 96 and #places[8] eq 708;
forkProvenField40:=RationalExtensionRepresentation(AlgorithmicFunctionField(FunctionField(C)));
assert assigned forkProvenField40`cld and forkProvenField40`cld`proven;
assert &*[n:n in Invariants(G)|n ne 0] eq 31579088486400;
load "C40/export_factorbase_geometry.m";

// Finish the cheap mixed support tests once, before copying the snapshot.
pairout40:=Open("C40/data/quartic_pair_classes_fork.py","w");
fprintf pairout40,"quartic_pair_tests = [\n";
paircount40:=0;pairpasses40:=0;paircandidates40:={G|};
for P40 in reps40[4] do
 first40:=Index(places[4],P40);
 for second40 in [1..32] do
  a40:=quartic40classes[first40]+quartic40classes[second40];member40:=a40 in HQ;
  paircount40+:=1;
  if member40 then pairpasses40+:=1;Include(~paircandidates40,a40);end if;
  fprintf pairout40,"[%o,%o,%o,%o],\n",first40,second40,Eltseq(a40),member40 select 1 else 0;
 end for;
end for;
assert paircount40 eq 128;
fprintf pairout40,"]\ncomplete = True\n";delete pairout40;
print "C40_QUARTIC_PAIR_FILTER_DONE",paircount40,pairpasses40,#paircandidates40;

// Close the inherited setup handle: all subsequent evidence is shard-local.
delete out40;
load "common/class_group_direct_RR_checked_v4.m";
snapshotG40:=G;snapshotRelations40:=RelationMatrix(G);
/* Compile the complete worker body before Fork: child processes do not\n   race while reading further source from a shared file descriptor. */
procedure C40RunForkShards()
forkChild40:=false;forkShard40:=3;forkPIDs40:=[];
for jj40 in [0..2] do
 pid40:=Fork();
 if pid40 eq 0 then
  forkChild40:=true;forkShard40:=jj40;break;
 end if;
 assert pid40 gt 0;Append(~forkPIDs40,pid40);
 print "C40_FORK_CHILD_CREATED",jj40,pid40;
end for;
if not forkChild40 then
 pidout40:=Open("C40/data/fork_processes.py","w");
 fprintf pidout40,"children = %o\n",forkPIDs40;delete pidout40;
end if;
forkIndices40:=[i:i in [1..96]|(i-1) mod 4 eq forkShard40];
assert #forkIndices40 eq 24;
forkPrefix40:=Sprintf("C40/logs/fork_shard%o",forkShard40);
forkLog40:=Open(forkPrefix40 cat ".log","w");
fprintf forkLog40,"C40_FORK_SHARD_STARTED shard=%o indices=%o master_sha=%o group_order=31579088486400\n",forkShard40,forkIndices40,C40FixedMasterSHA;Flush(forkLog40);
SetSeed(400800+forkShard40);
Gfork40,mrawfork40,phifork40,witnessfork40:=MDClassGroupDirectRRCheckedV4(AlgorithmicFunctionField(FunctionField(C)),forkPrefix40 cat "_identities.log");
assert Gfork40 cmpeq snapshotG40 and RelationMatrix(Gfork40) eq snapshotRelations40;
forkMap40:=map<DivisorGroup(C)->Gfork40|DD:->phifork40(FunctionFieldDivisor(DD))>;
forkAr40:=Open(Sprintf("C40/data/fork_shard%o_arithmetic.py",forkShard40),"w");
fprintf forkAr40,"shard = %o\nmaster_sha = \"%o\"\nexpected_indices = %o\nrelations = %o\nrational_cusp_generators = %o\ntests = [\n",forkShard40,C40FixedMasterSHA,forkIndices40,[Eltseq(r):r in Rows(RelationMatrix(Gfork40))],[Eltseq(a):a in qc];
forkDone40:=[];forkPass40:=0;
for idx40 in forkIndices40 do
 fprintf forkLog40,"C40_FORK_CASE_STARTED shard=%o index=%o/96 time=%o\n",forkShard40,idx40,Cputime();Flush(forkLog40);
 aa40:=forkMap40(Divisor(reps40[8][idx40])-8*Divisor(c0));
 member40:=aa40 in HQ;if member40 then forkPass40+:=1;end if;
 fprintf forkAr40,"[%o,%o,%o],\n",idx40,Eltseq(aa40),member40 select 1 else 0;Flush(forkAr40);
 Append(~forkDone40,idx40);
 fprintf forkLog40,"C40_FORK_CASE_DONE shard=%o index=%o/96 in_rational_image=%o time=%o\n",forkShard40,idx40,member40,Cputime();Flush(forkLog40);
end for;
assert forkDone40 eq forkIndices40;
fprintf forkAr40,"]\ncomplete = True\n";delete forkAr40;
fprintf witnessfork40,"]\n";delete witnessfork40;
fprintf forkLog40,"C40_FORK_SHARD_DONE shard=%o tested=24 admitted=%o complete=true exact_indices=%o master_sha=%o\n",forkShard40,forkPass40,forkDone40,C40FixedMasterSHA;delete forkLog40;
print "C40_FORK_SHARD_DONE",forkShard40,24,forkPass40;
if forkChild40 then quit;end if;
print "C40_FORK_PARENT_RETAINED_AT_PROMPT children=",forkPIDs40;

end procedure;
C40RunForkShards();
