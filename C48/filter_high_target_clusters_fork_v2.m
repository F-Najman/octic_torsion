SetQuitOnError(true);
load "common/class_group_direct_RR_checked_v4_fastfb.m";
/* Fixed-source, fixed-target-orbit high single-support class tests.
   Compile the complete procedure BEFORE Fork. Each worker owns its output
   handles and exact RR witnesses; the clean proved presentation is shared. */
load "C48/data/clean_run/high_target_aut_orbits.m";
assert C48TargetOrbitMasterSHA eq C48MasterSHA;
assert #C48TargetOrbitRepresentatives eq #C48TargetOrbitKeys;
assert FR`cld`proven and FR`cld`Cl cmpeq G;
snapshotG48:=G;snapshotRelations48:=RelationMatrix(G);
procedure C48RunHighClusterShards()
 child:=false;shard:=7;pids:=[];
 for j in [0..6] do
  pid:=Fork();
  if pid eq 0 then child:=true;shard:=j;break;end if;
  assert pid gt 0;Append(~pids,pid);print "C48_HIGH_FORK_CHILD",j,pid;
 end for;
 if not child then
  po:=Open("C48/data/clean_run_v2/high_fork_processes.py","w");fprintf po,"children = %o\n",pids;delete po;
 end if;
 ids:=[i:i in [1..#C48TargetOrbitRepresentatives]|(i-1) mod 8 eq shard];
 prefix:=Sprintf("C48/logs/clean_run_v2/high_shard%o",shard);
 log:=Open(prefix cat ".log","w");ar:=Open(Sprintf("C48/data/clean_run_v2/high_shard%o_arithmetic.py",shard),"w");
 fprintf log,"C48_HIGH_SHARD_START shard=%o master_sha=%o exact_cluster_ids=%o\n",shard,C48MasterSHA,ids;Flush(log);
 SetSeed(480800+shard);
 GC,_,phic,witnessc:=MDClassGroupDirectRRCheckedV4FastFB(F0,prefix cat "_identities.log");
 assert GC cmpeq snapshotG48 and RelationMatrix(GC) eq snapshotRelations48;
 fprintf ar,"shard = %o\nmaster_sha = \"%o\"\nexpected_cluster_ids = %o\nrelations = %o\nrational_cusp_generators = %o\ntests = [\n",shard,C48MasterSHA,ids,[Eltseq(r):r in Rows(RelationMatrix(GC))],[Eltseq(a):a in gens];Flush(ar);
 done:=[];admitted:=0;tt:=Cputime();
 for cid in ids do
  d,idx:=Explode(C48TargetOrbitRepresentatives[cid]);key:=C48TargetOrbitKeys[cid];
  assert key[1] eq d and #key eq 1+2*d;K:=sourceFields48[d];
  b:=&+[K!key[2+i]*K.1^i:i in [0..d-1]];
  c:=&+[K!key[2+d+i]*K.1^i:i in [0..d-1]];
  // The chosen point is an actual target level automorph of the source
  // image. The predicate is invariant under all16 such automorphisms.
  EE:=EllipticCurve([K!0,c,K!0,b*(1-b-c),K!0]);P2:=EE![0,0];Q:=EE![b,b];
  assert 24*Q eq EE!0 and 12*Q ne EE!0 and 8*Q ne EE!0 and 12*Q ne P2;
  yy:=C48PairPlaceLinear(Y`_coordinates,[b,c]);e:=Degree(yy);
  assert d mod e eq 0 and d div e in [1,2];
  fprintf log,"C48_HIGH_CASE_START cluster=%o degree=%o original_index=%o target_degree=%o time=%o\n",cid,d,idx,e,Cputime(tt);Flush(log);
  a:=(d div e)*phic(FunctionFieldDivisor(Divisor(yy)-e*Divisor(c0)));
  pass:=-pi(a) in comps48[8-d];if pass then admitted+:=1;end if;
  fprintf ar,"[%o,%o,%o,%o,%o,%o],\n",cid,d,idx,e,Eltseq(a),pass select 1 else 0;Flush(ar);
  Append(~done,cid);
  fprintf log,"C48_HIGH_CASE_DONE cluster=%o degree=%o original_index=%o admitted=%o time=%o\n",cid,d,idx,pass,Cputime(tt);Flush(log);
  if #done mod 50 eq 0 then print "C48_HIGH_SHARD_PROGRESS",shard,#done,#ids,admitted,Cputime(tt);end if;
 end for;
 assert done eq ids;
 fprintf ar,"]\ncomplete = True\n";delete ar;
 fprintf witnessc,"]\n";delete witnessc;
 fprintf log,"C48_HIGH_SHARD_DONE shard=%o tested=%o admitted=%o complete=true exact_ids=%o time=%o\n",shard,#done,admitted,done,Cputime(tt);delete log;
 print "C48_HIGH_SHARD_DONE",shard,#done,admitted,Cputime(tt);
 if child then quit;end if;
 print "C48_CLEAN_PARENT_RETAINED children",pids;
end procedure;
C48RunHighClusterShards();
