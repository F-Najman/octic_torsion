/* Exact target-automorphism equivalence for SINGLE high-degree supports.
   Source image need not be stable under these maps. The necessary target
   Q-image-plus-cusps predicate is stable. Mixed supports are not reduced. */
function C48M2NormalizedPair(E,P,Q)
 a1,a2,a3,a4,a6:=Explode(aInvariants(E));assert a1 eq 0 and a3 eq 0 and a6 eq 0;
 assert 2*P eq E!0 and P ne E!0 and P[2] eq 0;
 r:=P[1];xx:=Q[1]-r;yy:=Q[2];assert xx ne 0 and yy ne 0;
 b:=xx^3/yy^2;c:=(a2+3*r)*xx^2/yy^2;
 assert (a4+2*a2*r+3*r^2)*xx^4/yy^4 eq b*(1-b-c);
 return [b,c];
end function;
function C48TargetAutKey(bc,dd)
 K:=Parent(bc[1]);assert #K eq 5^dd;b,c:=Explode(bc);
 E0:=EllipticCurve([1-c,-b,-b,K!0,K!0]);P0:=E0![0,0];E2,L2:=C48TwoIsogenousModuli(E0,P0);
 P2:=L2`P;Q:=L2`Q;Q2:=12*Q;assert Q2 ne P2 and Q2 ne E2!0;
 keys:=[];
 for eps,eta in [0,1] do
  for aa in [1,5,7,11] do
   PP:=P2+eps*Q2;QQ:=aa*Q+eta*P2;
   assert 24*QQ eq E2!0 and 12*QQ ne E2!0 and 8*QQ ne E2!0 and PP ne 12*QQ;
   pair:=C48M2NormalizedPair(E2,PP,QQ);
   Append(~keys,MDM2FrobeniusKey(pair,K,5,dd));
  end for;
 end for;
 assert #keys eq 16;return [dd] cat Sort(keys)[1];
end function;
procedure C48BuildHighTargetOrbits()
 tt:=Cputime();seen:=AssociativeArray();representatives:=[];records:=[];representativeKeys:=[];
 log:=Open("C48/logs/clean_run/high_target_aut_orbits.log","w");
 ar:=Open("C48/data/clean_run/high_target_aut_orbits.py","w");
 fprintf ar,"master_sha = \"%o\"\ntarget_automorphisms = [[e,h,a] for e in [0,1] for h in [0,1] for a in [1,5,7,11]]\ncases = [\n",C48MasterSHA;
 for dd in [6..8] do
  before:=#representatives;
  for i in [1..#sourceRepBC48[dd]] do
   key:=C48TargetAutKey(sourceRepBC48[dd][i],dd);
   if not IsDefined(seen,key) then Append(~representatives,[dd,i]);Append(~representativeKeys,key);seen[key]:=#representatives;end if;
   cid:=seen[key];Append(~records,[dd,i,cid]);
   fprintf ar,"[%o,%o,%o,%o],\n",dd,i,cid,key;
   if i mod 500 eq 0 then fprintf log,"C48_TARGET_AUT_PROGRESS degree=%o input=%o total=%o clusters=%o time=%o\n",dd,i,#sourceRepBC48[dd],#representatives-before,Cputime(tt);Flush(log);end if;
  end for;
  fprintf log,"C48_TARGET_AUT_TYPE degree=%o input=%o clusters=%o time=%o\n",dd,#sourceRepBC48[dd],#representatives-before,Cputime(tt);Flush(log);
 end for;
 assert #records eq &+[#sourceRepBC48[d]:d in [6..8]];
 assert Seqset([r[3]:r in records]) eq {1..#representatives};
 fprintf ar,"]\nrepresentatives = %o\ncomplete = True\n",representatives;delete ar;
 fprintf log,"C48_TARGET_AUT_DONE inputs=%o clusters=%o complete=true time=%o\n",#records,#representatives,Cputime(tt);delete log;
 targetOrbitRepresentatives48:=representatives;
 Write("C48/data/clean_run/high_target_aut_orbits.m",Sprintf("C48TargetOrbitRepresentatives:=%o;C48TargetOrbitKeys:=%o;C48TargetOrbitLabels:=%o;C48TargetOrbitMasterSHA:=\"%o\";",representatives,representativeKeys,records,C48MasterSHA));
 print "C48_TARGET_AUT_DONE",#records,#representatives,Cputime(tt);
end procedure;
C48BuildHighTargetOrbits();
