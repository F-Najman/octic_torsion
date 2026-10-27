/* Fresh necessary class tests on all three X1(30) quotients of each
   oriented full2 candidate. Source/cusp complements remain unrestricted.
   Each rough class descent saves an exact two-order divisor identity. */
if not assigned C230Out then C230Out:="C2x30/logs/certified_run_v2";end if;
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
load "common/class_group_unflattened.m";
load "common/native_group_generators.m";
load "common/class_group_direct_RR_checked_v4_fastfb.m";
load "common/pair_af_place_linear.m";
// Exact fibre data are explicitly loaded by the literal-path wrapper.
assert C230FibresComplete and #C230Candidates eq 18 and #C230Points eq 18;
X:=MDX1(30,GF(7));C:=Curve(X);assert Genus(X) eq 9;
os:=CuspOrbitsQ(X);cs:=Cusps(X);c0:=Support([o:o in os|Degree(o) eq 1][1])[1];
F0,mp:=AlgorithmicFunctionField(FunctionField(C));F:=RationalExtensionRepresentation(F0);
c0div:=DivisorGroup(F)!FunctionFieldDivisor(Divisor(c0));
iv:=ClassGroupAbelianInvariants(F0:ReductionDivisor:=c0div,Proof:=true);
assert F`cld`proven;MDProvenClassGroupGenerators(F);
G,_,phi,wh:=MDClassGroupDirectRRCheckedV4FastFB(F0,C230Out cat "/class_identities.log");
gens:=[phi(FunctionFieldDivisor(o-Degree(o)*Divisor(c0))):o in os];
H:=sub<G|gens>;assert Invariants(H) eq [4,8160];GG,pi:=quo<G|H>;
coords:=[F!mp(FunctionField(C).i):i in [1,2]];
log:=Open(C230Out cat "/three_quotients.log","w");
ar:=Open(C230Out cat "/three_quotients.py","w");
fprintf ar,"relations = %o\nrational_cusp_generators = %o\n",[Eltseq(r):r in Rows(RelationMatrix(G))],[Eltseq(g):g in gens];
fprintf ar,"targets = [\n";
fprintf log,"THREE_QUOTIENTS_GROUP G=%o H=%o quotient=%o time=%o\n",Invariants(G),Invariants(H),Invariants(GG),Cputime(tt);Flush(log);

function PairDescriptionFlat(values)
 K:=Parent(values[1]);assert Parent(values[2]) cmpeq K;
 k:=PrimeField(K);p:=#k;e:=LCM([Degree(MinimalPolynomial(v,k)):v in values]);found:=false;
 for bits in [[1,0],[0,1]] cat [[1,a]:a in [1..p-1]] do
  gen:=&+[K!bits[i]*values[i]:i in [1,2]];pol:=MinimalPolynomial(gen,k);
  if Degree(pol) eq e then chosen:=bits;found:=true;break;end if;
 end for;
 assert found;M:=Matrix(k,[Eltseq(gen^i,k):i in [0..e-1]]);assert Rank(M) eq e;
 cc:=[Eltseq(Solution(M,Vector(k,Eltseq(v,k)))):v in values];
 assert &and[&+[K!cc[j][i+1]*gen^i:i in [0..e-1]] eq values[j]:j in [1,2]];
 return [e] cat chosen cat [Integers()!a:a in Coefficients(pol)] cat &cat[[Integers()!a:a in row]:row in cc];
end function;

targets:=[* *];targetindices:=AssociativeArray();sourceclasses:=[* *];targetrecords:=[* *];
for idx in [1..#C230Points] do
 de,bits,pol,cc:=Explode(C230Points[idx]);assert de eq 8;
 K<z>:=ext<GF(7)|Polynomial(GF(7),pol)>;
 bc:=[&+[K!row[i+1]*z^i:i in [0..de-1]]:row in cc];b,c:=Explode(bc);
 E:=EllipticCurve([K!0,c,K!0,b*(1-b-c),K!0]);P2:=E![0,0];Q:=E![b,b];
 assert Order(Q) eq 30 and Order(P2) eq 2 and P2 ne 15*Q;
 qs:=[Q,P2+16*Q,P2+Q];assert &and[Order(q) eq 30:q in qs];
 allclasses:=[];localrecords:=[];
 for variant in [1..3] do
  pair:=TateNormalForm_xy(E,qs[variant]);desc:=PairDescriptionFlat(pair);e:=desc[1];assert 8 mod e eq 0;
  if IsDefined(targetindices,desc) then ti:=targetindices[desc];gc:=targets[ti][2];
  else
   pp:=C48PairAFPlaceLinear(coords,pair);assert Degree(pp) eq e;
   assert Valuation(F!mp(jInvariant(EllipticCurve(X`_E))),pp) ge 0;
   gc:=(8 div e)*phi(Divisor(pp)-e*c0div);
   Append(~targets,[*desc,gc*]);ti:=#targets;targetindices[desc]:=ti;
   fprintf ar,"[%o,%o,%o,%o],\n",ti,desc,Eltseq(gc),pi(gc) eq GG!0 select 1 else 0;Flush(ar);
  end if;
  Append(~allclasses,gc);Append(~localrecords,[variant,ti,pi(gc) eq GG!0 select 1 else 0]);
  fprintf log,"THREE_QUOTIENT_POINT source=%o variant=%o target=%o targetdegree=%o passes=%o class=%o time=%o\n",idx,variant,ti,e,pi(gc) eq GG!0,Eltseq(gc),Cputime(tt);Flush(log);
 end for;
 Append(~sourceclasses,allclasses);Append(~targetrecords,localrecords);
end for;
fprintf ar,"]\nsource_maps = %o\ncandidate_tests = [\n",[x:x in targetrecords];
remaining:=[];
for idx in [1..#C230Candidates] do
 cand:=C230Candidates[idx];assert #cand eq 1 and cand[1][2] eq 1;
 pid:=cand[1][1];passes:=&and[pi(a) eq GG!0:a in sourceclasses[pid]];
 if passes then Append(~remaining,idx);end if;
 fprintf ar,"[%o,%o,%o],\n",idx,pid,passes select 1 else 0;
 fprintf log,"THREE_QUOTIENT_CANDIDATE candidate=%o source=%o passes=%o time=%o\n",idx,pid,passes,Cputime(tt);Flush(log);
end for;
fprintf ar,"]\nremaining = %o\ncomplete = True\nexclusion_claimed = False\n",remaining;delete ar;
fprintf wh,"]\ncomplete = True\n";delete wh;
fprintf log,"THREE_QUOTIENTS_DONE candidates=18 maps=54 distincttargets=%o remaining=%o complete=true exclusion_claimed=false time=%o\n",#targets,#remaining,Cputime(tt);delete log;
print "THREE_QUOTIENTS_DONE",#targets,remaining,Cputime(tt);quit;
