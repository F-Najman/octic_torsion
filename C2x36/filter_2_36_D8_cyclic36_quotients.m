/* Three rational cyclic36 quotients on the461 exact first-stage survivors.
   All target cusp complements are retained. No exclusion assumed. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
Attach("common/fast_hecke_cyclic_commonzero.m");
Attach("common/class_rational_reduction.m");
Y:=MDX1(36,GF(5));C:=Curve(Y);cs:=Cusps(Y);os:=CuspOrbitsQ(Y);
assert Genus(Y) eq 17 and Sort([Degree(o):o in os]) eq [1:i in [1..9]] cat [2:i in [1..6]] cat [3,4,6,6];
c0:=Support([o:o in os|Degree(o) eq 1][1])[1];
print "M2N36_CYCLIC36_SETUP_READY",Cputime(t0);
G,m1,m2:=ClassGroupRationalReduction(C,c0:Proof:=true);
assert &*[a:a in Invariants(G)|a ne 0] eq 5061852020736;
hg:=[m2(o-Degree(o)*Divisor(c0)):o in os];H:=sub<G|hg>;
assert Invariants(H) eq [12,252,7812];QH,qh:=quo<G|H>;
cc:=[m2(Divisor(c)-Divisor(c0)):c in cs|Degree(c) eq 1];
cq:=[m2(Divisor(c)-2*Divisor(c0)):c in cs|Degree(c) eq 2];
comps:=[* {QH!0},{qh(a):a in cc},{qh(a):a in cq} join {qh(a+b):a,b in cc} *];
groupdata:=[*Invariants(G),[Eltseq(r):r in Rows(RelationMatrix(G))],
 [Degree(o):o in os],[Eltseq(v):v in hg],Invariants(H),
 [[Eltseq(v):v in comp]:comp in comps]*];
print "M2N36_CYCLIC36_GROUP_READY",Invariants(G),Invariants(H),[#s:s in comps],Cputime(t0);
data:=eval Read("C2x36/data/filter_2_36_D8_quotient_arithmetic.txt");
singles:=data[3];small:=data[4];survivors:=data[6];assert #survivors eq 461;
fields:=data[2];own:=AssociativeArray();sm:=AssociativeArray();
for d in [3..8] do own[d]:=[r:r in singles|r[1] eq d];end for;
for d in [3..5] do sm[d]:=[r:r in small|r[1] eq d];end for;
R:=PolynomialRing(GF(5));cache:=AssociativeArray();images:=[* *];tests:=[* *];remaining:=[* *];
procedure Push(~answer,~cache,~images,row,mapno,Y,c0,m2)
 D:=row[1];ky:=row[3];cachekey:=[mapno,D] cat ky;
 if IsDefined(cache,cachekey) then answer:=cache[cachekey];return;end if;
 def:=[a[2]:a in fields|a[1] eq D][1];K<z>:=ext<GF(5)|R!def>;
 b:=Evaluate(R!ky[1..D],z);c:=Evaluate(R!ky[D+1..2*D],z);
 E:=EllipticCurve([K|0,c,0,b*(1-b-c),0]);P:=E![0,0];Q:=E![b,b];
 assert 36*Q eq E!0 and 18*Q ne P and 18*Q ne E!0 and 12*Q ne E!0;
 if mapno eq 1 then EG:=E;QG:=Q;
 else
  T:=mapno eq 2 select P else P+18*Q;
  EG,phi:=IsogenyFromKernel(E,PolynomialRing(K).1-T[1]);QG:=phi(Q);
  assert 36*QG eq EG!0 and 18*QG ne EG!0 and 12*QG ne EG!0;
 end if;
 y:=ModuliPointCyclicCommonZeros(Y,EG,rec<recformat<P:PtEll>|P:=QG>);
 e:=Degree(y);assert D mod e eq 0;w:=D div e;
 answer:=w*m2(Divisor(y)-e*Divisor(c0));cache[cachekey]:=answer;
 Append(~images,[*cachekey,e,w,Eltseq(answer)*]);
end procedure;
for k in [1..#survivors] do
 desc:=survivors[k];sp:=[* *];
 if desc[1] eq 1 then
  d:=desc[2];rr:=[r:r in own[d]|r[2] eq desc[3]];assert #rr eq 1;Append(~sp,rr[1]);
 else
  Append(~sp,own[desc[2]][desc[3]]);Append(~sp,sm[desc[4]][desc[5]]);
 end if;
 ncdeg:=&+[r[1]:r in sp];cd:=8-ncdeg;assert cd in [0,1,2];passed:=true;checks:=[* *];
 for mapno in [1..3] do
  cl:=G!0;
  for r in sp do a:=G!0;Push(~a,~cache,~images,r,mapno,Y,c0,m2);cl+:=a;end for;
  hits:=[Eltseq(a):a in comps[cd+1]|qh(cl)+a eq QH!0];
  Append(~checks,[*mapno,Eltseq(cl),hits*]);
  if #hits eq 0 then passed:=false;break;end if;
 end for;
 Append(~tests,[*k,desc,cd,checks,passed select 1 else 0*]);
 if passed then Append(~remaining,desc);end if;
 if k mod 20 eq 0 then
  Write("C2x36/data/filter_2_36_D8_cyclic36_arithmetic.txt",Sprint([*groupdata,images,tests,remaining*]):Overwrite:=true);
  print "M2N36_CYCLIC36_PROGRESS",k,#remaining,Cputime(t0);
 end if;
end for;
Write("C2x36/data/filter_2_36_D8_cyclic36_arithmetic.txt",Sprint([*groupdata,images,tests,remaining*]):Overwrite:=true);
print "M2N36_CYCLIC36_FILTER_DONE",#tests,#remaining,"complete=true",Cputime(t0);quit;
