/* Three rational2-isogeny maps to the SAME genus5 X1(21).
   All288 inputs are the complete two-small-quotient survivor set. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
if not assigned Shard then Shard:=0;end if;
if Type(Shard) eq MonStgElt then Shard:=StringToInteger(Shard);end if;assert Shard in [0..3];
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/class_rational_reduction.m");
load "common/pair_place_linear.m";
Z:=MDX1(21,GF(5));C:=Curve(Z);FF:=FunctionField(C);cs:=Cusps(Z);os:=CuspOrbitsQ(Z);
c0:=Support([o:o in os|Degree(o) eq 1][1])[1];
G,_,m2:=ClassGroupRationalReduction(C,c0:Proof:=true);
hg:=[m2(o-Degree(o)*Divisor(c0)):o in os];H:=sub<G|hg>;
assert Invariants(G) eq [2184,0] and Invariants(H) eq [364];
cc:=[m2(Divisor(c)-Degree(c)*Divisor(c0)):c in cs];
assert &and[cc[i] in H:i in [1..#cs]|Degree(cs[i]) le 2];
groupdata:=[*Invariants(G),[Eltseq(r):r in Rows(RelationMatrix(G))],
 [Degree(o):o in os],[Eltseq(a):a in hg],Invariants(H),[Degree(c):c in cs],[Eltseq(a):a in cc]*];
QG,qg:=quo<G|H>;compG:=[*{G!0}*];comps:=[*{QG!0}*];
for dd in [1..4] do
 ss:={G|};
 for j in [1..#cs] do
  ee:=Degree(cs[j]);if ee gt dd then continue;end if;
  ss join:={cc[j]+a:a in compG[dd-ee+1]};
 end for;
 Append(~compG,ss);Append(~comps,{qg(a):a in ss});
end for;
Append(~groupdata,[*Invariants(QG),[Eltseq(r):r in Rows(RelationMatrix(QG))],
 [Eltseq(qg(G.i)):i in [1..Ngens(G)]],[[Eltseq(a):a in ss]:ss in compG]*]);
input:=[* *];
for sh in [0..3] do
 a:=eval Read(Sprintf("C2x42/data/filter_2_42_D8_small_quotients_shard%o_arithmetic.txt",sh));
 input cat:=a[6];
end for;assert #input eq 288;
R:=PolynomialRing(GF(5));fieldobjects:=AssociativeArray();masters:=AssociativeArray();own:=AssociativeArray();fields:=[* *];
for D in [4..8] do
 a:=eval Read(Sprintf("C2x42/data/M2N42_D8_direct_full%o_arithmetic.txt",D));
 assert a[1] eq D and a[3];fieldobjects[D]:=ext<GF(5)|R!a[2]>;Append(~fields,[*D,a[2]*]);
 masters[D]:=a[9];own[D]:=[row:row in a[9]|row[1] eq D];
end for;
aa:=eval Read("C2x42/data/M2N42_D8_direct_small_all_points.txt");assert #aa eq 1 and aa[1][1] eq 4;
all4:=[[*4,[Integers()!a:a in row[1]] cat [Integers()!a:a in row[2]]*]:row in aa[1][3]];
function SourceCurve(D,key)
 K:=fieldobjects[D];b:=Evaluate(R!key[1..D],K.1);c:=Evaluate(R!key[D+1..2*D],K.1);
 E:=EllipticCurve([K|0,c,0,b*(1-b-c),0]);P:=E![0,0];Q:=E![b,b];
 assert 42*Q eq E!0 and 21*Q ne P and 21*Q ne E!0 and 14*Q ne E!0 and 6*Q ne E!0;
 return E,P,Q;
end function;
function IsogPoint(E,P,Q,mapno)
 ts:=[P,21*Q,P+21*Q];T:=ts[mapno];K:=BaseRing(E);
 EG,phi:=IsogenyFromKernel(E,PolynomialRing(K).1-T[1]);A:=phi(2*Q);
 assert 21*A eq EG!0 and 7*A ne EG!0 and 3*A ne EG!0;return EG,A;
end function;
// Fresh original-library comparison on every one of the three maps.
row:=own[4][1];E,P,Q:=SourceCurve(4,row[2]);
for mapno in [1..3] do
 EG,A:=IsogPoint(E,P,Q,mapno);xy:=TateNormalForm_xy(EG,A);
 fast:=C48PairPlaceLinear([FF.1,FF.2],xy);
 ref:=ModuliPoint(Z,EG,rec<recformat<P:PtEll>|P:=A>);assert fast eq ref;
 print "M2N42_ISOG21_NATIVE_REFERENCE",Shard,mapno,Degree(fast),Cputime(t0);
end for;
cache:=AssociativeArray();images:=[* *];tests:=[* *];remaining:=[* *];
procedure Push(~cl,~cache,~images,D,key,mapno)
 ck:=[mapno,D] cat key;
 if IsDefined(cache,ck) then cl:=cache[ck];return;end if;
 E,P,Q:=SourceCurve(D,key);EG,A:=IsogPoint(E,P,Q,mapno);xy:=TateNormalForm_xy(EG,A);
 y:=C48PairPlaceLinear([FF.1,FF.2],xy);e:=Degree(y);assert D mod e eq 0;w:=D div e;
 cl:=w*m2(Divisor(y)-e*Divisor(c0));cache[ck]:=cl;
 Append(~images,[*mapno,D,key,e,w,Eltseq(cl)*]);
end procedure;
for k in [1..#input] do
 if (k-1) mod 4 ne Shard then continue;end if;
 desc:=input[k];sp:=[* *];
 if desc[1] eq 1 then D:=desc[2];row:=masters[D][desc[3]];assert row[1] eq D;Append(~sp,row);
 else
  assert desc[2] eq 4 and desc[4] eq 4;Append(~sp,own[4][desc[3]]);Append(~sp,all4[desc[5]]);
 end if;
 ncdeg:=&+[r[1]:r in sp];assert ncdeg in [4,6,7,8];checks:=[* *];yes:=true;
 for mapno in [1..3] do
  cl:=G!0;
  for row in sp do temp:=G!0;Push(~temp,~cache,~images,row[1],row[2],mapno);cl+:=temp;end for;
  flag:=-qg(cl) in comps[9-ncdeg];Append(~checks,[*mapno,Eltseq(cl),flag select 1 else 0*]);
  if not flag then yes:=false;break;end if;
 end for;
 Append(~tests,[*k,desc,checks,yes select 1 else 0*]);if yes then Append(~remaining,desc);end if;
 if #tests mod 20 eq 0 then print "M2N42_ISOG21_PROGRESS",Shard,#tests,#remaining,Cputime(t0);end if;
end for;
Write(Sprintf("C2x42/data/filter_2_42_D8_three_2isogenies_v2_shard%o_arithmetic.txt",Shard),Sprint([*groupdata,fields,images,tests,remaining*]):Overwrite:=true);
print "M2N42_ISOG21_FILTER_DONE",Shard,#tests,#remaining,"complete=true",Cputime(t0);quit;
