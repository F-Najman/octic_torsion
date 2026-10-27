/* Fresh necessary-condition filters through X1(2,14) and X1(21).
   Their full Q-Jacobians are rational cuspidal. No source-generation claim. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
if not assigned Shard then Shard:=0;end if;
if Type(Shard) eq MonStgElt then Shard:=StringToInteger(Shard);end if;assert Shard in [0..3];
OutputFile:=Sprintf("C2x42/data/filter_2_42_D8_small_quotients_shard%o_arithmetic.txt",Shard);
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
Attach("common/fast_hecke_cyclic_commonzero.m");
Attach("common/class_rational_reduction.m");
load "common/pair_place_linear.m";
Y:=MDX11(2,14,GF(5));CY:=Curve(Y);csy:=Cusps(Y);cy:=Rep({c:c in csy|Degree(c) eq 1});
GY,_,my:=ClassGroupRationalReduction(CY,cy:Proof:=true);
gy:=[my(Divisor(c)-Degree(c)*Divisor(cy)):c in csy];HY:=sub<GY|gy>;
assert Invariants(GY) eq [2,2,12,36,0] and Invariants(HY) eq [2,2,6,18];
Z:=MDX1(21,GF(5));CZ:=Curve(Z);csz:=Cusps(Z);osz:=CuspOrbitsQ(Z);
cz:=Support([o:o in osz|Degree(o) eq 1][1])[1];
GZ,_,mz:=ClassGroupRationalReduction(CZ,cz:Proof:=true);
gz:=[mz(o-Degree(o)*Divisor(cz)):o in osz];HZ:=sub<GZ|gz>;
assert Invariants(GZ) eq [2184,0] and Invariants(HZ) eq [364];
QZ,qz:=quo<GZ|HZ>;
ccz:=[mz(Divisor(c)-Degree(c)*Divisor(cz)):c in csz];
comps:=[*{QZ!0}*];
for d in [1..4] do
 s:={QZ|};
 for j in [1..#csz] do
  e:=Degree(csz[j]);if e gt d then continue;end if;
  s join:={qz(ccz[j])+a:a in comps[d-e+1]};
 end for;Append(~comps,s);
end for;
groupdata:=[*[*Invariants(GY),[Eltseq(r):r in Rows(RelationMatrix(GY))],[Degree(c):c in csy],[Eltseq(a):a in gy],Invariants(HY)*],
 [*Invariants(GZ),[Eltseq(r):r in Rows(RelationMatrix(GZ))],[Degree(o):o in osz],[Eltseq(a):a in gz],Invariants(HZ),[[Eltseq(a):a in s]:s in comps]*]*];
print "M2N42_SMALL_TARGETS_READY",Invariants(GY),Invariants(HY),Invariants(GZ),Invariants(HZ),[#s:s in comps],Cputime(t0);
R:=PolynomialRing(GF(5));fieldobjects:=AssociativeArray();fields:=[* *];
for D in [4..8] do
 a:=eval Read(Sprintf("C2x42/data/M2N42_D8_direct_full%o_arithmetic.txt",D));
 assert a[1] eq D and a[3];fieldobjects[D]:=ext<GF(5)|R!a[2]>;Append(~fields,[*D,a[2]*]);
end for;
single:=[* *];small:=[* *];pairs:=[* *];survivors:=[* *];own:=AssociativeArray();
cachey:=AssociativeArray();cachez:=AssociativeArray();
function Key(values,D)
 keys:=[];for i in [1..D] do Append(~keys,&cat[[Integers()!x:x in Eltseq(a,GF(5))]:a in values]);values:=[a^5:a in values];end for;return Sort(keys)[1];
end function;
procedure Push(~ay,~az,~ed,~cachey,~cachez,D,key)
 K:=fieldobjects[D];zz:=K.1;b:=Evaluate(R!key[1..D],zz);c:=Evaluate(R!key[D+1..2*D],zz);
 E:=EllipticCurve([K|0,c,0,b*(1-b-c),0]);P:=E![0,0];Q:=E![b,b];
 assert 42*Q eq E!0 and 21*Q ne E!0 and 14*Q ne E!0 and 6*Q ne E!0 and 21*Q ne P;
 A:=3*Q;bc:=[A[1]^3/A[2]^2,c*A[1]^2/A[2]^2];ky:=[D] cat Key(bc,D);
 if IsDefined(cachey,ky) then ay:=cachey[ky][1];ey:=cachey[ky][2];
 else
  y:=C48PairPlaceLinear(Y`_coordinates,bc);ey:=Degree(y);assert D mod ey eq 0;
  ay:=(D div ey)*my(Divisor(y)-ey*Divisor(cy));cachey[ky]:=<ay,ey>;
 end if;
 A:=2*Q;xy:=TateNormalForm_xy(E,A);kz:=[D] cat Key(xy,D);
 if IsDefined(cachez,kz) then az:=cachez[kz][1];ez:=cachez[kz][2];
 else
  z:=C48PairPlaceLinear([FunctionField(CZ).1,FunctionField(CZ).2],xy);ez:=Degree(z);assert D mod ez eq 0;
  az:=(D div ez)*mz(Divisor(z)-ez*Divisor(cz));cachez[kz]:=<az,ez>;
 end if;
 ed:=[ey,D div ey,ez,D div ez];
end procedure;
for D in [4..8] do
 if D lt 8 and Shard ne 0 then continue;end if;
 a:=eval Read(Sprintf("C2x42/data/M2N42_D8_direct_full%o_arithmetic.txt",D));assert a[1] eq D and a[3];
 def:=a[2];assert fields[D-3][2] eq def;records:=a[9];own[D]:=[* *];passed:=0;
 for i in [1..#records] do
  d,key:=Explode(records[i]);if d ne D then continue;end if;
  if D eq 8 and (i-1) mod 4 ne Shard then continue;end if;
  ay:=GY!0;az:=GZ!0;ed:=[];Push(~ay,~az,~ed,~cachey,~cachez,D,key);
  yes:=ay in HY and -qz(az) in comps[9-D];
  row:=[*D,i,key,ed,Eltseq(ay),Eltseq(az),yes select 1 else 0*];Append(~single,row);Append(~own[D],row);
  if yes then passed+:=1;Append(~survivors,[*1,D,i*]);end if;
  if #own[D] mod 100 eq 0 then print "M2N42_SMALL_SINGLE_PROGRESS",D,#own[D],passed,Cputime(t0);end if;
 end for;
 print "M2N42_SMALL_SINGLE_TYPE",D,#own[D],passed,Cputime(t0);
 Write(OutputFile,Sprint([*groupdata,fields,single,small,pairs,survivors*]):Overwrite:=true);
end for;
if Shard eq 0 then
allsmall:=eval Read("C2x42/data/M2N42_D8_direct_small_all_points.txt");assert #allsmall eq 1 and allsmall[1][1] eq 4;
for i in [1..#allsmall[1][3]] do
 vals:=allsmall[1][3][i];key:=[Integers()!a:a in vals[1]] cat [Integers()!a:a in vals[2]];
 ay:=GY!0;az:=GZ!0;ed:=[];Push(~ay,~az,~ed,~cachey,~cachez,4,key);
 Append(~small,[*4,i,key,ed,Eltseq(ay),Eltseq(az)*]);
end for;
assert #small eq 126;
for i in [1..#own[4]] do for j in [1..#small] do
 ay:=GY!own[4][i][5]+GY!small[j][5];az:=GZ!own[4][i][6]+GZ!small[j][6];
 yes:=ay in HY and az in HZ;Append(~pairs,[*i,j,Eltseq(ay),Eltseq(az),yes select 1 else 0*]);
 if yes then Append(~survivors,[*2,4,i,4,j*]);end if;
end for;end for;
assert #pairs eq 630;end if;
Write(OutputFile,Sprint([*groupdata,fields,single,small,pairs,survivors*]):Overwrite:=true);
print "M2N42_SMALL_QUOTIENT_FILTER_DONE",Shard,#single,#pairs,#survivors,"complete=true",Cputime(t0);quit;
