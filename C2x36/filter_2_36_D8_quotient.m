/* Exact quotient necessary condition. NO global exclusion assumed. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/fast_hecke.m");Attach("common/fast_hecke_11.m");
Attach("common/class_rational_reduction.m");
Y:=MDX11(2,18,GF(5));C:=Curve(Y);cs:=Cusps(Y);c0:=Rep({c:c in cs|Degree(c) eq 1});
G,m1,m2:=ClassGroupRationalReduction(C,c0:Proof:=true);hg:=[m2(Divisor(c)-Degree(c)*Divisor(c0)):c in cs];H:=sub<G|hg>;
assert Invariants(G) eq [6,84,252,0] and Invariants(H) eq [2,42,126];
groupdata:=[*Invariants(G),[Eltseq(r):r in Rows(RelationMatrix(G))],[Degree(c):c in cs],[Eltseq(v):v in hg],Invariants(H)*];
function BcKey(bc,D)
 keys:=[];for i in [1..D] do Append(~keys,&cat[[Integers()!x:x in Eltseq(z,GF(5))]:z in bc]);bc:=[z^5:z in bc];end for;return Sort(keys)[1];
end function;
procedure PushClass(~cache,~answer,~edata,b,c,D,Y,c0,m2)
 dd:=b*(1-b-c);s:=(3*b^2+2*c*b+dd)/(2*b);xx:=s^2-c-2*b;yy:=-b+s*(b-xx);assert xx ne 0 and yy ne 0;
 bc:=[xx^3/yy^2,c*xx^2/yy^2];ky:=[D] cat BcKey(bc,D);
 if IsDefined(cache,ky) then answer:=cache[ky][1];edata:=cache[ky][2];return;end if;
 y:=MDPlaceCommonZeros(Y`_coordinates,bc);e:=Degree(y);assert e ge 3 and D mod e eq 0 and D div e in [1,2];
 w:=D div e;answer:=w*m2(Divisor(y)-e*Divisor(c0));edata:=[e,w];cache[ky]:=<answer,edata>;
end procedure;
R:=PolynomialRing(GF(5));fields:=[* *];single:=[* *];small:=[* *];pairs:=[* *];survivors:=[* *];own:=AssociativeArray();cache:=AssociativeArray();
for D in [3..8] do
 a:=eval Read(Sprintf("C2x36/data/M2N36_D8_direct_full%o_arithmetic.txt",D));
 if D eq 8 then def:=a[1];records:=a[8];else assert a[1] eq D;def:=a[2];records:=a[9];end if;
 K<z>:=ext<GF(5)|R!def>;Append(~fields,[*D,def*]);own[D]:=[* *];passed:=0;
 for i in [1..#records] do
  d,key:=Explode(records[i]);if d ne D then continue;end if;
  b:=Evaluate(R!key[1..D],z);c:=Evaluate(R!key[D+1..2*D],z);cl:=G!0;ed:=[];
  PushClass(~cache,~cl,~ed,b,c,D,Y,c0,m2);yes:=cl in H;row:=[*D,i,key,ed[1],ed[2],Eltseq(cl),yes select 1 else 0*];Append(~single,row);Append(~own[D],row);
  if yes then passed+:=1;Append(~survivors,[*1,D,i*]);end if;
  if #own[D] mod 100 eq 0 then print "M2N36_D8_QUOTIENT_SINGLE_PROGRESS",D,#own[D],passed,Cputime(t0);end if;
 end for;
 print "M2N36_D8_QUOTIENT_SINGLE_TYPE",D,#own[D],passed,Cputime(t0);
 Write("C2x36/data/filter_2_36_D8_quotient_arithmetic.txt",Sprint([*groupdata,fields,single,small,pairs,survivors*]):Overwrite:=true);
end for;
assert #single eq 2712;
allsmall:=eval Read("C2x36/data/M2N36_D8_direct_small_all_points.txt");smalls:=AssociativeArray();
for data in allsmall do
 D:=data[1];def:=[Integers()!x:x in data[2]];K<z>:=ext<GF(5)|R!def>;smalls[D]:=[* *];
 for i in [1..#data[3]] do
  vals:=data[3][i];b:=Evaluate(R![Integers()!x:x in vals[1]],z);c:=Evaluate(R![Integers()!x:x in vals[2]],z);cl:=G!0;ed:=[];
  PushClass(~cache,~cl,~ed,b,c,D,Y,c0,m2);key:=[Integers()!x:x in vals[1]] cat [Integers()!x:x in vals[2]];
  row:=[*D,i,key,ed[1],ed[2],Eltseq(cl)*];Append(~small,row);Append(~smalls[D],row);
 end for;
 print "M2N36_D8_QUOTIENT_ALL_SMALL",D,#smalls[D],Cputime(t0);
end for;
assert #small eq 888;
for typ in [[3,3],[3,4],[3,5],[4,4]] do
 D,E:=Explode(typ);passed:=0;tested:=0;
 for i in [1..#own[D]] do
  cl1:=G!own[D][i][6];
  for j in [1..#smalls[E]] do
   cl:=cl1+G!smalls[E][j][6];yes:=cl in H;tested+:=1;Append(~pairs,[*D,i,E,j,Eltseq(cl),yes select 1 else 0*]);
   if yes then passed+:=1;Append(~survivors,[*2,D,i,E,j*]);end if;
  end for;
 end for;
 print "M2N36_D8_QUOTIENT_PAIR_TYPE",D,E,tested,passed,Cputime(t0);
end for;
assert #pairs eq 2376;
Write("C2x36/data/filter_2_36_D8_quotient_arithmetic.txt",Sprint([*groupdata,fields,single,small,pairs,survivors*]):Overwrite:=true);
print "M2N36_D8_QUOTIENT_FILTER_DONE",#single,#pairs,#survivors,"complete=true",Cputime(t0);quit;
