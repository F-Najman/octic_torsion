SetColumns(0);SetQuitOnError(true);SetMemoryLimit(12*10^9);tt54:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/class_rational_reduction.m");
function CuspDivisors54(cs,d,C)
 if d eq 0 then return [DivisorGroup(C)!0];end if;
 if #cs eq 0 then return [];end if;
 P:=cs[1];rest:=cs[2..#cs];out:=[];
 for k in [0..d div Degree(P)] do
  out cat:=[k*Divisor(P)+a:a in CuspDivisors54(rest,d-k*Degree(P),C)];
 end for;
 return out;
end function;
X2:=MDX1(27,GF(2));C2:=Curve(X2);assert Genus(C2) eq 13;
os2:=[Type(o) eq PlcCrvElt select Divisor(o) else o:o in CuspOrbitsQ(X2)];
sig2:=[CuspSignature(X2,Support(o)[1]):o in os2];assert #Seqset(sig2) eq #os2;
c02:=Support([o:o in os2|Degree(o) eq 1][1])[1];baseSig:=CuspSignature(X2,c02);
G2,a2,b2:=ClassGroupRationalReduction(C2,c02:Proof:=true);qc2:=[b2(o-Degree(o)*Divisor(c02)):o in os2];
H2:=sub<G2|qc2>;assert Invariants(G2) eq [3,3,52497,0] and Invariants(H2) eq [3,3,52497];
print "OCTIC54_F2_CLASSES",Invariants(G2),Invariants(H2),Cputime(tt54);
ps2:=[];
for d in [1..8] do
 Append(~ps2,NoncuspidalPlaces(X2,d));print "OCTIC54_F2_NC_DEGREE",d,#ps2[d],Cputime(tt54);
end for;
assert [#v:v in ps2[1..7]] eq [0,0,0,0,0,9,9];
candidate2:={G2|};coverage2:=[];
for d in [6..8] do
 for placeIndex in [1..#ps2[d]] do
  p:=ps2[d][placeIndex];cs:=CuspDivisors54(Cusps(X2),8-d,C2);
  for cc in cs do
   v:=b2(Divisor(p)+cc-8*Divisor(c02));assert v in H2;Include(~candidate2,v);
  end for;
 end for;
 Append(~coverage2,<d,#ps2[d],#CuspDivisors54(Cusps(X2),8-d,C2)>);
end for;
print "OCTIC54_F2_CANDIDATES",#candidate2,coverage2,Cputime(tt54);
X5:=MDX1(27,GF(5));C5:=Curve(X5);assert Genus(C5) eq 13;
raw5:=[Type(o) eq PlcCrvElt select Divisor(o) else o:o in CuspOrbitsQ(X5)];sig5:=[CuspSignature(X5,Support(o)[1]):o in raw5];
assert Seqset(sig5) eq Seqset(sig2) and #Seqset(sig5) eq #raw5;os5:=[raw5[Index(sig5,s)]:s in sig2];
c05:=Support(os5[Index(sig2,baseSig)])[1];assert Degree(c05) eq 1;
G5,a5,b5:=ClassGroupRationalReduction(C5,c05:Proof:=true);qc5:=[b5(o-Degree(o)*Divisor(c05)):o in os5];H5:=sub<G5|qc5>;
assert Invariants(G5) eq [3,3,342175446,0] and Invariants(H5) eq [3,3,52497];
L:=FreeAbelianGroup(#os2);h2:=hom<L->G2|qc2>;h5:=hom<L->G5|qc5>;assert Kernel(h2) eq Kernel(h5);
print "OCTIC54_F5_CLASSES",Invariants(G5),Invariants(H5),"cuspidal_relations_equal=true",Cputime(tt54);
comps:=AssociativeArray();
for k in [0..5] do
 comps[k]:={b5(c-k*Divisor(c05)):c in CuspDivisors54(Cusps(X5),k,C5)};
 print "OCTIC54_F5_CUSP_CLASSES",k,#comps[k],Cputime(tt54);
end for;
print "OCTIC54_CROSSPRIME_PREFLIGHT_DONE",Cputime(tt54);

Attach("common/fast_hecke.m");
Attach("common/fast_hecke_11.m");
Attach("common/fast_hecke_cyclic_commonzero.m");
tl54:=Cputime();allps54:=[[]:d in [1..6]];good54:=[[]:d in [1..6]];inert54:=[];root54:=[];
for d in [2..6] do
 allps54[d]:=NoncuspidalPlaces(X5,d);nrlist:=[];
 for P in allps54[d] do
  EP:=EllipticCurve(X5,P);b2e,b4e,b6e,b8e:=Explode(bInvariants(EP));RE<T>:=PolynomialRing(BaseRing(EP));
  nr:=#Roots(4*T^3+b2e*T^2+2*b4e*T+b6e);assert nr in [0,1,3];Append(~nrlist,nr);
  if nr gt 0 then Append(~good54[d],P);elif d eq 2 then Append(~inert54,P);end if;
 end for;
 Append(~root54,nrlist);print "OCTIC54_LOWER_POINTS",d,#allps54[d],#good54[d],Cputime(tl54);
end for;
assert #good54[2] eq 0 and #inert54 eq 9;
checked54:=0;members54:=0;passing54:=0;hits54:=[];
if assigned TestDataFile54 then testData54:=TestDataFile54;else testData54:="C54/data/C54_lower_classes.txt";end if;
testout54:=Open(testData54,"w");delete testout54; // start a fresh export (Write appends)
procedure TestClass54(a,e,desc,~checked,~members,~passing,~hits)
 for cc in comps[8-e] do
  q:=a+cc;checked+:=1;member:=q in H5;passes:=false;
  if member then
   members+:=1;v:=q@@h5;assert h5(v) eq q;
   if h2(v) in candidate2 then passes:=true;passing+:=1;Append(~hits,<desc,Eltseq(v)>);end if;
  end if;
  Write(testData54,Sprint([Integers()!z:z in Eltseq(q)] cat [member select 1 else 0,passes select 1 else 0]));
 end for;
end procedure;
for d in [3..6] do
 reps:=PlacesUpToDiamondFast(X5,good54[d]);print "OCTIC54_LOWER_SINGLE_START",d,#reps,Cputime(tl54);
 for i in [1..#reps] do TestClass54(b5(Divisor(reps[i])-d*Divisor(c05)),d,<d,i>,~checked54,~members54,~passing54,~hits54);end for;
 print "OCTIC54_LOWER_SINGLE_DONE",d,checked54,members54,passing54,Cputime(tl54);
end for;
for P in PlacesUpToDiamondFast(X5,inert54) do TestClass54(b5(3*Divisor(P)-6*Divisor(c05)),6,<2,3>,~checked54,~members54,~passing54,~hits54);end for;
lowerQ54:=AssociativeArray();for d in [3,4,5] do lowerQ54[d]:=[b5(Divisor(P)-d*Divisor(c05)):P in good54[d]];end for;
for pair in [<3,3>,<3,4>,<3,5>,<4,4>] do
 d,e:=Explode(pair);reps:=PlacesUpToDiamondFast(X5,good54[d]);
 for i in [1..#reps] do
  a:=b5(Divisor(reps[i])-d*Divisor(c05));
  for b in lowerQ54[e] do TestClass54(a+b,d+e,<d,e,i>,~checked54,~members54,~passing54,~hits54);end for;
 end for;
 print "OCTIC54_LOWER_MIXED_DONE",pair,checked54,members54,passing54,Cputime(tl54);
end for;
printf "OCTIC54_LOWER_SUPPORTS_DONE checked=%o members=%o passing=%o hits=%o complete=true time=%o\n",checked54,members54,passing54,hits54,Cputime(tl54);

Attach("C54/C54_coordinate_place.m");
tv54:=Cputime();data7:=eval Read("C54/data/F5_degree7_pointdata.txt");def7,records7:=Explode(data7);
assert #records7 eq 828;R7:=PolynomialRing(GF(5));f7:=R7!def7;assert IsIrreducible(f7) and Degree(f7) eq 7;
K7<a7>:=ext<GF(5)|f7>;FF5:=FunctionField(C5);good7fresh54:=[];weighted7fresh54:=0;jkeys7fresh54:={};
for r in records7 do
 xycoeff,n,r2,jcoeff:=Explode(r);assert n mod 54 eq 0 and r2 in [1,3];
 xy:=[Evaluate(R7!co,a7):co in xycoeff];jj:=Evaluate(R7!jcoeff,a7);assert Degree(MinimalPolynomial(jj,GF(5))) eq 7;
 jkey:=Eltseq(MinimalPolynomial(jj,GF(5)));assert jkey notin jkeys7fresh54;Include(~jkeys7fresh54,jkey);
 P:=MDPlaceCoordinateFullDegree([FF5.1,FF5.2],xy);EP:=EllipticCurve(X5,P);PP:=LevelStructure(X5,P)`P;
 assert 27*PP eq EP!0 and 9*PP ne EP!0;
 b2e,b4e,b6e,b8e:=Explode(bInvariants(EP));RR<TT>:=PolynomialRing(BaseRing(EP));assert #Roots(4*TT^3+b2e*TT^2+2*b4e*TT+b6e) eq r2;
 assert MinimalPolynomial(jInvariant(EP),GF(5)) eq MinimalPolynomial(jj,GF(5));
 Append(~good7fresh54,P);weighted7fresh54+:=9*r2;
end for;
assert #good7fresh54 eq 828 and #Seqset(good7fresh54) eq 828 and #jkeys7fresh54 eq 828 and weighted7fresh54 eq 11412;
for i in [1..#good7fresh54] do
 TestClass54(b5(Divisor(good7fresh54[i])-7*Divisor(c05)),7,<7,i>,~checked54,~members54,~passing54,~hits54);
end for;
printf "OCTIC54_DEGREE7_REVALIDATED representatives=828 weighted_source_places=11412 checked=%o members=%o passing=%o complete=true time=%o\n",checked54,members54,passing54,Cputime(tv54);

assert #candidate2 eq 531 and Kernel(h2) eq Kernel(h5);
if assigned BaseDataFile54 then data54out:=BaseDataFile54;else data54out:="C54/data/C54_class_arithmetic_base.py";end if;
Write(data54out,Sprintf("relations2 = %o\nrelations5 = %o\ncusp2 = %o\ncusp5 = %o\ncandidate2 = %o\n",[Eltseq(r):r in Rows(RelationMatrix(G2))],[Eltseq(r):r in Rows(RelationMatrix(G5))],[Eltseq(v):v in qc2],[Eltseq(v):v in qc5],[Eltseq(v):v in candidate2]):Overwrite:=true);
rawtests54:=Split(Read(testData54),"\n");fulltests54:=[* *];
for line in rawtests54 do
 if #line eq 0 then continue;end if;
 row:=eval line;assert #row eq 6;q:=G5!row[1..4];member:=q in H5;assert (member select 1 else 0) eq row[5];im2:=[];preimage:=[];passes:=false;
 if member then
  v:=q@@h5;assert h5(v) eq q;im2:=Eltseq(h2(v));preimage:=Eltseq(v);passes:=h2(v) in candidate2;
 end if;
 assert (passes select 1 else 0) eq row[6];Append(~fulltests54,[*row[1..4],row[5],row[6],im2,preimage*]);
end for;
if not assigned FullDataFile54 then FullDataFile54:="C54/data/C54_replayed_testclasses.py";end if;classout54:=Open(FullDataFile54,"w");fprintf classout54,"tests = [\n";
for row in fulltests54 do fprintf classout54,"[%o,%o,%o,%o,%o],\n",row[1],row[2],row[3],row[4],row[5];end for;
fprintf classout54,"]\ncomplete = True\n";delete classout54;
printf "OCTIC54_CLASS_ARITHMETIC_EXPORTED tests=%o members=%o\n",#fulltests54,#[r:r in fulltests54|r[2] eq 1];

print "OCTIC54_QUARTIC_SPLITTING",#[r:r in root54[3]|r eq 1],#[r:r in root54[3]|r eq 3];
quit;
