SetColumns(0);SetQuitOnError(false);SetMemoryLimit(12*10^9);tt54:=Cputime();
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
