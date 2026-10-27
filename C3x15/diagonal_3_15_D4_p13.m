/* Third-prime refinement of the COMPLETE two-prime arithmetic export.
   The12 global predictions for3[D-4c0] are stable under the AFFINE
   action a -> g(a)+12[g(c0)-c0] of all48 modular automorphisms.
   This follows from completeness of the two previous divisor sets,
   equivariant cusp matching, and stability of the global cusp subgroup.
   Thus only the first NC support is reduced by this group here. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/class_rational_reduction.m");
Attach("common/fast_full_aut_hash_m3_15.m");
raw:=eval Read("C3x15/data/diagonal_3_15_D4_arithmetic.txt");
assert #raw[13] eq 5102 and #raw[14] eq 6264 and #raw[15] eq 60;
G7:=AbelianGroup(raw[5]);assert [Eltseq(r):r in Rows(RelationMatrix(G7))] eq raw[6];
Free:=FreeAbelianGroup(20);h7:=hom<Free->G7|[G7!v:v in raw[8]]>;H7:=Image(h7);
assert Invariants(H7) eq [3,3,24,96,96];
indices4:=Sort(Setseq({p[1]:p in raw[15]}));assert #indices4 eq 60;
pred7vectors:=Sort(Setseq({raw[13][i][5]:i in indices4}));assert #pred7vectors eq 12;
for pair in raw[15] do
 assert G7!raw[13][pair[1]][5] eq 3*(G7!raw[14][pair[2]][4]);
end for;
X:=MDX11(3,15,GF(13));C:=Curve(X);FF:=FunctionField(C);
E:=EllipticCurve([FF!f:f in X`_E]);P:=E![FF!f:f in X`_P];Q:=E![FF!f:f in X`_Q];
ai:=aInvariants(E);assert ai[2] eq 0 and ai[4] eq 0 and ai[5] eq 0;
fs:=[ai[3],ai[1]^3-27*ai[3]];
for k in [1..7] do
 for j in [0..2] do
  A:=k*Q+j*P;assert A[3] ne 0 and A[1] ne 0 and A[2] ne 0;
  Append(~fs,A[1]);Append(~fs,A[2]);
 end for;
 A:=k*Q;B:=A+P;assert A ne B and A ne -B;Append(~fs,A[1]-B[1]);
end for;
groups:=AssociativeArray();cs:=Cusps(X);
assert {*Degree(c):c in cs*} eq {*1^^16,4^^4*};
for c in cs do
 sig:=[Valuation(f,c):f in fs];
 if IsDefined(groups,sig) then Append(~groups[sig],c);else groups[sig]:=[c];end if;
end for;
sigs:=Sort(Setseq(Keys(groups)));assert sigs eq raw[1];
blocks:=[&+[Divisor(c):c in groups[sig]]:sig in sigs];
ds:=[Degree(o):o in blocks];assert ds eq raw[2];c0:=Support(blocks[Index(ds,1)])[1];
print "M3N15_DIAGONAL13_CUSPS",ds,Cputime(t0);
G13,m1,m2:=ClassGroupRationalReduction(C,c0:Proof:=true);
J:=JH(45,[16]);fp:=FrobeniusPolynomial(J,13);assert Dimension(J) eq 9;
assert &*[n:n in Invariants(G13)|n ne 0] eq Integers()!Evaluate(fp,1);
h13:=hom<Free->G13|[m2(blocks[i]-ds[i]*Divisor(c0)):i in [1..20]]>;H13:=Image(h13);
assert Invariants(H13) eq Invariants(H7);
for v in Generators(Kernel(h7)) do assert h13(Free!v) eq G13!0;end for;
for v in Generators(Kernel(h13)) do assert h7(Free!v) eq G7!0;end for;
ims:=[h13(Free!((G7!H7.i)@@h7)):i in [1..Ngens(H7)]];
transport:=hom<H7->G13|ims>;
for i in [1..20] do assert transport(H7!h7(Free.i)) eq h13(Free.i);end for;
pred13:=[transport(H7!(G7!v)):v in pred7vectors];S:=AssociativeArray();
for i in [1..#pred13] do
 ky:=Eltseq(pred13[i]);assert not IsDefined(S,ky);S[ky]:=i;
end for;
print "M3N15_DIAGONAL13_GROUP",Invariants(G13),Invariants(H13),#pred13,Cputime(t0);
places:=[NoncuspidalPlaces(X,d):d in [1..4]];assert #places[1] eq 0;
rr:=[* places[1] *];labels:=[* [] *];
for d in [2..4] do
 rs,ll:=PlacesUpToFullAutM3FifteenHash(X,places[d]);Append(~rr,rs);Append(~labels,ll);
 print "M3N15_DIAGONAL13_ORBITS",d,#places[d],#rs,Cputime(t0);
end for;
pc:=[* [G13|] *];
for d in [2..4] do
 cls:=[];
 for i in [1..#rr[d]] do
  Append(~cls,m2(Divisor(rr[d][i])-d*Divisor(c0)));
  if i mod 25 eq 0 then print "M3N15_DIAGONAL13_CLASS_PROGRESS",d,i,#rr[d],Cputime(t0);end if;
 end for;
 Append(~pc,cls);
end for;
rat:=[c:c in cs|Degree(c) eq 1];rcls:=[m2(Divisor(c)-Divisor(c0)):c in rat];
function WeightedVectors(ds,d)
 if #ds eq 0 then return d eq 0 select [[]] else [];end if;
 ans:=[];
 for k in [0..d div ds[1]] do
  for v in WeightedVectors(ds[2..#ds],d-k*ds[1]) do Append(~ans,[k] cat v);end for;
 end for;
 return ans;
end function;
states13:=[* *];matches:=[* *];counts:=[];
for d in [2..4] do
 vs:=WeightedVectors([1:i in [1..16]],4-d);cnt:=0;
 for i in [1..#rr[d]] do
  for j in [1..#vs] do
   v:=vs[j];cl:=pc[d][i]+&+[v[k]*rcls[k]:k in [1..16]];
   Append(~states13,[*d,i,j,Eltseq(cl)*]);cnt+:=1;key:=Eltseq(3*cl);
   if IsDefined(S,key) then Append(~matches,[*#states13,S[key]*]);end if;
  end for;
 end for;
 Append(~counts,cnt);print "M3N15_DIAGONAL13_PATTERN",d,cnt,#matches,Cputime(t0);
end for;
nc2:=[m2(Divisor(x)-2*Divisor(c0)):x in places[2]];cnt:=0;
for i in [1..#rr[2]] do
 for j in [1..#places[2]] do
  cl:=pc[2][i]+nc2[j];Append(~states13,[*22,i,j,Eltseq(cl)*]);cnt+:=1;key:=Eltseq(3*cl);
  if IsDefined(S,key) then Append(~matches,[*#states13,S[key]*]);end if;
 end for;
end for;
Append(~counts,cnt);
Write("C3x15/data/diagonal_3_15_D4_p13_arithmetic.txt",Sprint([*
 sigs,ds,Invariants(G7),[Eltseq(v):v in Rows(RelationMatrix(G7))],
 Invariants(G13),[Eltseq(v):v in Rows(RelationMatrix(G13))],
 [Eltseq(h7(Free.i)):i in [1..20]],[Eltseq(h13(Free.i)):i in [1..20]],
 pred7vectors,[Eltseq(a):a in pred13],states13,matches,counts,
 [#ss:ss in places],[#ss:ss in rr],labels *]):Overwrite:=true);
printf "M3N15_DIAGONAL13_DONE states=%o patterncounts=%o predictions=%o matches=%o complete=true time=%o\n",#states13,counts,#pred13,#matches,Cputime(t0);quit;
