/* Fourth-prime refinement of the COMPLETE two-prime arithmetic export.
   The12 global predictions for3[D-4c0] are stable under the AFFINE
   action a -> g(a)+12[g(c0)-c0] of all48 modular automorphisms.
   This follows from completeness of the two initial complete divisor sets,
   equivariant cusp matching, and stability of the global cusp subgroup.
   Thus only the first NC support is reduced by this group here. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("common/class_rational_reduction.m");
Attach("common/fast_full_aut_hash_m3_15.m");
Attach("C3x15/fast_full_aut_hash_m3_15_native.m");
Attach("C3x15/fast_full_aut_hash_m3_15_keytable.m");
Attach("common/fast_hecke.m");
Attach("common/fast_hecke_11.m");
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
X:=MDX11(3,15,GF(19));C:=Curve(X);FF:=FunctionField(C);
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
assert {*Degree(c):c in cs*} eq {*1^^16,2^^8*};
for c in cs do
 sig:=[Valuation(f,c):f in fs];
 if IsDefined(groups,sig) then Append(~groups[sig],c);else groups[sig]:=[c];end if;
end for;
sigs:=Sort(Setseq(Keys(groups)));assert sigs eq raw[1];
blocks:=[&+[Divisor(c):c in groups[sig]]:sig in sigs];
ds:=[Degree(o):o in blocks];assert ds eq raw[2];c0:=Support(blocks[Index(ds,1)])[1];
print "M3N15_DIAGONAL19_CUSPS",ds,Cputime(t0);
AFg,mg:=AlgorithmicFunctionField(FF);
print "M3N15_DIAGONAL19_DEGREE3_CLASSGROUP_START",Cputime(t0);
G19,_,phi19:=ClassGroup(AFg:ReductionDivisor:=FunctionFieldDivisor(Divisor(c0)),DegreeBound:=3,Proof:=true);
m2:=map<DivisorGroup(C)->G19|D:->phi19(FunctionFieldDivisor(D)),g:->CurveDivisor(C,g@@phi19)>;
print "M3N15_DIAGONAL19_DEGREE3_CLASSGROUP_DONE",Invariants(G19),Cputime(t0);
J:=JH(45,[16]);fp:=FrobeniusPolynomial(J,19);assert Dimension(J) eq 9;
assert &*[n:n in Invariants(G19)|n ne 0] eq Integers()!Evaluate(fp,1);
h19:=hom<Free->G19|[m2(blocks[i]-ds[i]*Divisor(c0)):i in [1..20]]>;H19:=Image(h19);
assert Invariants(H19) eq Invariants(H7);
for v in Generators(Kernel(h7)) do assert h19(Free!v) eq G19!0;end for;
for v in Generators(Kernel(h19)) do assert h7(Free!v) eq G7!0;end for;
ims:=[h19(Free!((G7!H7.i)@@h7)):i in [1..Ngens(H7)]];
transport:=hom<H7->G19|ims>;
for i in [1..20] do assert transport(H7!h7(Free.i)) eq h19(Free.i);end for;
pred19:=[transport(H7!(G7!v)):v in pred7vectors];S:=AssociativeArray();
for i in [1..#pred19] do
 ky:=Eltseq(pred19[i]);assert not IsDefined(S,ky);S[ky]:=i;
end for;
print "M3N15_DIAGONAL19_GROUP",Invariants(G19),Invariants(H19),#pred19,Cputime(t0);
AF,mma:=AlgorithmicFunctionField(FF);afcs:=[FunctionFieldPlace(c):c in cs];
afcoords:=[mma(f):f in X`_coordinates];
places:=[* [x:x in Places(AF,d)|not x in afcs]:d in [1..3] *];
inputkeys:=[];
for shard in [1..6] do
 tab:=eval Read("C3x15/data/M3N15_p19_native_keys_shard" cat IntegerToString(shard) cat ".txt");
 assert tab[1] eq shard and tab[2] eq 1+(shard-1)*5500 and tab[3] eq Minimum(shard*5500,32076);
 assert tab[4] eq 32076 and tab[5] eq Coefficients(fp);
 assert #tab[6] eq tab[3]-tab[2]+1 and [r[1]:r in tab[6]] eq [tab[2]..tab[3]];
 inputkeys cat:= [r[2]:r in tab[6]];
end for;
assert #inputkeys eq 32076 and #Seqset(inputkeys) eq 32076;
Append(~places,[1..#inputkeys]);assert #places[1] eq 0;
traces:=[];closed:=[];
for n in [1..4] do
 sn:=-n*Coefficient(fp,18-n);
 for j in [1..n-1] do sn-:=Coefficient(fp,18-j)*traces[n-j];end for;
 Append(~traces,Integers()!sn);count:=19^n+1-Integers()!sn;
 for e in Divisors(n) do if e lt n then count-:=e*closed[e];end if;end for;
 assert count mod n eq 0;Append(~closed,count div n);
 assert closed[n] eq #places[n]+#[c:c in cs|Degree(c) eq n];
end for;
print "M3N15_DIAGONAL19_FROBENIUS_COUNTS",closed,[#ss:ss in places],Cputime(t0);
rr:=[* [] *];labels:=[* [] *];
for d in [2..4] do
 if d lt 4 then
  rsAF,ll:=PlacesUpToFullAutM3FifteenNative(afcoords,places[d],19);
  rs:=[CurvePlace(C,x):x in rsAF];
 else
  repkeys,ll:=PartitionFullAutM3FifteenKeyTable(inputkeys,19,4);K:=GF(19^4);rs:=[];
  for ky in repkeys do
   pair:=M3FifteenPairFromNativeKey(ky,K,19,4);
   y:=MDPlaceCommonZeros(X`_coordinates,pair);assert Degree(y) eq 4;
   assert M3FifteenNativeKey([Evaluate(f,y):f in X`_coordinates],19,4) eq ky;
   Append(~rs,y);
  end for;
  // Compare the first full orbit against the original inverse explicitly.
  x:=rs[1];EE:=EllipticCurve(X,x);Ls:=LevelStructure(X,x);PP:=Ls`P;QQ:=Ls`Q;
  Q3:=5*QQ;Q5:=3*QQ;actual:={};
  for mat in M3FifteenAutMatrices() do
   a,b,c,d1:=Explode(mat);
   for u in [1,2] do
    La:=rec<recformat<P:PtEll,Q:PtEll>|P:=a*PP+b*Q3,Q:=2*(c*PP+d1*Q3)+2*u*Q5>;
    y:=ModuliPointFast(X,EE,La);assert Degree(y) eq 4;
    Include(~actual,M3FifteenNativeKey([Evaluate(f,y):f in X`_coordinates],19,4));
   end for;
  end for;
  assert actual eq {inputkeys[i]:i in [1..#inputkeys]|ll[i] eq 1};
  print "M3N15_DIAGONAL19_KEYTABLE_REFERENCE_EQUAL",#actual,Cputime(t0);
 end if;
 Append(~rr,rs);Append(~labels,ll);
 print "M3N15_DIAGONAL19_ORBITS",d,#places[d],#rs,Cputime(t0);
end for;
pc:=[* [G19|] *];
for d in [2..4] do
 cls:=[];
 for i in [1..#rr[d]] do
  print "M3N15_DIAGONAL19_CLASS_BEFORE",d,i,#rr[d],Cputime(t0);
  Append(~cls,m2(Divisor(rr[d][i])-d*Divisor(c0)));
  print "M3N15_DIAGONAL19_CLASS_AFTER",d,i,#rr[d],Cputime(t0);
  if i mod 25 eq 0 then print "M3N15_DIAGONAL19_CLASS_PROGRESS",d,i,#rr[d],Cputime(t0);end if;
 end for;
 Append(~pc,cls);
end for;
rcls:=[m2(Divisor(c)-Degree(c)*Divisor(c0)):c in cs];
assert #rcls eq 24;
// Degree-two complements include all eight quadratic cusps, as well
// as the136 unordered rational-cusp pairs. No cusp complement is omitted.
function WeightedVectors(ds,d)
 if #ds eq 0 then return d eq 0 select [[]] else [];end if;
 ans:=[];
 for k in [0..d div ds[1]] do
  for v in WeightedVectors(ds[2..#ds],d-k*ds[1]) do Append(~ans,[k] cat v);end for;
 end for;
 return ans;
end function;
states19:=[* *];matches:=[* *];counts:=[];
for d in [2..4] do
 vs:=WeightedVectors([Degree(c):c in cs],4-d);cnt:=0;
 assert #vs eq [144,16,1][d-1];
 for i in [1..#rr[d]] do
  for j in [1..#vs] do
   v:=vs[j];cl:=pc[d][i]+&+[v[k]*rcls[k]:k in [1..#cs]];
   Append(~states19,[*d,i,j,Eltseq(cl)*]);cnt+:=1;key:=Eltseq(3*cl);
   if IsDefined(S,key) then Append(~matches,[*#states19,S[key]*]);end if;
  end for;
 end for;
 Append(~counts,cnt);print "M3N15_DIAGONAL19_PATTERN",d,cnt,#matches,Cputime(t0);
end for;
nc2:=[m2(Divisor(CurvePlace(C,x))-2*Divisor(c0)):x in places[2]];cnt:=0;
for i in [1..#rr[2]] do
 for j in [1..#places[2]] do
  cl:=pc[2][i]+nc2[j];Append(~states19,[*22,i,j,Eltseq(cl)*]);cnt+:=1;key:=Eltseq(3*cl);
  if IsDefined(S,key) then Append(~matches,[*#states19,S[key]*]);end if;
 end for;
end for;
Append(~counts,cnt);
Write("C3x15/data/diagonal_3_15_D4_p19_arithmetic.txt",Sprint([*
 sigs,ds,Invariants(G7),[Eltseq(v):v in Rows(RelationMatrix(G7))],
 Invariants(G19),[Eltseq(v):v in Rows(RelationMatrix(G19))],
 [Eltseq(h7(Free.i)):i in [1..20]],[Eltseq(h19(Free.i)):i in [1..20]],
 pred7vectors,[Eltseq(a):a in pred19],states19,matches,counts,
 [#ss:ss in places],[#ss:ss in rr],labels *]):Overwrite:=true);
printf "M3N15_DIAGONAL19_DONE states=%o patterncounts=%o predictions=%o matches=%o complete=true time=%o\n",#states19,counts,#pred19,#matches,Cputime(t0);quit;
